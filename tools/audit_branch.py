#!/usr/bin/env python3
"""Read-only structural/privacy check for the combined KSSOLV/WI source package."""
import hashlib
import json
import re
import sys
from pathlib import Path
from urllib.parse import unquote
ROOT=Path(__file__).resolve().parents[1]
SCRATCH={'_validation','_work','_output','__pycache__'}
AUTHOR_PATH=re.compile(r'(?i)(?:(?<![A-Za-z0-9%])[A-Z]:[\\/](?!Program Files)|/(?:public/home|home|Users)/[A-Za-z0-9_.-]+)')
LINK=re.compile(r'(?<!!)\[[^\]\n]+\]\(([^\s)]+)(?:\s+"[^"]*")?\)')
SECRETS=[re.compile(r'-{5}BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-{5}'),re.compile(r'\bgh[pousr]_[A-Za-z0-9]{30,}\b'),re.compile(r'\bgithub_pat_[A-Za-z0-9_]{30,}\b'),re.compile(r'\bsk-(?:proj-)?[A-Za-z0-9_-]{30,}\b'),re.compile(r'\bAKIA[A-Z0-9]{16}\b')]
def audit():
    errors=[]
    manifest=json.loads((ROOT/'docs/wi/SOURCE_MANIFEST.json').read_text(encoding='utf-8'))['files']
    expected={e['destination']:e for e in manifest}
    for rel,e in expected.items():
        f=(ROOT/rel).resolve()
        if ROOT not in f.parents or not f.is_file():errors.append('Missing/unsafe copied file: '+rel);continue
        if hashlib.sha256(f.read_bytes()).hexdigest()!=e['release_sha256']:errors.append('Hash mismatch: '+rel)
    files=[];scratch_count=0
    for f in ROOT.rglob('*'):
        if not f.is_file():continue
        parts=f.relative_to(ROOT).parts
        if '.git' in parts:continue
        if any(p in SCRATCH for p in parts):scratch_count+=1;continue
        files.append(f);rel=f.relative_to(ROOT).as_posix()
        if f.stat().st_size>10*1024*1024:errors.append('Unexpected file larger than 10 MiB: '+rel)
        if f.suffix.lower() in {'.mat','.log','.asv','.swp','.pyc','.amn','.mmn','.eig','.cube','.chk'} or f.name.startswith('UNK') or '.bak_' in f.name:
            errors.append('Generated or backup file: '+rel)
        if f.suffix.lower() in {'.pem','.key','.ppk'} or f.name.startswith('.env') or 'askpass' in f.name.lower():errors.append('Credential/local filename: '+rel)
        raw=f.read_bytes()
        # Binary assets are allowed only when inherited and explicitly manifested.
        try:text=raw.decode('utf-8-sig')
        except UnicodeDecodeError:
            if rel not in expected:errors.append('Unmanifested binary: '+rel)
            if re.search(rb'/(?:public/home|home|Users)/[A-Za-z0-9_.-]+',raw):errors.append('Embedded author path in binary: '+rel)
            continue
        if AUTHOR_PATH.search(text):errors.append('Author-specific absolute path: '+rel)
        if any(pattern.search(text) for pattern in SECRETS):errors.append('Possible credential: '+rel)
        if f.suffix.lower()=='.m':
            if re.search(r'\bsetup_wi\s*(?:\(|;)',text):errors.append('Standalone setup call remains: '+rel)
            if re.search(r"\b(?:getenv|setenv)\s*\(\s*['\"](?:KSSOLV_ROOT|WI_ROOT)['\"]",text):errors.append('External-root environment dependency remains: '+rel)
        if f.suffix.lower()=='.md':
            for m in LINK.finditer(text):
                target=m.group(1).strip('<>')
                if re.match(r'^[A-Za-z][A-Za-z0-9+.-]*:',target) or target.startswith('#'):continue
                d=(f.parent/unquote(target.split('#',1)[0])).resolve()
                if ROOT not in d.parents or not d.exists():errors.append('Broken/outside Markdown link: '+rel+' -> '+target)
    if scratch_count:errors.append('Runtime/scratch files remain: '+str(scratch_count))
    if (ROOT/'setup_wi.m').exists():errors.append('Standalone setup file unexpectedly included.')
    wi=list((ROOT/'src/WI').glob('*.m'))
    if len(wi)!=27:errors.append('Expected 27 WI numerical source files.')
    for f in wi:
        e=expected.get(f.relative_to(ROOT).as_posix(),{})
        if not e.get('identical_to_source'):errors.append('WI numerical source was modified: '+f.name)
    validation=json.loads((ROOT/'docs/wi/VALIDATION_RESULTS.json').read_text(encoding='utf-8'))
    for rel in validation['static_dependency_closure']['files']:
        if not (ROOT/rel).is_file():errors.append('Missing recorded static dependency: '+rel)
    return {'status':'PASS' if not errors else 'FAIL','files':len(files),'bytes':sum(f.stat().st_size for f in files),'copied_files_checked':len(manifest),'wi_core_files':len(wi),'errors':errors,'publication_rights_cleared':False,'note':'Structural PASS does not settle source licenses, PP redistribution or complete material/MLWF validation.'}
if __name__=='__main__':
    report=audit();print(json.dumps(report,ensure_ascii=False,indent=2));sys.exit(0 if report['status']=='PASS' else 1)
