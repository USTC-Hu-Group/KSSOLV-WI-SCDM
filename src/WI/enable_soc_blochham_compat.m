function info = enable_soc_blochham_compat
%ENABLE_SOC_BLOCHHAM_COMPAT Fix the known SOC branch in a session-local copy.
% Call after KSSOLV_startup, BEFORE creating/loading any BlochHam objects.
% No KSSOLV source, saved MATLAB path or checkpoint is modified. The copy
% preserves the class name/properties so existing checkpoints remain usable.
persistent cached source_snapshot;
pptype_file = which('kssolvpptype');
if isempty(pptype_file)
    error('enable_soc_blochham_compat:MissingKSSOLV', 'Run KSSOLV_startup first.');
end
root = fileparts(fileparts(pptype_file));
source_dir = fullfile(root, 'src', 'BasisClass', '@BlochHam');
source_file = fullfile(source_dir, 'BlochHam.m');
if ~isfile(source_file) || isfile(fullfile(source_dir, 'BlochHam.p'))
    error('enable_soc_blochham_compat:UnsupportedClass', ...
        'Expected BlochHam.m without a shadowing .p file at %s.', source_file);
end
source = fileread(source_file);
if ~isempty(cached)
    if ~strcmp(cached.source_file, source_file) || ~strcmp(source_snapshot, source)
        error('enable_soc_blochham_compat:SourceChanged', ...
            'KSSOLV changed during this MATLAB session. Start a new MATLAB process.');
    end
    if cached.runtime_copy
        addpath(cached.path_root, '-begin');
    end
    assert_resolved(cached.class_file);
    info = cached;
    return;
end

% Only accept the known guards (or their already-fixed equivalents).
fixed = rewrite_guard(source, ...
    '(?m)^([ \t]*)if[ \t]+(cry\.nspin[^\r\n]*)(\r?\n[ \t]+BH\.vnp[ \t]*=[ \t]*vhart[ \t]*\+[ \t]*vxc;)', ...
    'cry.nspin==1', 'cry.nspin==1||(cry.nspin==4&&~cry.domag)', ...
    'if cry.nspin == 1 || (cry.nspin == 4 && ~cry.domag)');
fixed = rewrite_guard(fixed, ...
    '(?m)^([ \t]*)elseif[ \t]+(cry\.nspin[ \t]*==[ \t]*4[^\r\n]*)(\r?\n[ \t]+BH\.vnp[ \t]*=[ \t]*cell\(cry\.nspin,1\);)', ...
    'cry.nspin==4', 'cry.nspin==4&&cry.domag', ...
    'elseif cry.nspin == 4 && cry.domag');
info = struct('source_file', source_file, 'class_file', source_file, ...
    'runtime_copy', false, 'path_root', '');
if ~strcmp(fixed, source)
    [~, ~, loaded_classes] = inmem;
    if any(strcmp(loaded_classes, 'BlochHam'))
        error('enable_soc_blochham_compat:ClassAlreadyLoaded', ...
            ['BlochHam is already loaded. Start a new MATLAB process and call ', ...
             'enable_soc_blochham_compat BEFORE loading checkpoints or constructing Hamiltonians.']);
    end
    info.path_root = tempname;
    mkdir(info.path_root);
    destination = fullfile(info.path_root, '@BlochHam');
    [ok, message] = copyfile(source_dir, destination);
    if ~ok
        error('enable_soc_blochham_compat:CopyFailed', '%s', message);
    end
    info.class_file = fullfile(destination, 'BlochHam.m');
    fid = fopen(info.class_file, 'w');
    if fid < 0
        error('enable_soc_blochham_compat:WriteFailed', 'Cannot write %s.', info.class_file);
    end
    cleanup = onCleanup(@() fclose(fid));
    count = fprintf(fid, '%s', fixed);
    if count < numel(fixed)
        error('enable_soc_blochham_compat:WriteFailed', 'Incomplete compatibility copy.');
    end
    clear cleanup;
    info.runtime_copy = true;
    addpath(info.path_root, '-begin');
    rehash;
end
assert_resolved(info.class_file);
if ~strcmp(fileread(source_file), source)
    error('enable_soc_blochham_compat:SourceChanged', 'KSSOLV source changed during activation.');
end
fprintf('SOC compatibility source (unchanged): %s\n', source_file);
fprintf('SOC compatibility active class: %s\n', info.class_file);
cached = info;
source_snapshot = source;
end

function output = rewrite_guard(source, pattern, old_guard, new_guard, replacement)
[first, last, tokens] = regexp(source, pattern, 'start', 'end', 'tokens');
if numel(first) ~= 1
    error('enable_soc_blochham_compat:UnknownVersion', ...
        'Unrecognized BlochHam vnp block. No KSSOLV source was modified.');
end
parts = tokens{1};
guard = regexprep(parts{2}, '\s+', '');
if strcmp(guard, new_guard)
    output = source;
elseif strcmp(guard, old_guard)
    output = [source(1:first-1), parts{1}, replacement, parts{3}, source(last+1:end)];
else
    error('enable_soc_blochham_compat:UnknownVersion', ...
        'Unrecognized BlochHam spin condition. No KSSOLV source was modified.');
end
end

function assert_resolved(expected)
resolved = which('BlochHam');
if ~strcmpi(strrep(resolved, '\', '/'), strrep(expected, '\', '/'))
    error('enable_soc_blochham_compat:ShadowedClass', ...
        'Expected BlochHam at %s, but resolved %s. Start a fresh MATLAB session.', expected, resolved);
end
end
