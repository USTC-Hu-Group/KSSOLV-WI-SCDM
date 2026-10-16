function test_scdm_provenance
% Synthetic branch tests; these data are not manuscript benchmark results.
repo=fileparts(fileparts(fileparts(mfilename('fullpath')))); addpath(fullfile(repo, 'src', 'WI'), '-begin');
parent=fullfile(repo,'tests','wi','_work'); if ~isfolder(parent), mkdir(parent); end
work=tempname(parent); mkdir(work);
old=pwd; cleanup=onCleanup(@() cd(old)); cd(work);
crystal=struct('nkpts',1,'nbnd',4,'kpts',[0 0 0], ...
    'supercell',eye(3)*8,'xyzlist',[0 0 0],'alist',1, ...
    'atoms',struct('symbol','Si','anum',14,'venum',4));
write_input_files(crystal,[-1;0;1;2],2,'erfc',0,2);
files=dir(fullfile(work,'scdm_provenance','run_*','metadata.json'));
assert(numel(files)==1);
meta=jsondecode(fileread(fullfile(files(1).folder,files(1).name)));
assert(strcmp(meta.status,'success') && meta.manual_parameters && meta.num_wann==2);
assert(meta.mu_eV==0 && meta.sigma_eV==2);
original=fileread('scdm.in');
try
    write_input_files(crystal,[-1;0;1;2],5,'erfc',0,2);
    error('test_scdm_provenance:ExpectedFailure','Expected invalid target to fail.');
catch err
    assert(strcmp(err.identifier,'write_input_files:TooManyWannierFunctions'));
end
assert(strcmp(original,fileread('scdm.in')));
files=dir(fullfile(work,'scdm_provenance','run_*','metadata.json'));
assert(numel(files)==2);
states=cell(1,2);
for i=1:2
    meta=jsondecode(fileread(fullfile(files(i).folder,files(i).name)));
    states{i}=meta.status;
    if strcmp(meta.status,'failed')
        assert(strcmp(meta.error_identifier,'write_input_files:TooManyWannierFunctions'));
    end
end
assert(all(ismember({'success','failed'},states)));

record=struct('status','success','energy_eV',[0;1;2;3], ...
    'projectability',[1;.8;.3;0],'fit_result',[1.5,1], ...
    'pao_info',struct('raw_projectability',[1;.8;.3;0], ...
       'singular_values_by_k',{{[2;1;1e-14]}},'rank_tolerance_by_k',2e-10));
folder=save_scdm_provenance(record,fullfile(work,'synthetic_fit_archive'));
assert(isfile(fullfile(folder,'metadata.json')));
tab=readtable(fullfile(folder,'fit_input_output.csv'));
assert(height(tab)==4 && all(tab.used_in_fit));
tab=readtable(fullfile(folder,'pao_singular_values.csv'));
assert(isequal(tab.retained,[1;1;0]));
fprintf('PASS: provenance success, failure/rethrow, immutable runs, CSV and SVD records.\n');
fprintf('Synthetic test artifacts: %s\n',work);
end
