function folder = save_scdm_provenance(record, archive_root)
%SAVE_SCDM_PROVENANCE Save one immutable postprocessing attempt per directory.
% Logging failures warn without replacing the numerical result/exception.
if nargin < 2, archive_root = fullfile(pwd, 'scdm_provenance'); end
folder = '';
try
    if ~isfolder(archive_root), mkdir(archive_root); end
    [~, token] = fileparts(tempname);
    folder = fullfile(archive_root, ['run_', token]);
    mkdir(folder);
    record.schema_version = 1;
    record.matlab_version = version;
    save(fullfile(folder,'record.mat'),'record','-v7');

    meta = record;
    arrays = {'energy_eV','projectability','valid_fit_mask','fit_residual'};
    for n=1:numel(arrays)
        if isfield(meta,arrays{n}), meta=rmfield(meta,arrays{n}); end
    end
    if isfield(meta,'pao_info')
        for name = {'channels_by_species','raw_projectability','singular_values_by_k','pw_indices'}
            if isfield(meta.pao_info,name{1}), meta.pao_info=rmfield(meta.pao_info,name{1}); end
        end
    end
    meta.data_note = 'Full arrays and PAO radial channels are in record.mat; CSVs are readable exports.';
    manifest = struct('role',{},'source',{},'file',{},'sha256',{});
    if isfield(record,'pao_info') && isfield(record.pao_info,'pp_files')
        ppdir=fullfile(folder,'pseudopotentials'); mkdir(ppdir);
        for n=1:numel(record.pao_info.pp_files)
            source=record.pao_info.pp_files{n};
            [~,base,ext]=fileparts(source);
            relative=fullfile('pseudopotentials',sprintf('%03d_%s%s',n,base,ext));
            copyfile(source,fullfile(folder,relative));
            manifest(end+1)=entry('pseudopotential',source,relative,folder); %#ok<AGROW>
        end
    end
    for name={'compute_pao_projectability','write_input_files','save_scdm_provenance','ylmr','bessel','simpson','Ggrid','kssolvpptype'}
        source=which(name{1});
        if ~isempty(source) && isfile(source)
            [~,base,ext]=fileparts(source);
            relative=['source_',base,ext];
            copyfile(source,fullfile(folder,relative));
            manifest(end+1)=entry('code',source,relative,folder); %#ok<AGROW>
        end
    end
    if strcmp(record.status,'success')
        for name={'scdm.in','scdm.eig'}
            if isfile(name{1})
                copyfile(name{1},fullfile(folder,name{1}));
                manifest(end+1)=entry('output',fullfile(pwd,name{1}),name{1},folder); %#ok<AGROW>
            end
        end
    end
    meta.files=manifest;
    if isfield(record,'energy_eV') && isfield(record,'projectability') && ...
            isequal(size(record.energy_eV),size(record.projectability)) && ...
            ~isempty(record.projectability)
        [nb,nk]=size(record.energy_eV);
        energy=record.energy_eV(:); projectability=record.projectability(:);
        raw=nan(size(energy));
        if isfield(record.pao_info,'raw_projectability') && ...
                numel(record.pao_info.raw_projectability)==numel(energy)
            raw=record.pao_info.raw_projectability(:);
        end
        mask=isfinite(energy)&isfinite(projectability);
        predicted=nan(size(energy));
        if isfield(record,'fit_result')
            fit=record.fit_result;
            predicted=0.5*erfc((energy-fit(1))/fit(2));
        end
        tab=table(repmat((1:nb)',nk,1),reshape(repelem(1:nk,nb),[],1),energy,raw, ...
            projectability,mask,predicted,predicted-projectability, ...
            'VariableNames',{'band_index','k_index','energy_eV','projectability_raw', ...
            'projectability','used_in_fit','fit_prediction','fit_minus_observation'});
        writetable(tab,fullfile(folder,'fit_input_output.csv'));
    end
    if isfield(record,'pao_info') && isfield(record.pao_info,'singular_values_by_k')
        rows=zeros(0,5);
        for k=1:numel(record.pao_info.singular_values_by_k)
            sv=record.pao_info.singular_values_by_k{k}(:);
            tau=record.pao_info.rank_tolerance_by_k(k);
            rows=[rows;[repmat(k,numel(sv),1),(1:numel(sv))',sv, ...
                repmat(tau,numel(sv),1),double(sv>tau)]]; %#ok<AGROW>
        end
        writetable(array2table(rows,'VariableNames', ...
            {'k_index','singular_index','singular_value','threshold','retained'}), ...
            fullfile(folder,'pao_singular_values.csv'));
    end
    fid=fopen(fullfile(folder,'metadata.json'),'w');
    if fid<0, error('Cannot write provenance metadata.'); end
    cleanup=onCleanup(@() fclose(fid));
    fprintf(fid,'%s\n',jsonencode(meta,'PrettyPrint',true));
    clear cleanup;
    fprintf('SCDM provenance saved: %s\n',folder);
catch err
    warning('save_scdm_provenance:ArchiveFailed', ...
        'SCDM provenance could not be fully saved: %s',err.message);
end
end

function item=entry(role,source,relative,folder)
fid=fopen(fullfile(folder,relative),'rb');
if fid<0, error('Cannot hash %s.',relative); end
cleanup=onCleanup(@() fclose(fid));
bytes=fread(fid,inf,'*uint8');
digest=java.security.MessageDigest.getInstance('SHA-256');
digest.update(bytes);
hash=lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[]));
item=struct('role',role,'source',source,'file',relative,'sha256',hash);
end
