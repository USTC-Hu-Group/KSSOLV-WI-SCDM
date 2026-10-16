function audit_wi_scdm_f_inverse
% Numerical characterization only; no material interpolation is rerun here.
repo = fileparts(fileparts(fileparts(mfilename('fullpath'))));
old_path = path; old_warnings = warning;
cleanup = onCleanup(@() restore_environment(old_path, old_warnings));
addpath(fullfile(repo, 'src', 'WI'), '-begin');
out = fullfile(repo, 'tests', 'wi', '_work', 'inverse');
if ~isfolder(out), mkdir(out); end
rows = struct('n',{},'a_eV',{},'relative_gap',{},'energy_eV',{}, ...
    'transformed_eV',{},'inverse_eV',{},'physical_error_eV',{}, ...
    'transformed_residual_eV',{},'absolute_inverse_condition',{},'warning_id',{});
for n = [0.5,3]
    for a = [0.1,1,10]
        params = wi_scdm_f_params([-2,0;1,2], 'TransformN', n, ...
            'TransformA', a, 'TransformEpsilon', 0);
        for gap = [2,1,0.1,0.01,1e-4,1e-6,1e-8]
            E = -a*gap;
            lambda = wi_scdm_f_transform(E, params);
            lastwarn('');
            [recovered, diagnostics] = wi_scdm_f_inverse(lambda, params);
            [~,warning_id] = lastwarn;
            derivative = 1;
            if gap < 1
                derivative = 0.5-erf(n*(0.5-gap))/(2*erf(n/2));
            end
            residual = abs(wi_scdm_f_transform(recovered,params)-lambda);
            assert(abs(diagnostics.max_residual_eV - ...
                abs(wi_scdm_f_transform(recovered,params)-min(lambda,0)))<1e-20);
            assert(diagnostics.converged && residual <= 1.1e-14);
            assert(isfinite(recovered) && recovered<=1e-12);
            if gap>=0.1
                assert(abs(recovered-E)<1e-8*max(a,1));
            end
            rows(end+1) = struct('n',n,'a_eV',a,'relative_gap',gap, ...
                'energy_eV',E,'transformed_eV',lambda,'inverse_eV',recovered, ...
                'physical_error_eV',abs(recovered-E), ...
                'transformed_residual_eV',residual, ...
                'absolute_inverse_condition',1/derivative,'warning_id',warning_id); %#ok<AGROW>
        end
    end
end
writetable(struct2table(rows),fullfile(out,'round_trip_conditioning.csv'));

lastwarn('');
flat = wi_scdm_f_params([-1,-1;2,2]);
[~,flat_warning] = lastwarn;
assert(strcmp(flat_warning,'wi_scdm_f_params:FlatTopBand'));
assert(abs(flat.a-0.15)<1e-14);
lastwarn('');
near = wi_scdm_f_params([-1,-1;2,2+1e-12]);
[~,near_warning] = lastwarn;
assert(isempty(near_warning) && near.a>0 && near.a<1e-10);
params = wi_scdm_f_params([-2,0;1,2], 'TransformA',1,'TransformEpsilon',0);
lastwarn(''); small_positive = wi_scdm_f_inverse(1e-9,params);
[~,small_warning] = lastwarn;
assert(small_positive==0 && isempty(small_warning));
lastwarn(''); large_positive = wi_scdm_f_inverse(1e-4,params);
[~,large_warning] = lastwarn;
assert(large_positive==0 && strcmp(large_warning, ...
    'wi_scdm_f_inverse:PositiveTransformedEnergy'));
assert(all(wi_scdm_f_transform([0,1,2],params)==0));
assert(wi_scdm_f_inverse(0,params)==0);
mixed_params=wi_scdm_f_params([-2,0;1,2], 'TransformN',3, ...
    'TransformA',0.1,'TransformEpsilon',0);
mixed_E=[-0.2,-0.01,-1e-7,0];
mixed_lambda=wi_scdm_f_transform(mixed_E,mixed_params);
[mixed_return,mixed_info]=wi_scdm_f_inverse(mixed_lambda,mixed_params);
mixed_residual=abs(wi_scdm_f_transform(mixed_return,mixed_params)-mixed_lambda);
assert(mixed_info.converged && max(mixed_residual)<1.1e-14);
assert(abs(mixed_return(3)-mixed_E(3))<=1e-7 && mixed_return(4)==0);
edges = struct('flat_width_eV',flat.a,'flat_warning',flat_warning, ...
    'nearly_flat_width_eV',near.a,'nearly_flat_warning',near_warning, ...
    'small_positive_target_eV',1e-9,'small_positive_return_eV',small_positive, ...
    'small_positive_warning',small_warning,'large_positive_target_eV',1e-4, ...
    'large_positive_return_eV',large_positive,'large_positive_warning',large_warning);
edges.mixed_vector_residual_eV=max(mixed_residual);
summary = struct('matlab_version',version,'round_trip_probes',numel(rows), ...
    'scope','Scalar-function audit, not a material sensitivity or accuracy benchmark', ...
    'max_physical_error_eV',max([rows.physical_error_eV]), ...
    'max_transformed_residual_eV',max([rows.transformed_residual_eV]),'edges',edges);
fid=fopen(fullfile(out,'summary.json'),'w'); closefile=onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(summary,'PrettyPrint',true)); clear closefile
fprintf('PASS: %d scalar probes and flat/near-flat/plateau/clipping checks.\n',numel(rows));
fprintf('Max physical error %.6g eV, transformed residual %.6g eV.\n', ...
    summary.max_physical_error_eV,summary.max_transformed_residual_eV);
end

function restore_environment(old_path,old_warnings)
path(old_path); warning(old_warnings);
end
