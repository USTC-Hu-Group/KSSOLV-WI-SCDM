function result = run_wannier90(varargin)
%RUN_WANNIER90 Prepare and optionally run a Wannier90 job after run_scdm.m.
%
% This template consumes files produced by the current KSSOLV-SCDM workflow:
%   nscf_checkpoint.mat  -> crystal structure and k-point information
%   scdm.in              -> [num_bands, nkpts, num_wann]
%   scdm.eig             -> eigenvalues in Wannier90 .eig format
%   scdm.amn             -> SCDM projection matrix in Wannier90 .amn format
%   tmp/UNKxxxxx.*       -> periodic wave functions used to build .mmn
%
% Example:
%   result = run_wannier90();
%   result = run_wannier90('SeedName','si_scdm','Run',false);
%   result = run_wannier90('W90Command','wannier90.x','OutputDirectory','w90_scdm');
%   result = run_wannier90('NumIter',10000,'DisNumIter',0);
%
% Note:
%   For Run=true, this routine calls `wannier90.x -pp seedname` first, parses
%   the generated .nnkp neighbor list, writes seedname.mmn from KSSOLV UNK
%   files, and then calls `wannier90.x seedname`.

parser = inputParser;
parser.FunctionName = mfilename;
addParameter(parser, 'WorkingDirectory', pwd, @(x) ischar(x) || isstring(x));
addParameter(parser, 'SeedName', 'wannier90', @(x) ischar(x) || isstring(x));
addParameter(parser, 'OutputDirectory', 'wannier90_post', @(x) ischar(x) || isstring(x));
addParameter(parser, 'W90Command', default_w90_command(), @(x) ischar(x) || isstring(x));
addParameter(parser, 'Run', true, @(x) islogical(x) || isnumeric(x));
addParameter(parser, 'GenerateMMN', true, @(x) islogical(x) || isnumeric(x));
addParameter(parser, 'RunPreprocess', [], @(x) isempty(x) || islogical(x) || isnumeric(x));
addParameter(parser, 'AddRandomProjections', true, @(x) islogical(x) || isnumeric(x));
addParameter(parser, 'Force', true, @(x) islogical(x) || isnumeric(x));
addParameter(parser, 'ErrorOnFailure', false, @(x) islogical(x) || isnumeric(x));
addParameter(parser, 'NumIter', 10000, ...
    @(x) isnumeric(x) && isscalar(x) && isreal(x) && isfinite(x) && x >= 0);
addParameter(parser, 'DisNumIter', 0, ...
    @(x) isnumeric(x) && isscalar(x) && isreal(x) && isfinite(x) && x >= 0);
addParameter(parser, 'ConvTol', [], ...
    @(x) isempty(x) || (isnumeric(x) && isscalar(x) && isreal(x) && isfinite(x) && x > 0));
addParameter(parser, 'ConvWindow', 3, ...
    @(x) isnumeric(x) && isscalar(x) && isreal(x) && isfinite(x) && x >= 1);
addParameter(parser, 'UseWSDistance', true, @(x) islogical(x) || isnumeric(x));
addParameter(parser, 'KmeshTol', 1e-8, ...
    @(x) isnumeric(x) && isscalar(x) && isreal(x) && isfinite(x) && x >= 0);
addParameter(parser, 'ExtraWinLines', {}, @(x) iscell(x) || isstring(x) || ischar(x));
addParameter(parser, 'WannierPlot', false, @(x) islogical(x) || isnumeric(x));
addParameter(parser, 'WannierPlotSupercell', 1, ...
    @(x) isnumeric(x) && isscalar(x) && isreal(x) && isfinite(x) && x >= 1);
addParameter(parser, 'WannierPlotSpinorMode', 'total', ...
    @(x) ischar(x) || (isstring(x) && isscalar(x)));
parse(parser, varargin{:});

case_dir = char(parser.Results.WorkingDirectory);
seedname = char(parser.Results.SeedName);
out_dir = char(parser.Results.OutputDirectory);
w90_command = char(parser.Results.W90Command);
do_run = logical(parser.Results.Run);
generate_mmn = logical(parser.Results.GenerateMMN);
run_preprocess = parser.Results.RunPreprocess;
if isempty(run_preprocess)
    run_preprocess = do_run;
else
    run_preprocess = logical(run_preprocess);
end
add_random_projections = logical(parser.Results.AddRandomProjections);
force_overwrite = logical(parser.Results.Force);
error_on_failure = logical(parser.Results.ErrorOnFailure);
num_iter = round(double(parser.Results.NumIter));
dis_num_iter = round(double(parser.Results.DisNumIter));
conv_tol = double(parser.Results.ConvTol);
conv_window = round(double(parser.Results.ConvWindow));
use_ws_distance = logical(parser.Results.UseWSDistance);
kmesh_tol = double(parser.Results.KmeshTol);
extra_win_lines = normalize_extra_win_lines(parser.Results.ExtraWinLines);
wannier_plot = logical(parser.Results.WannierPlot);
wannier_plot_supercell = round(double(parser.Results.WannierPlotSupercell));
wannier_plot_spinor_mode = lower(strtrim(char(parser.Results.WannierPlotSpinorMode)));
if ~ismember(wannier_plot_spinor_mode, {'total', 'up', 'down'})
    error('run_wannier90:BadWannierPlotSpinorMode', ...
        'WannierPlotSpinorMode must be total, up, or down.');
end
if wannier_plot
    extra_win_lines = [extra_win_lines(:); { ...
        'wannier_plot = true'; ...
        'wannier_plot_format = cube'; ...
        'wannier_plot_mode = crystal'; ...
        sprintf('wannier_plot_supercell = %d', wannier_plot_supercell); ...
        sprintf('wannier_plot_spinor_mode = %s', wannier_plot_spinor_mode); ...
        'wannier_plot_spinor_phase = false'}];
end

if ~isfolder(case_dir)
    error('run_wannier90:MissingWorkingDirectory', ...
        'WorkingDirectory does not exist: %s', case_dir);
end
if ~is_absolute_path(out_dir)
    out_dir = fullfile(case_dir, out_dir);
end
if ~isfolder(out_dir)
    mkdir(out_dir);
end

checkpoint_file = fullfile(case_dir, 'nscf_checkpoint.mat');
scdm_in_file = fullfile(case_dir, 'scdm.in');
scdm_eig_file = fullfile(case_dir, 'scdm.eig');
scdm_amn_file = fullfile(case_dir, 'scdm.amn');

require_file(checkpoint_file, 'nscf checkpoint');
require_file(scdm_in_file, 'SCDM input');
require_file(scdm_eig_file, 'Wannier90 eigenvalue input');
require_file(scdm_amn_file, 'Wannier90 AMN projection input');

loaded = load(checkpoint_file);
if ~isfield(loaded, 'crystal')
    error('run_wannier90:MissingCrystal', ...
        'nscf_checkpoint.mat must contain variable crystal.');
end
crystal = loaded.crystal;
metadata = struct();
if isfield(loaded, 'metadata')
    metadata = loaded.metadata;
end

header = read_scdm_header(scdm_in_file);
crystal = set_crystal_nbnd_for_w90(crystal, header.num_bands);
natom = infer_atom_count(crystal);
if isempty(conv_tol)
    conv_tol = 1e-5 * max(1, natom);
end

fprintf('Preparing Wannier90 input in %s\n', out_dir);
fprintf('SCDM header: num_bands=%d, nkpts=%d, num_wann=%d\n', ...
    header.num_bands, header.nkpts, header.num_wann);
fprintf(['Wannier90 settings: num_iter=%d, dis_num_iter=%d, ', ...
    'conv_tol=%s, conv_window=%d, use_ws_distance=%d, kmesh_tol=%s\n'], ...
    num_iter, dis_num_iter, format_w90_real(conv_tol), conv_window, ...
    use_ws_distance, format_w90_real(kmesh_tol));

copy_seed_file(scdm_eig_file, fullfile(out_dir, [seedname, '.eig']), force_overwrite);
copy_seed_file(scdm_amn_file, fullfile(out_dir, [seedname, '.amn']), force_overwrite);

win_prefix = fullfile(out_dir, seedname);
win_file = write_w90_win(crystal, header.num_wann, win_prefix);
append_w90_runtime_options(win_file, metadata, extra_win_lines, ...
    add_random_projections, num_iter, dis_num_iter, conv_tol, conv_window, ...
    use_ws_distance, kmesh_tol);
if wannier_plot
    stage_wannier_plot_unk_files(case_dir, out_dir);
end

nnkp_file = fullfile(out_dir, [seedname, '.nnkp']);
mmn_file = fullfile(out_dir, [seedname, '.mmn']);
mmn_info = [];
preprocess_status = NaN;
preprocess_command_output = '';
preprocess_log_file = fullfile(out_dir, [seedname, '.w90_pp.log']);

if generate_mmn
    if ~isfile(nnkp_file)
        if run_preprocess
            if isempty(strtrim(w90_command))
                error('run_wannier90:MissingW90Command', ...
                    'W90Command is empty. Cannot run Wannier90 preprocessing to create .nnkp.');
            end
            [preprocess_status, preprocess_command_output] = run_w90_preprocess( ...
                out_dir, w90_command, seedname, preprocess_log_file);
            if preprocess_status ~= 0
                msg = sprintf(['Wannier90 preprocessing returned nonzero status %d. ', ...
                    'See log file:\n%s'], preprocess_status, preprocess_log_file);
                if error_on_failure
                    error('run_wannier90:Wannier90PreprocessFailed', '%s', msg);
                else
                    warning('run_wannier90:Wannier90PreprocessFailed', '%s', msg);
                end
            end
        else
            warning('run_wannier90:NNKPNotGenerated', ...
                ['RunPreprocess=false and %s is missing. ', ...
                 'MMN generation is skipped until an nnkp file is available.'], nnkp_file);
        end
    end

    if isfile(nnkp_file)
        [~, mmn_info] = write_w90_mmn(seedname, ...
            'SCDMInput', scdm_in_file, ...
            'NNKPFile', nnkp_file, ...
            'OutputFile', mmn_file, ...
            'NumBands', header.num_bands);
    end
end

result = struct();
result.case_dir = case_dir;
result.output_dir = out_dir;
result.seedname = seedname;
result.win_file = win_file;
result.eig_file = fullfile(out_dir, [seedname, '.eig']);
result.amn_file = fullfile(out_dir, [seedname, '.amn']);
result.mmn_file = mmn_file;
result.nnkp_file = nnkp_file;
result.mmn_info = mmn_info;
result.preprocess_status = preprocess_status;
result.preprocess_command_output = preprocess_command_output;
result.preprocess_log_file = preprocess_log_file;
result.header = header;
result.num_iter = num_iter;
result.dis_num_iter = dis_num_iter;
result.conv_tol = conv_tol;
result.conv_window = conv_window;
result.use_ws_distance = use_ws_distance;
result.kmesh_tol = kmesh_tol;
result.wannier_plot = wannier_plot;
result.wannier_plot_supercell = wannier_plot_supercell;
result.wannier_plot_spinor_mode = wannier_plot_spinor_mode;
result.command = '';
result.status = NaN;
result.log_file = fullfile(out_dir, [seedname, '.w90_driver.log']);

if do_run
    if isempty(strtrim(w90_command))
        error('run_wannier90:MissingW90Command', ...
            'W90Command is empty. Set W90Command or environment variable WANNIER90.');
    end

    old_dir = pwd;
    cleanup_obj = onCleanup(@() cd(old_dir));
    cd(out_dir);

    if generate_mmn && ~isfile(mmn_file)
        warning('run_wannier90:MissingMMN', ...
            ['%s was not generated. Wannier90 may fail if it needs overlap ', ...
             'matrices for localization, disentanglement, or spread minimization.'], ...
            mmn_file);
    end

    command = sprintf('%s %s > %s 2>&1', ...
        quote_command(w90_command), quote_shell(seedname), quote_shell([seedname, '.w90_driver.log']));
    result.command = command;
    fprintf('Running Wannier90 command: %s\n', command);
    [status, cmdout] = system(command);
    result.status = status;
    result.command_output = cmdout;

    if status ~= 0
        msg = sprintf(['Wannier90 returned nonzero status %d. See log file:\n%s\n', ...
            'If the log reports a missing or inconsistent .mmn file, check ', ...
            'the preprocessing log and generated nnkp/mmn files.'], ...
            status, result.log_file);
        if error_on_failure
            error('run_wannier90:Wannier90Failed', '%s', msg);
        else
            warning('run_wannier90:Wannier90Failed', '%s', msg);
        end
    else
        fprintf('Wannier90 completed successfully. Log: %s\n', result.log_file);
    end
else
    fprintf('Run=false: Wannier90 input files were generated but not executed.\n');
end

save(fullfile(out_dir, 'wannier90_checkpoint.mat'), 'result');
end

function stage_wannier_plot_unk_files(case_dir, out_dir)
source_dir = fullfile(case_dir, 'tmp');
files = dir(fullfile(source_dir, 'UNK*'));
files = files(~[files.isdir]);
if isempty(files)
    error('run_wannier90:MissingUNKForPlot', ...
        'WannierPlot=true requires UNK files in %s.', source_dir);
end

for index = 1:numel(files)
    source = fullfile(files(index).folder, files(index).name);
    target = fullfile(out_dir, files(index).name);
    if isfile(target)
        continue;
    end
    if isunix
        [status, output] = system(sprintf('ln -sfn %s %s', ...
            quote_shell(source), quote_shell(target)));
        if status ~= 0
            error('run_wannier90:CannotLinkUNK', ...
                'Cannot link %s into Wannier90 directory: %s', source, output);
        end
    else
        [ok, message] = copyfile(source, target);
        if ~ok
            error('run_wannier90:CannotCopyUNK', ...
                'Cannot copy %s into Wannier90 directory: %s', source, message);
        end
    end
end
fprintf('Staged %d UNK files for Wannier90 Cube output in %s\n', ...
    numel(files), out_dir);
end

function command = default_w90_command()
command = getenv('WANNIER90');
if isempty(strtrim(command))
    command = 'wannier90.x';
end
end

function require_file(filename, description)
if ~isfile(filename)
    error('run_wannier90:MissingFile', ...
        'Missing %s file: %s. Run run_scdm.m first, or copy the file into this case directory.', ...
        description, filename);
end
end

function header = read_scdm_header(filename)
fid = fopen(filename, 'r');
if fid < 0
    error('run_wannier90:CannotOpenSCDMInput', ...
        'Cannot open SCDM input file: %s', filename);
end
cleanup_obj = onCleanup(@() fclose(fid));
line = fgetl(fid);
values = sscanf(line, '%f');
if numel(values) < 3
    error('run_wannier90:InvalidSCDMHeader', ...
        'The first line of scdm.in must contain num_bands, nkpts, and num_wann.');
end
header = struct();
header.num_bands = round(values(1));
header.nkpts = round(values(2));
header.num_wann = round(values(3));
end

function crystal = set_crystal_nbnd_for_w90(crystal, num_bands)
try
    crystal = set(crystal, 'nbnd', num_bands);
catch
    try
        crystal.nbnd = num_bands;
    catch
        warning('run_wannier90:CannotSetNumBands', ...
            ['Could not set crystal.nbnd to %d before write_w90_win. ', ...
             'The generated .win may use the original crystal.nbnd value.'], num_bands);
    end
end
end

function natom = infer_atom_count(crystal)
natom = [];
try
    xyzlist = crystal.xyzlist;
catch
    xyzlist = [];
end
if ~isempty(xyzlist) && size(xyzlist, 2) == 3
    natom = size(xyzlist, 1);
    return;
end

try
    atomlist = crystal.atomlist;
catch
    atomlist = [];
end
if ~isempty(atomlist)
    try
        natom = numel(atomlist);
        return;
    catch
    end
end

try
    natom = double(get(crystal, 'natom'));
catch
    natom = [];
end
if isempty(natom) || ~isscalar(natom) || ~isfinite(natom) || natom <= 0
    warning('run_wannier90:UnknownAtomCount', ...
        ['Cannot infer atom count from crystal. Using natom=1 for automatic ', ...
         'Wannier90 conv_tol; pass ConvTol manually to override.']);
    natom = 1;
end
natom = round(natom);
end

function text = format_w90_real(value)
text = sprintf('%.12e', double(value));
text = strrep(text, 'e', 'd');
text = strrep(text, 'E', 'd');
end

function text = logical_to_w90(value)
if value
    text = 'true';
else
    text = 'false';
end
end

function copy_seed_file(source_file, target_file, force_overwrite)
if isfile(target_file) && ~force_overwrite
    fprintf('Keeping existing file: %s\n', target_file);
    return;
end
[ok, msg] = copyfile(source_file, target_file, 'f');
if ~ok
    error('run_wannier90:CopyFailed', ...
        'Failed to copy %s to %s: %s', source_file, target_file, msg);
end
end

function [status, cmdout] = run_w90_preprocess(out_dir, w90_command, seedname, log_file)
old_dir = pwd;
cleanup_obj = onCleanup(@() cd(old_dir));
cd(out_dir);

[~, log_name, log_ext] = fileparts(log_file);
log_name = [log_name, log_ext];
command = sprintf('%s -pp %s > %s 2>&1', ...
    quote_command(w90_command), quote_shell(seedname), quote_shell(log_name));
fprintf('Running Wannier90 preprocess command: %s\n', command);
[status, cmdout] = system(command);
if status == 0
    fprintf('Wannier90 preprocessing completed. Log: %s\n', log_file);
end
end

function append_w90_runtime_options(win_file, metadata, extra_win_lines, ...
        add_random_projections, num_iter, dis_num_iter, conv_tol, ...
        conv_window, use_ws_distance, kmesh_tol)
fid = fopen(win_file, 'a');
if fid < 0
    error('run_wannier90:CannotAppendWin', ...
        'Cannot append Wannier90 options to %s', win_file);
end
cleanup_obj = onCleanup(@() fclose(fid));

fprintf(fid, '\n! Additional options appended by run_wannier90.m\n');
fprintf(fid, '! SCDM supplies the disentangled subspace; Wannier90 only refines the gauge.\n');
fprintf(fid, 'num_iter = %d\n', num_iter);
fprintf(fid, 'dis_num_iter = %d\n', dis_num_iter);
fprintf(fid, 'conv_tol = %s\n', format_w90_real(conv_tol));
fprintf(fid, 'conv_window = %d\n', conv_window);
fprintf(fid, 'use_ws_distance = %s\n', logical_to_w90(use_ws_distance));
fprintf(fid, 'kmesh_tol = %s\n', format_w90_real(kmesh_tol));
fprintf(fid, 'write_hr = true\n');
fprintf(fid, 'write_xyz = true\n');

if add_random_projections
    fprintf(fid, '\n! Ordinary projection block for Wannier90 preprocessing.\n');
    fprintf(fid, '! The actual SCDM projections are supplied by seedname.amn.\n');
    fprintf(fid, 'begin projections\n');
    fprintf(fid, 'random\n');
    fprintf(fid, 'end projections\n');
end

if has_path_kpoints(metadata)
    fprintf(fid, 'bands_plot = true\n');
    fprintf(fid, 'begin kpoint_path\n');
    path_kpoints = metadata.path_kpoints;
    for i = 1:(size(path_kpoints, 1) - 1)
        label_a = sanitize_kpath_label(path_kpoints{i, 5});
        label_b = sanitize_kpath_label(path_kpoints{i + 1, 5});
        ka = cell2mat(path_kpoints(i, 1:3));
        kb = cell2mat(path_kpoints(i + 1, 1:3));
        fprintf(fid, '%-8s %16.10f %16.10f %16.10f   %-8s %16.10f %16.10f %16.10f\n', ...
            label_a, ka(1), ka(2), ka(3), label_b, kb(1), kb(2), kb(3));
    end
    fprintf(fid, 'end kpoint_path\n');
end

for i = 1:numel(extra_win_lines)
    fprintf(fid, '%s\n', extra_win_lines{i});
end
end

function tf = has_path_kpoints(metadata)
tf = isstruct(metadata) && isfield(metadata, 'path_kpoints') && ...
    iscell(metadata.path_kpoints) && size(metadata.path_kpoints, 1) >= 2 && ...
    size(metadata.path_kpoints, 2) >= 5;
end

function label = sanitize_kpath_label(label)
label = char(string(label));
if strcmp(label, '\Gamma') || strcmpi(label, 'Gamma') || strcmp(label, char(915))
    label = 'G';
end
label = regexprep(label, '\s+', '');
if isempty(label)
    label = 'K';
end
end

function lines = normalize_extra_win_lines(lines)
if ischar(lines)
    lines = cellstr(lines);
elseif isstring(lines)
    lines = cellstr(lines(:));
elseif isempty(lines)
    lines = {};
elseif ~iscell(lines)
    error('run_wannier90:InvalidExtraWinLines', ...
        'ExtraWinLines must be a char, string array, or cell array of char.');
end
end

function tf = is_absolute_path(pathname)
pathname = char(pathname);
if ispc
    tf = numel(pathname) >= 2 && pathname(2) == ':';
else
    tf = startsWith(pathname, filesep);
end
end

function quoted = quote_command(command)
command = char(command);
if isfile(command)
    quoted = quote_shell(command);
else
    quoted = command;
end
end

function quoted = quote_shell(text)
text = char(text);
if isempty(text)
    quoted = '""';
    return;
end
quoted = ['"', strrep(text, '"', '\"'), '"'];
end
