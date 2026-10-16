function tests = test_write_w90_mmn_from_unk
tests = functiontests(localfunctions);
end

function test_identity_overlap_for_zero_g(testCase)
workdir = create_toy_case(testCase, [0 0 0]);
output_file = fullfile(workdir, 'toy.mmn');

write_w90_mmn('toy', ...
    'SCDMInput', fullfile(workdir, 'scdm.in'), ...
    'NNKPFile', fullfile(workdir, 'toy.nnkp'), ...
    'OutputFile', output_file);

blocks = read_mmn_blocks(output_file, 2, 2, 1);
verifyEqual(testCase, blocks{1}, eye(2), 'AbsTol', 1e-12);
verifyEqual(testCase, blocks{2}, eye(2), 'AbsTol', 1e-12);
end

function test_phase_overlap_for_nonzero_g(testCase)
workdir = create_toy_case(testCase, [1 0 0]);
output_file = fullfile(workdir, 'toy.mmn');

write_w90_mmn('toy', ...
    'SCDMInput', fullfile(workdir, 'scdm.in'), ...
    'NNKPFile', fullfile(workdir, 'toy.nnkp'), ...
    'OutputFile', output_file);

blocks = read_mmn_blocks(output_file, 2, 2, 1);
expected = diag([1, -1]);
verifyEqual(testCase, blocks{1}, expected, 'AbsTol', 1e-12);
verifyEqual(testCase, blocks{2}, expected, 'AbsTol', 1e-12);
end

function test_complex_spinor_overlap_and_periodic_phase(testCase)
% Both spin components contribute, including the reciprocal-boundary phase.
workdir = create_toy_case(testCase, [1 0 0]);
psi1 = [1 0; 0 1; 1i 0; 0 1i] / sqrt(2);
rotation = [1 1i; 1i 1] / sqrt(2);
psi2 = psi1 * rotation;
write_toy_spinor_unk(fullfile(workdir, 'tmp', 'UNK00001.NC'), 1, psi1);
write_toy_spinor_unk(fullfile(workdir, 'tmp', 'UNK00002.NC'), 2, psi2);
output_file = fullfile(workdir, 'toy.mmn');
write_w90_mmn('toy', 'SCDMInput', fullfile(workdir, 'scdm.in'), ...
    'NNKPFile', fullfile(workdir, 'toy.nnkp'), 'OutputFile', output_file);
blocks = read_mmn_blocks(output_file, 2, 2, 1);
phase = diag([1 -1 1 -1]);
verifyEqual(testCase, blocks{1}, psi1' * phase * psi2, 'AbsTol', 1e-10);
verifyEqual(testCase, blocks{2}, psi2' * phase * psi1, 'AbsTol', 1e-10);
end

function write_toy_spinor_unk(filename, ik, psi)
fid = fopen(filename, 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%5d %5d %5d %5d %5d\n', 2, 1, 1, ik, 2);
up = psi(1:2, :);
down = psi(3:4, :);
values = [real(up(:)), imag(up(:)), real(down(:)), imag(down(:))];
fprintf(fid, '%.16e %.16e %.16e %.16e\n', values.');
end

function workdir = create_toy_case(testCase, gvec)
workdir = char(testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder);
mkdir(fullfile(workdir, 'tmp'));

write_toy_scdm_input(fullfile(workdir, 'scdm.in'));
write_toy_nnkp(fullfile(workdir, 'toy.nnkp'), gvec);
write_toy_unk(fullfile(workdir, 'tmp', 'UNK00001.1'), 1);
write_toy_unk(fullfile(workdir, 'tmp', 'UNK00002.1'), 2);
end

function write_toy_scdm_input(filename)
fid = fopen(filename, 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '2 2 1\n');
fprintf(fid, 'scdm.eig isolated\n');
fprintf(fid, '0.0 0.0 0.0\n');
fprintf(fid, '0.5 0.0 0.0\n');
end

function write_toy_nnkp(filename, gvec)
fid = fopen(filename, 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'begin nnkpts\n');
fprintf(fid, '1\n');
fprintf(fid, '1 2 %d %d %d\n', gvec(1), gvec(2), gvec(3));
fprintf(fid, '2 1 %d %d %d\n', gvec(1), gvec(2), gvec(3));
fprintf(fid, 'end nnkpts\n');
end

function write_toy_unk(filename, ik)
fid = fopen(filename, 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%5d %5d %5d %5d %5d\n', 2, 1, 1, ik, 2);
values = [1 0; 0 0; 0 0; 1 0].';
fprintf(fid, '%16.10E %16.10E\n', values);
end

function blocks = read_mmn_blocks(filename, nbnd, nk, nntot)
fid = fopen(filename, 'r');
cleanup = onCleanup(@() fclose(fid));
fgetl(fid);
header = sscanf(fgetl(fid), '%d');
assert(isequal(header(:).', [nbnd nk nntot]));

blocks = cell(nk * nntot, 1);
for iblock = 1:numel(blocks)
    pair = sscanf(fgetl(fid), '%d');
    assert(numel(pair) == 5);
    values = fscanf(fid, '%g', [2, nbnd * nbnd]).';
    fgetl(fid);
    matrix = complex(zeros(nbnd, nbnd));
    row = 1;
    for n = 1:nbnd
        for m = 1:nbnd
            matrix(m, n) = complex(values(row, 1), values(row, 2));
            row = row + 1;
        end
    end
    blocks{iblock} = matrix;
end
end
