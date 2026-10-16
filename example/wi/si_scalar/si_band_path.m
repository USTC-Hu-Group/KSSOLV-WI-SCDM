function [path_kpoints, k_frac, path_distance] = si_band_path(C)
% FCC reciprocal fractional coordinates for C=a/2*[0 1 1;1 0 1;1 1 0].
% Distances are physical reciprocal lengths (Bohr^-1), not fractional norms.
path_kpoints = {
    0,     0,     0,    31, 'G';
    0.5,   0,     0.5,  31, 'X';
    0.5,   0.25,  0.75, 31, 'W';
    0.375, 0.375, 0.75, 31, 'K';
    0,     0,     0,    31, 'G';
    0.5,   0.5,   0.5,  31, 'L'};
nodes = cell2mat(path_kpoints(:, 1:3));
k_frac = nodes(1, :);
for i = 1:size(nodes, 1)-1
    t = (1:path_kpoints{i, 4}-1)' / (path_kpoints{i, 4}-1);
    k_frac = [k_frac; nodes(i, :) + t .* (nodes(i+1, :) - nodes(i, :))]; %#ok<AGROW>
end
reciprocal = 2*pi*inv(C)';
k_cart = k_frac * reciprocal;
path_distance = [0; cumsum(vecnorm(diff(k_cart), 2, 2))];
end
