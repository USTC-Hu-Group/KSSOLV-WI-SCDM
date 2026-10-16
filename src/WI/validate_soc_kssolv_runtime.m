function validate_soc_kssolv_runtime(crystal)
%VALIDATE_SOC_KSSOLV_RUNTIME Detect the known nonmagnetic SOC vnp bug early.
% domag=0 has scalar XC/local potentials but still two-component spinors.
% Do not switch domag to 1: that would enable a different magnetic model.
if crystal.nspin ~= 4 || crystal.domag
    return;
end
filename = which('BlochHam');
fprintf('SOC runtime: BlochHam -> %s\n', filename);
source = regexprep(fileread(filename), '\s+', '');
if contains(source, 'ifcry.nspin==1BH.vnp=vhart+vxc;')
    error('validate_soc_kssolv_runtime:NonmagneticSOCBranchBug', ...
        ['This BlochHam lacks its nspin=4/domag=0 scalar-potential branch. ', ...
         'Start a new MATLAB process and call enable_soc_blochham_compat ', ...
         'after KSSOLV_startup, before loading checkpoints or creating Hamiltonians. ', ...
         'Affected file: %s'], filename);
end
end
