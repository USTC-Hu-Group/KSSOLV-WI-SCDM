function report = run_analytic_demo
%RUN_ANALYTIC_DEMO Small scalar/spinor file workflow with analytic energies.
% Uses synthetic states, manual filtering and no KSSOLV or Wannier90 solver.
% This is a software demonstration, not a material or paper benchmark.
case_dir = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(fileparts(case_dir)));
old_dir = pwd; old_path = path; old_rng = rng;
restore = onCleanup(@() restore_environment(old_dir,old_path,old_rng));
addpath(root,'-begin'); KSSOLV_startup;
parent=fullfile(case_dir,'_output'); if ~isfolder(parent), mkdir(parent); end
work=tempname(parent); mkdir(work);
rng(17,'twister');
k_frac=[0 0 0; 1/3 0 0; 2/3 0 0];
E=[-3+0.3*cos(2*pi*k_frac(:,1)).'; ...
   -1+0.2*sin(2*pi*k_frac(:,1)).'; ...
    8+0.1*cos(2*pi*k_frac(:,1)).'];
records=struct('mode',{},'orthogonality_error',{}, ...
    'node_error_wi_scdm_eV',{},'node_error_wi_scdm_f_eV',{});
for spinor=0:1
    mode='scalar'; if spinor, mode='spinor'; end
    folder=fullfile(work,mode); mkdir(folder); mkdir(fullfile(folder,'tmp')); cd(folder);
    n_grid=4; nb=3; nk=size(k_frac,1); nw=2;
    [basis,~]=qr(randn((1+spinor)*n_grid,nb)+1i*randn((1+spinor)*n_grid,nb),0);
    for ik=1:nk
        angle=0.2*sin(2*pi*k_frac(ik,1));
        rotation=[cos(angle),-sin(angle),0; sin(angle),cos(angle),0; 0,0,1];
        psi=basis*rotation;
        write_toy_unk(ik,psi,n_grid,spinor);
    end
    C=8*eye(3);
    crystal=struct('nkpts',nk,'nbnd',nb,'nkxyz',[3 1 1], ...
        'kpts',k_frac*(2*pi*inv(C).'),'supercell',C, ...
        'xyzlist',[0 0 0],'alist',1,'efermi',0,'nspin',1+3*spinor, ...
        'atoms',struct('symbol','Si','anum',14,'venum',4));
    % Manual target and filter skip PAO construction and nonlinear fitting.
    write_input_files(crystal,E/27.211386246,nw,'erfc',0,0.5);
    [wan_r,Amn,cube_files]=scdm_m('scdm.in','scdm.amn');
    orthogonality_error=0;
    for ik=1:nk
        orthogonality_error=max(orthogonality_error,norm(Amn{ik}'*Amn{ik}-eye(nw),'fro'));
    end
    query=[k_frac; linspace(0,1,25).',zeros(25,2)];
    bands=compute_scdm_only_bands('scdm.amn','scdm.eig','scdm.in',query,0,'Method','wi-scdm');
    bands_f=compute_scdm_only_bands('scdm.amn','scdm.eig','scdm.in',query,0, ...
        'Method','wi-scdm-f','TransformA',4,'TransformEpsilon',0);
    expected=sort(E(1:nw,:),1).';
    node_error=max(abs(bands(1:nk,:)-expected),[],'all');
    node_error_f=max(abs(bands_f(1:nk,:)-expected),[],'all');
    assert(orthogonality_error<1e-10 && node_error<1e-8 && node_error_f<1e-8);
    % Exercise the actual overlap writer on a small, explicit neighbour list.
    % This single-neighbour fixture is not a complete W90 finite-difference mesh.
    fid=fopen('toy.nnkp','w');
    fprintf(fid,'begin nnkpts\n1\n1 2 0 0 0\n2 3 0 0 0\n3 1 1 0 0\nend nnkpts\n');
    fclose(fid);
    write_w90_mmn('toy','SCDMInput','scdm.in','NNKPFile','toy.nnkp');
    write_w90_win(crystal,nw,'toy','KpointCoordinates','cart','AtomCoordinates','cart');
    calculate_scdm_spread('CubePattern','scdm_*.cube','OutputPrefix','spread');
    save('demo_results.mat','wan_r','Amn','cube_files','query','bands','bands_f','E');
    records(end+1)=struct('mode',mode,'orthogonality_error',orthogonality_error, ...
        'node_error_wi_scdm_eV',node_error,'node_error_wi_scdm_f_eV',node_error_f); %#ok<AGROW>
end
report=struct('kind','synthetic software demonstration','matlab_version',version, ...
    'output_directory',work,'records',records);
fid=fopen(fullfile(work,'summary.json'),'w');
closer=onCleanup(@() fclose(fid)); fprintf(fid,'%s\n',jsonencode(report,'PrettyPrint',true)); clear closer;
fprintf('PASS: scalar/spinor SCDM, WI-SCDM, WI-SCDM-f and overlap export.\nOutputs: %s\n',work);
end
function write_toy_unk(ik,psi,n_grid,spinor)
if spinor, suffix='NC'; else, suffix='1'; end
fid=fopen(fullfile('tmp',sprintf('UNK%05d.%s',ik,suffix)),'w');
closer=onCleanup(@() fclose(fid));
fprintf(fid,'%d 1 1 %d %d\n',n_grid,ik,size(psi,2));
up=psi(1:n_grid,:);
if spinor
    down=psi(n_grid+1:end,:);
    values=[real(up(:)),imag(up(:)),real(down(:)),imag(down(:))];
else
    values=[real(up(:)),imag(up(:))];
end
fprintf(fid,[repmat('%.17g ',1,size(values,2)),'\n'],values.');
end
function restore_environment(old_dir,old_path,old_rng)
cd(old_dir); path(old_path); rng(old_rng);
end
