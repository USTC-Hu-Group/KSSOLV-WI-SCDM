function report = test_scdm_shifted_kgrid()
% Test the real file-to-band path against a dense trigonometric interpolant.
% The reference uses cardinal weights, not FFT indexing or origin correction.
repo = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(repo, 'src', 'WI'), '-begin');
if ~isfolder(fullfile(repo, 'tests', 'wi', '_work')), mkdir(fullfile(repo, 'tests', 'wi', '_work')); end
work = tempname(fullfile(repo, 'tests', 'wi', '_work'));
mkdir(work);
rng(903, 'twister');
specs = { ...
    'gamma_even', [4 1 1], [0 0 0]; ...
    'even_mp', [4 1 1], [-.125 0 0]; ...
    'shifted_odd', [5 3 1], [.08 -.07 .23]; ...
    'shifted_mixed_3d', [4 3 2], [.125 .1 .25]; ...
    'nyquist_3d', [2 2 2], [.25 -.25 .25]; ...
    'single_shifted', [1 1 1], [.2 .1 -.3]};
report = struct('case_results', [], 'work_directory', work);
for t = 1:size(specs,1)
    name = specs{t,1}; dims = specs{t,2}; delta = specs{t,3};
    [x,y,z] = ndgrid((0:dims(1)-1)/dims(1), ...
        (0:dims(2)-1)/dims(2), (0:dims(3)-1)/dims(3));
    k = [x(:),y(:),z(:)] + delta;
    k = k - floor(k+.5);
    nk = size(k,1); nb = 3;
    E = zeros(nb,nk); A = zeros(nb,nb,nk); H = zeros(nb,nb,nk);
    for ik=1:nk
        Z=randn(nb)+1i*randn(nb);
        H(:,:,ik)=diag([0,2,4])+.3*(Z+Z');
        [V,D]=eig(H(:,:,ik));
        E(:,ik)=diag(D); A(:,:,ik)=V';
    end
    query=[k; rand(9,3)-.5; delta; delta+1./(2*dims)];
    folder=fullfile(work,name); mkdir(folder);
    write_case(folder,k,E,A,dims,[]);
    [kr,dr,indices,shift]=read_scdm_uniform_kgrid(fullfile(folder,'scdm.in'),nk);
    assert(isequal(dr,dims) && max(abs(kr-k),[],'all')<1e-12);
    assert(all(shift.*dims>=-.5-1e-12) && all(shift.*dims<.5));
    assert(numel(unique(indices))==nk);

    for transformed=0:1
        method='wi-scdm'; Href=H;
        args={'Method',method};
        if transformed
            method='wi-scdm-f';
            params=wi_scdm_f_params(E,'TransformA',10,'TransformEpsilon',8);
            T=wi_scdm_f_transform(E,params);
            for ik=1:nk, Href(:,:,ik)=A(:,:,ik)'*diag(T(:,ik))*A(:,:,ik); end
            args={'Method',method,'TransformParams',params};
        end
        evalc('[bands,HR,info]=compute_scdm_only_bands(fullfile(folder,''scdm.amn''),fullfile(folder,''scdm.eig''),fullfile(folder,''scdm.in''),query,0,args{:});');
        expected=zeros(size(bands));
        for iq=1:size(query,1)
            w=ones(nk,1);
            for axis=1:3
                w=w.*cardinal_weight(query(iq,axis)-k(:,axis),dims(axis));
            end
            HQ=sum(Href.*reshape(w,1,1,[]),3);
            values=eig((HQ+HQ')/2);
            if transformed, values=wi_scdm_f_inverse(values,params); end
            expected(iq,:)=sort(real(values));
        end
        band_error=max(abs(bands-expected),[],'all');
        node_error=max(abs(bands(1:nk,:)-sort(E,1).'),[],'all');
        assert(band_error<2e-8 && node_error<2e-8);
        assert(isequal(info.kgrid_dims,dims));
        assert(max(abs(info.kgrid_shift_fractional-shift))<1e-12);

        % Independently integrate the actual k coordinates for each returned R.
        [r1,r2,r3]=ndgrid(centered_indices(dims(1)), ...
            centered_indices(dims(2)),centered_indices(dims(3)));
        vectors=[r1(:),r2(:),r3(:)];
        returned=reshape(HR,nk,nb,nb); coefficient_error=0;
        for ir=1:nk
            phase=exp(-2i*pi*(k*vectors(ir,:)'))/nk;
            direct=sum(Href.*reshape(phase,1,1,[]),3);
            coefficient_error=max(coefficient_error, ...
                norm(squeeze(returned(ir,:,:))-direct,'fro'));
        end
        assert(coefficient_error<1e-10);

        % Reorder positive/negative half-step points, including the first row.
        order=randperm(nk);
        shuffled=fullfile(folder,['shuffled_',method]); mkdir(shuffled);
        write_case(shuffled,k(order,:),E(:,order),A(:,:,order),dims,[]);
        evalc('[reordered,HRreordered]=compute_scdm_only_bands(fullfile(shuffled,''scdm.amn''),fullfile(shuffled,''scdm.eig''),fullfile(shuffled,''scdm.in''),query,0,args{:});');
        order_error=max(abs(reordered-bands),[],'all');
        assert(order_error<2e-8 && max(abs(HRreordered-HR),[],'all')<1e-10);
        entry=struct('name',name,'method',method,'node_error_eV',node_error, ...
            'dense_reference_error_eV',band_error,'H_R_error',coefficient_error, ...
            'reordering_error_eV',order_error);
        report.case_results=[report.case_results;entry]; %#ok<AGROW>
    end
end

% Legacy Cartesian reciprocal coordinates on a nonorthogonal lattice.
dims=[4 3 2]; delta=[-.125,.07,-.25];
[x,y,z]=ndgrid((0:3)/4,(0:2)/3,(0:1)/2);
k=[x(:),y(:),z(:)]+delta;
cell_bohr=[4 0 0;1 5 0;.2 .3 6]; reciprocal=2*pi*inv(cell_bohr).';
folder=fullfile(work,'cartesian_input'); mkdir(folder);
write_case(folder,k*reciprocal,zeros(1,24),ones(1,1,24),dims,cell_bohr);
[recovered,~,~,shift]=read_scdm_uniform_kgrid(fullfile(folder,'scdm.in'),24);
assert(max(abs(recovered-k),[],'all')<1e-12);
assert(max(abs(shift-delta))<1e-12);

% Invalid duplicate/missing points remain rejected.
bad=k; bad(end,:)=bad(1,:);
folder=fullfile(work,'duplicate_points'); mkdir(folder);
write_case(folder,bad,zeros(1,24),ones(1,1,24),dims,[]);
rejected=false;
try
    read_scdm_uniform_kgrid(fullfile(folder,'scdm.in'),24);
catch err
    rejected=strcmp(err.identifier,'read_scdm_uniform_kgrid:InvalidUniformKGrid');
end
assert(rejected);
report.cartesian_conversion_passed=true;
report.invalid_grid_rejected=true;
fprintf('PASS: %d band/FFT cases, reordered meshes, Cartesian conversion, invalid-grid rejection.\n',numel(report.case_results));
end

function w=cardinal_weight(x,N)
w=ones(size(x));
for r=1:floor((N-1)/2), w=w+2*cos(2*pi*r*x); end
if mod(N,2)==0, w=w+cos(pi*N*x); end
w=w/N;
end

function r=centered_indices(N)
r=0:N-1; r(r>=N/2)=r(r>=N/2)-N;
end

function write_case(folder,k,E,A,dims,cell_bohr)
[nb,nk]=size(E); nw=size(A,2);
fid=fopen(fullfile(folder,'scdm.in'),'w');
fprintf(fid,'%d %d %d\nscdm.eig isolated\n',nb,nk,nw);
fprintf(fid,'%.17g %.17g %.17g\n',k.');
fprintf(fid,'KGRID %d %d %d\n',dims);
if ~isempty(cell_bohr), fprintf(fid,'CELL_BOHR\n');fprintf(fid,'%.17g %.17g %.17g\n',cell_bohr.');end
fclose(fid);
fid=fopen(fullfile(folder,'scdm.eig'),'w');
for ik=1:nk
    for ib=1:nb, fprintf(fid,'%d %d %.17g\n',ib,ik,E(ib,ik));end
end
fclose(fid);
fid=fopen(fullfile(folder,'scdm.amn'),'w');
fprintf(fid,'Analytic grid regression\n%d %d %d\n',nb,nk,nw);
for ik=1:nk
    for iw=1:nw
        for ib=1:nb
            fprintf(fid,'%d %d %d %.17g %.17g\n',ib,iw,ik,real(A(ib,iw,ik)),imag(A(ib,iw,ik)));
        end
    end
end
fclose(fid);
end
