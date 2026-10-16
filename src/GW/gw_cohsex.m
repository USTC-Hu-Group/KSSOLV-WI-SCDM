function [Esx_x, Ech, Ex] = gw_cohsex(ksinfo, omega, eta, nvbands, ncbands, NBANDS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ksinfo is a structure that contains ground state calculation results
% and parameters of the molecules
%   nv      --- number of valence states
%   Z       --- contains the eigenvecgtors from the KSDFT calculation or Quantum-Espresso calculation
%   ev      --- contains the corresponding eigenvalues
%   vol     --- volume of the unit cell
%   ntot    --- total number of grid points on which the wavefunction is sampled
% 
% eta     --- Lorentzian broadening factor that turns a Dirac delta into a 
%             a smoother peak
% nvbands --- the number of valence bands (KS orbitals) from the Fermi 
%             level included in the kernel calculation. 
% ncbands --- the number of empty bands (KS orbitals) included in 
%             the kernel calculation
% NBANDS  --- the number of empty bands (KS orbitals) included in
%             the Coulomb-hole (Ech) calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialization
nv    = ksinfo.nv;
Z     = ksinfo.Z;
ev    = ksinfo.ev;
F     = ksinfo.F;
vol   = ksinfo.vol;
Vxc   = ksinfo.Vxc;
nbands = nvbands + ncbands;
ry2ev = 13.60569253;

im = sqrt(-1);
[ng,nr] = size(F);

fprintf('nr = %d, ng = %d, nbands = %d\n', nr, ng, nbands);

% Normalize the wavefunction in the Fourier space
% 1/Vol sum_G |Psi(G)|^2 = 1
tic

if nargin < 6
  for iv = 1:nbands
    Z(:,iv) = Z(:,iv) * sqrt(vol) / (norm(Z(:,iv)));
  end
else
  for iv = 1:NBANDS
    Z(:,iv) = Z(:,iv) * sqrt(vol) / (norm(Z(:,iv)));
  end
end

timeforZ = toc

if(0)

tic

for n = 1:nbands
    psir(:,n) = F'*Z(:,nv-nvbands+n);
end

timeforC2R = toc

chi0 = zeros(ng);

tic

for iv = nv-nvbands+1:nv
  for jc = nv+1:nv-nvbands+nbands
    eden = 2.0/sqrt(ev(iv)-ev(jc));
    Mg(:,jc-nvbands) = F*(conj(psir(:,iv)).*psir(:,jc)*eden);
  end
  chi0 = chi0 - Mg*Mg'/vol;
end

timeforchi0 = toc

tic

I = eye(ng);
Dcoul = spdiags(ksinfo.coulG(:,4), 0, ng, ng);
inveps = inv(I - Dcoul*chi0);

timeforepsilon = toc

Esx_x = zeros(nbands, nbands);
Ex = zeros(nbands, nbands);
Ex_tmp = zeros(nbands,nbands);
tic

Dcoul(1,1) = ksinfo.coulG0;
inveps_1 = inveps - I;
W1 = inveps_1*Dcoul;
for n1 = nv-nvbands+1:nv-nvbands+nbands
  for n2 = nv-nvbands+1:nv-nvbands+nbands
    for iocc = nv-nvbands+1:nv
      psioccn1 = F*(psir(:,iocc).*conj(psir(:,n1)));
      psioccn2 = F*(psir(:,iocc).*conj(psir(:,n2)));
      Esx_x(n1, n2) = Esx_x(n1, n2) + psioccn1'*(W1*psioccn2);
      Ex_tmp = psioccn1'*(Dcoul*psioccn2);
      Ex(n1, n2) = Ex(n1, n2) + psioccn1'*(Dcoul*psioccn2);
    end
  end
end

timeforXandSX = toc

tic

if nargin < 6 
%
% Using eq.(29) in arXiv:1111.4429v3 to calculate the Coulomb-hole operator(Sigma_CH)
% < n k | \Sigma_{CH} (r, r`; 0) | m k > =
% \frac{1}{2} \sum_{q G G`}
% < n k | e^{i (G - G`) \cdot r} | m k >
% [\eps_{G G`}^{-1} (q; 0) - \delta_{G G`}] v (q + G`)
%
  Ech = zeros(nbands,1);
  for n=nv-nvbands+1:nv-nvbands+nbands
    aqsch(:,n) = F*(conj(psir(:,n)).*psir(:,n));
  end

  for i=1:ng
    schx_tmp = zeros(nbands, 1);
    for j=1:ng
      if abs(inveps_1(i,j)) > 1E-6
        G=ksinfo.coulG(i,1:3)-ksinfo.coulG(j,1:3);
        for k=1:ng
          if G==ksinfo.coulG(k,1:3)
                schx(:) = aqsch(k,:) * inveps_1(j,i);
                schx_tmp(:) = schx_tmp(:) + schx(:);
          end
        end
      end
    end
    Ech(:) = Ech(:) + schx_tmp(:) * ksinfo.coulG(i,4);
  end

else
%
% Using eq.(28) in arXiv:1111.4429v3 to calculate the Coulomb-hole operator(Sigma_CH)
% < n k | \Sigma_{CH} (r, r`; 0) | m k > =
% \frac{1}{2} \sum_{n''} \sum_{q G G`} < n k | e^{i (q + G ) \cdot r} | n'' k >
% < n'' k | e^{-i (q + G' ) \cdot r} | m k >
% [\eps_{G G`}^{-1} (q; 0) - \delta_{G G`}] v (q + G`)
%
  Ech = zeros(nbands,nbands);
  for n1 = 1:nbands  
    psin1r = conj(F'*Z(:,n1));
    for n2 = 1:nbands
      psin2r = conj(F'*Z(:,n2));
      for n3 = 1:NBANDS
        psioccr = F'*Z(:,n3);
        psioccn1 = F*(psioccr.*psin1r);
        psioccn2 = F*(psioccr.*psin2r);
        Ech(n1,n2) = Ech(n1,n2) + psioccn1'*(W1*psioccn2);
      end
    end
  end

end 

timeforCH = toc

end %%%if(0)

if(1)

tic

if nargin < 6
  psir = zeros(ng,nbands);
  psir = F'*Z(:,nv-nvbands+1:nbands);
else 
  psir = zeros(ng,NBANDS);
  psir = F'*Z(:,nv-nvbands+1:NBANDS);
end

timeforC2R = toc

chi0 = zeros(ng);
Mg = zeros(ng,nvbands*ncbands);

tic

psivc = prod_states(conj(psir(:,nv-nvbands+1:nv)),psir(:,nv+1:nbands),ev,nv);
Mg = F*psivc;
chi0 = Mg*Mg'/vol;
timeforchi0 = toc

tic

I = eye(ng);
Dcoul = spdiags(ksinfo.coulG(:,4), 0, ng, ng);
inveps = inv(I + Dcoul*chi0);

timeforinveps = toc

Esx_x = zeros(nbands, nbands);
Ech = zeros(nbands, nbands);
Ex = zeros(nbands, nbands);

tic

Dcoul(1,1) = ksinfo.coulG0;
inveps_1 = inveps - I;
W1 = inveps_1*Dcoul;

for iocc = 1:nvbands
  psioccn = F*(psir(:,iocc).*conj(psir(:,1:nbands)));
  Ex = Ex + psioccn'*(Dcoul*psioccn);
  Esx_x = Esx_x + psioccn'*(W1*psioccn);
end

timeforXandSX = toc

tic

for iocc = 1:nbands
  psioccn = F*(psir(:,iocc).*conj(psir(:,1:nbands)));
  Ech = Ech + psioccn'*(W1*psioccn);
end

timeforCH = toc

end %%if(1)

fac = 1.00/vol;

ev = ev * ry2ev
Esx_x = real(diag(0.0 - fac * Esx_x) *ry2ev)
Ech = real(diag(0.5 * fac * Ech * ry2ev))
Ex = real(diag(0.0 - fac * Ex) * ry2ev)
Sig = real(Ex + Esx_x + Ech)
Vxc
eqp = ev - Vxc + Sig
