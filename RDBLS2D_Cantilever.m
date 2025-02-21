function RDBLS2D_Cantilever()
%% Step 1: initialize parameters
[tau,dt,dv,nx,ny,Vmax,E,nu] = deal(1,2,0.9,200,100,0.3,1,0.3);tic;
N = [1 23 23 1; 23 529 529 23; 23 529 529 23;1 23 23 1]/2304;
PHIcoeff = ones(ny+3,nx+3);
phie=conv2(PHIcoeff,N,'valid');
x = phie>0;
V0=mean(x(:));
%% Step 2: prepare for FEM
[KE] = lk(E,nu);
nodenrs = reshape(1:(1+nx)*(1+ny),1+ny,1+nx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)-1,nx*ny,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*(ny+1)+[0 1 2 3] 2 3],nx*ny,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nx*ny,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nx*ny,1);
F = sparse(2*(ny+1)*(nx+1),1);
U = zeros(2*(ny+1)*(nx+1),1);
F(2*(nx*(ny+1)+round(ny/2)+1)) = -1; % Cantilever (right middle point)
freedofs = true(size(F));
freedofs([1:2*(ny+1)]) = false;      % Cantilever (left line)
%% STEP 3: assemble reaction diffusion matrices: NN dNdN
[NNe,dNdNe]=getNNe(1,1);
Bnrs=reshape(1:(nx+3)*(ny+3),ny+3,nx+3);
BVec=reshape(Bnrs(1:ny,1:nx),nx*ny,1);
BMat = zeros(nx*ny,16);
BMat(:,[1:4:end 2:4:end 3:4:end 4:4:end])=repmat(BVec,1,16)+repmat([0:3 ny+3:ny+6  2*ny+6:2*ny+9  3*ny+9:3*ny+12],nx*ny,1);
iN = reshape(kron(BMat,ones(16,1))',256*nx*ny,1);
jN = reshape(kron(BMat,ones(1,16))',256*nx*ny,1);
sNN = reshape(NNe(:)*ones(1,nx*ny),256*nx*ny,1);
NN = sparse(iN,jN,sNN);
sNNdif = reshape(dNdNe(:)*ones(1,nx*ny),256*nx*ny,1);
dNdN = sparse(iN,jN,sNNdif);
T = decomposition(NN + dt*tau*dNdN, 'chol', 'lower');
for iterNum = 1:30
    V0=max(Vmax,V0*dv);
%% Step 4: FEM & sensitivity analysis
    x=max(1e-4,x);
    sK = reshape(KE(:)*x(:)',64*nx*ny,1);
    K = sparse(iK,jK,sK);
    K = (K+K')/2;
    U(freedofs) = K(freedofs,freedofs) \ F(freedofs);
    dc = -x.*reshape(sum((U(edofMat)*KE).*U(edofMat),2),ny,nx);
    compEng(iterNum) = -sum(dc(:));
    disp(['It.: ' num2str(iterNum) ' Compl.: ' num2str(compEng(iterNum))  ' Vol.: ' num2str(mean(x(:))) ' Time:' num2str(toc)]);
    contourf(phie(end:-1:1,:),[0,0],'g','edgecolor','n'); axis equal; axis tight; axis off; drawnow;
    if iterNum>5 && (abs(mean(x(:))-Vmax)<0.005) && all(abs(compEng(end)-compEng(end-5:end-1))< 0.01*abs(compEng(end)))
        return; % Check for convergence
    end
%% Step 5: update the design variables
    dc=[dc(1,1) dc(1,:) dc(1,end); dc(:,1) dc dc(:,end); dc(end,1) dc(end,:) dc(end,end)]; %(ny+2)*(nx+2)
    dc = 0.25*(dc(1:end-1,1:end-1)+dc(2:end,1:end-1)+dc(1:end-1,2:end)+dc(2:end,2:end));%(ny+1)*(nx+1)
    dc=[dc(1,1) dc(1,:) dc(1,end); dc(:,1) dc dc(:,end) ; dc(end,1) dc(end,:) dc(end,end)]; %(ny+3)*(nx+3) %expanded outward one round
    l1 = 0;l2 = 1;
    while (l2-l1>1e-6)
        lambda = (l1+l2)/2; %bi-section Lagrangian method
        dcTest = dc + lambda;%TDN0: -UKU+¦Ë (nx+1)*(ny+1)
        dccoeff = generate_coeff_by_fast_interpolation(dcTest);
        Y = NN*(PHIcoeff(:) - dccoeff(:)*dt);%(ny+3)*(nx+3)        
        PHINew = T\Y;%column matrix
        PHINew = reshape(PHINew,ny+3,nx+3);   %(ny+3)*(nx+3) 
        phie = conv2(PHINew,N,'valid');
        x = phie>=0;
        if mean(x(:)) - V0 > 0 l1 = lambda; else l2 = lambda; end
    end
    PHIcoeff = PHINew/2;
end
function [c] = generate_coeff_by_fast_interpolation(s)
c = zeros(size(s));
for y = 1:size(s,1) c(y, :) = bspline_fast_interpolation(s(y, :)); end
s = c'; c = c';
for y = 1:size(s,1) c(y, :) = bspline_fast_interpolation(s(y, :)); end
c = c';
function [c] = bspline_fast_interpolation(s)
n = length(s)+1;
z1 = 3^0.5-2;
cm = zeros(n+1,1);
cp(1,1)= -1/(1-z1^(2*n))*dot(s,repmat(z1,1,n-1).^(1:n-1).*(1-repmat(z1,1,n-1).^(2*(n-1):-2:2)));
for k = 1:n-1 cp(k+1) = s(k)+z1*cp(k); end
for k = n:-1:2 cm(k) = z1*(cm(k+1)-cp(k)); end
c = 6*cm(2:end-1)';
function [KE] = lk(E,nu)
A11 = [12 3 -6 -3; 3 12 3 0; -6 3 12 -3; -3 0 -3 12];
A12 = [-6 -3 0 3; -3 -6 -3 -6; 0 -3 -6 3; 3 -6 3 -6];
B11 = [-4 3 -2 9; 3 -4 -9 4; -2 -9 -4 -3; 9 4 -3 -4];
B12 = [ 2 -3 4 -9; -3 2 9 -2; 4 9 2 3; -9 -2 3 2];
KE = E/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
function [NNe,dNdNe]=getNNe(ea,eb)
if isempty(which('NNe_func.m')) || isempty(which('dNdNe_func.m'))
    syms u r s a b
    B = [-(u-1)^3/6, u^3/2 - u^2 + 2/3, -u^3/2 + u^2/2 + u/2 + 1/6, u^3/6];
    N = reshape(subs(B,u,s).'*subs(B,u,r),1,[]);
    dN = [diff(N,r)/a; diff(N,s)/b];
    matlabFunction(a*b*int(int(N.'*N,'r','0','1'),'s','0','1'),'file','NNe_func.m');
    matlabFunction(a*b*int(int(dN.'*dN,'r','0','1'),'s','0','1'),'file','dNdNe_func.m');
end
NNe = NNe_func(ea,eb);
dNdNe = dNdNe_func(ea,eb);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This MATLAB code is written by Cong Wang & Shiwei Zhou                   %
% Centre for Innovative Structures and Materials, RMIT University          %
% Please send your comments to s3858767@student.rmit.edu.au                %
%                                                                          %
% The program is proposed for educational purposes.                        %
%                                                                          %
% Disclaimer:                                                              %
% The authors reserve all rights but do not guarantee that the code is free%
% from errors. Furthermore, we shall not be liable in any event.           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%