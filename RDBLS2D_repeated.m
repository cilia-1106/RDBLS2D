function RDBLS2D_repeated(tau,dt,dv,nx,ny,Vmax)
%% Step 1: initialize parameters
[E,nu] = deal(1,0.3);
[N,Rx,Ry,iRx,iRy] = transform(nx,ny);
PHIcoeff = ones(ny+3,nx+3);
phie = getphie(nx,ny,PHIcoeff,Rx,Ry,N);
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
F(2,1) = -1;                                      % MBB
freedofs = true(size(F));
freedofs([1:2:2*(ny+1),2*(ny+1)*(nx+1)]) = false; % MBB
%% STEP 3: assemble reaction diffusion matrices: NN dNdN
load NN.mat; % [NNe, NNdife]=getNNe_repeated(1,1);
Bnrs=reshape(1:(nx+3)*(ny+3),ny+3,nx+3);
BVec=reshape(Bnrs(1:ny,1:nx),nx*ny,1);
BMat = zeros(nx*ny,16);
BMat(:,[1:4:end 2:4:end 3:4:end 4:4:end])=repmat(BVec,1,16)+repmat([0:3 ny+3:ny+6  2*ny+6:2*ny+9  3*ny+9:3*ny+12],nx*ny,1);
iN = reshape(kron(BMat,ones(16,1))',256*nx*ny,1);
jN = reshape(kron(BMat,ones(1,16))',256*nx*ny,1);
sNN = allNN(NNe,nx,ny);
NN = sparse(iN,jN,sNN(:));
sNNdif = allNN(NNdife,nx,ny);
dNdN = sparse(iN,jN,sNNdif(:));
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
    disp(['It.: ' num2str(iterNum) ' Compl.: ' num2str(compEng(iterNum))  ' Vol.: ' num2str(mean(x(:)))]);
    contourf(phie(end:-1:1,:),[0,0],'g','edgecolor','n'); axis equal; axis tight; axis off; drawnow;
    if iterNum>5 && (abs(mean(x(:))-Vmax)<0.005) && all(abs(compEng(end)-compEng(end-5:end-1))< 0.01*abs(compEng(end)))
        return; % Check for convergence
    end
%% Step 5: update the design variables
    dc=[dc(1,1) dc(1,:) dc(1,end); dc(:,1) dc dc(:,end); dc(end,1) dc(end,:) dc(end,end)]; %(ny+2)*(nx+2)
    dc = 0.25*(dc(1:end-1,1:end-1)+dc(2:end,1:end-1)+dc(1:end-1,2:end)+dc(2:end,2:end));%(ny+1)*(nx+1)
    l1 = 0;l2 = 1;
    while (l2-l1>1e-6)
        lambda = (l1+l2)/2; %bi-section Lagrangian method
        dcTest = dc + lambda;%TDN0: -UKU+¦Ë (nx+1)*(ny+1)        
        c = [zeros(1,nx+3); [zeros(ny+1,1),generate_coeff_by_fast_interpolation(dcTest),zeros(ny+1,1)]; zeros(1,nx+3)];
        dccoeff = iRy*c*iRx';
        Y = NN*(PHIcoeff(:) - dccoeff(:)*dt);%(ny+3)*(nx+3)
        PHINew = T\Y;%column matrix
        PHINew = reshape(PHINew,ny+3,nx+3);   %(ny+3)*(nx+3) 
        [phie] = getphie(nx,ny,PHINew,Rx,Ry,N);
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
function [N,Rx,Ry,iRx,iRy] = transform(nx,ny)
R = [6 -6 1; 0 1.5 -0.5; 0 0 1];
[Rx,Ry] = deal(eye(nx+3),eye(ny+3));% Transformation matrix Rx & Ry
Rx([1:3,end-2:end],[1:3,end-2:end]) = [R zeros(3);zeros(3) rot90(R,2)];
Ry([1:3,end-2:end],[1:3,end-2:end]) = [R zeros(3);zeros(3) rot90(R,2)];
iRx=inv(Rx);iRy=inv(Ry);
N = [1 23 23 1; 23 529 529 23; 23 529 529 23;1 23 23 1]/2304;
function [sP]=allNN(Pe,nx,ny)
P = cell(ny,nx);
P([1:3,end-2:end],[1:3,end-2:end]) = Pe;% 4 corners
P([1:3,ny-2:ny],4:nx-3) = repmat(Pe(:,3),1,nx-6);% top & bottom
P(4:ny-3,[1:3 nx-2:nx]) = repmat(Pe(3,:),ny-6,1);% left & right
P(4:ny-3,4:nx-3) = repmat(Pe(3,3),ny-6,nx-6); % inside
P = cellfun(@(x) reshape(x, 256, 1), reshape(P,1,nx*ny), 'UniformOutput', false);
sP = cell2mat(P);
function [phie] = getphie(nx,ny,C,Rx,Ry,N)
%% obtain levelset function by basis functions: phie = c1*N1+...+c16*N16
for i = 1:nx
    for j = 1:ny
        CB = (Ry(j:j+3,j:j+3).'*N*Rx(i:i+3,i:i+3)) .* C(j:j+3,i:i+3);%4*4
        phie(j,i) = sum(CB(:));
    end
end

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