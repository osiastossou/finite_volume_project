function [ CFL, A ] = TransMatrix( Grid, q, VX, VY, VZ, maxdfdS)
%TransMatrix creates a sparse-matrix of the transmissibilities for
%oil/water flow in reservior
%   function uses specified grid (Nx,Ny,Nz) injection/production wells, and
%   edge fluxes to compute a sparse matrix A filled with upwind fractional
%   flow inside reservior. CFL time-step constraint used.

Nx=Grid.Nx;    %import grid
Ny=Grid.Ny;    
Nz=Grid.Nz;

%N-value
N = Nx.*Ny.*Nz;

%Define time/dimension limiters
XP=max(VX,0); XN=min(VX,0); 
YP=max(VY,0); YN=min(VY,0); 
ZP=max(VZ,0); ZN=min(VZ,0); 

%define CFL time step
Vi = XP(1:Nx,:,:)+YP(:,1:Ny,:)+ZP(:,:,1:Nz)-XN(2:Nx+1,:,:)-YN(:,2:Ny+1,:)-ZN(:,:,2:Nz+1);
%All vectors size N
pm = min(Grid.PV./(Vi(:)+max(q,0)));
%multiply with scalar
CFL = pm/maxdfdS;

x1=reshape(XN(1:Nx,:,:),N,1); %Ve- on left x-edges
y1=reshape(YN(:,1:Ny,:),N,1); % Ve- on bottom y-edge
z1=reshape(ZN(:,:,1:Nz),N,1); % Ve- on bottom z edge

x2=reshape(XP(2:Nx+1,:,:),N,1); % Ve+ on right x-edges
y2=reshape(YP(:,2:Ny+1,:),N,1); % Ve+ on top y-edge
z2=reshape(ZP(:,:,2:Nz+1),N,1); % Ve+ on top z-edge

DiagVecs=[z2, y2, x2, min(q,0)+x1-x2+y1-y2+z1-z2, -x1, -y1,-z1]; %6 daigs
DiagIndx=[-Nx*Ny, -Nx, -1, 0, 1, Nx, Nx*Ny]; % location of diags
A=spdiags(DiagVecs, DiagIndx, N, N);
end

