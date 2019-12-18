function [ P, VX, VY, VZ, Pvec ] = TPFAedge( Grid, Keff, q, BHP )
%Two point Flux approximation 
%   Summarizes inputs of grid size and lengths into flux approximations
%   between two points using arithmetic/harmonic averages
%
%requires permeability tensor inpput[K] & injection well vector [q]
%interval sizes for each dimension (dx, dy, dz) are calculated from side
%lengths [Lx, Ly, Lz]
%
%outputs final Pressures [in matrix] an Velocities in respective directions


%Define Simulation parameters
Nx=Grid.Nx;    Lx=Grid.Lx;
Ny=Grid.Ny;    Ly=Grid.Ly;
Nz=Grid.Nz;    Lz=Grid.Lz;
K = Keff;
invK = K.^-1; %for transmissibilities

%Initialize edge permeabilities
KedgeX = zeros(Nx+1,Ny,Nz);
KedgeY = zeros(Nx,Ny+1,Nz);
KedgeZ = zeros(Nx,Ny,Nz+1);

%initialize transmissibilities
TX = zeros(Nx+1,Ny,Nz);
TY = zeros(Nx,Ny+1,Nz);
TZ = zeros(Nx, Ny, Nz+1);

%define delta-t
dx = Lx/Nx;
dy = Ly/Ny;
dz = Lz/Nz;

%Arithmetic average of permeabilities on edges
KedgeX(2:Nx,1:Ny,1:Nz) = (K(1,2:Nx,1:Ny,1:Nz)+(K(1,1:Nx-1,1:Ny,1:Nz))/.2);
KedgeY(1:Nx,2:Ny,1:Nz) = ((K(2,1:Nx,2:Ny,1:Nz)+(K(2,1:Nx,1:Ny-1,1:Nz))./2));
KedgeZ(1:Nx,1:Ny,2:Nz) = ((K(3,1:Nx,1:Ny,2:Nz)+(K(3,1:Nx,1:Ny,1:Nz-1)))./2);


%redefine K eqauations for ease of coding
Kright = (invK(1,2:Nx,1:Ny,1:Nz));
Kleft = (invK(1,1:Nx-1,1:Ny,1:Nz));
Kfront = (invK(2,1:Nx,2:Ny,1:Nz));
Kback = invK(2,1:Nx,1:Ny-1,1:Nz);
Kup = invK(3,1:Nx,1:Ny,2:Nz);
Kdown = invK(3,1:Nx,1:Ny,1:Nz-1);

%Transmissibility harmonic average
TX(2:Nx,1:Ny,1:Nz) = ((2*dy*dz)/dx)./((Kright)+(Kleft));
TY(1:Nx,2:Ny,1:Nz) = ((2*dx*dz)/dy)./((Kfront)+(Kback));
TZ(1:Nx,1:Ny,2:Nz) = ((2*dx*dy)/dz)./((Kup)+(Kdown));



%N-value
N = Nx.*Ny.*Nz;

x1 = reshape(TX(1:Nx,:,:),N,1);
x2 = reshape(TX(2:Nx+1,:,:),N,1);
y1 = reshape(TY(:,1:Ny,:),N,1);
y2 = reshape(TY(:,2:Ny+1,:),N,1);
z1 = reshape(TZ(:,:,1:Nz),N,1);
z2 = reshape(TZ(:,:,2:Nz+1),N,1);
DiagVecs = [-z2,-y2,-x2,x1+x2+y1+y2+z1+z2,-x1,-y1,-z1];
DiagIndx = [-Nx*Ny,-Nx,-1,0,1,Nx,Nx*Ny];
A = spdiags(DiagVecs, DiagIndx, N, N);
A(1,1) = A(1,1)+sum(K(:,1,1,1));
full(A)

Pvec = A \ q;

P = reshape(Pvec, Nx,Ny,Nz);

P = P + BHP;

%Define edge flux
VX = zeros(Nx+1,Ny,Nz);
VY = zeros(Nx,Ny+1,Nz);
VZ = zeros(Nx, Ny, Nz+1);

%Compute averaged edge flux
VX(2:Nx,1:Ny,1:Nz) = -(P(2:Nx,1:Ny,1:Nz)-(P(1:Nx-1,1:Ny,1:Nz))).*TX(2:Nx,:,:);
VY(1:Nx,2:Ny,1:Nz) = -((P(1:Nx,2:Ny,1:Nz)-(P(1:Nx,1:Ny-1,1:Nz)))).*TY(:,2:Ny,:);
VZ(1:Nx,1:Ny,2:Nz) = -((P(1:Nx,1:Ny,2:Nz)-(P(1:Nx,1:Ny,1:Nz-1)))).*TZ(:,:,2:Nz);



end

