function [ P, Vx, Vy, Vz, Pvec ] = MFEedge(Grid, Keff, q )
%MFE This function computes the pressure and velocity fields using
% the mixed-finite-element method (MFE). The routine heavily relies on 
% structured grids, such that the grid size and grid numbering are known
% without even specifying a full grid with grid-coordinates etc.
% Input: 
% Nx, Ny, Nz: number of grid elements in each direction
% hx, hy, hz: domain sizes in x, y, and z directions (widths/heights)
%               deltax, deltay, deltaz will be computed from hx/Nx etc.
% K: tensor permeability field
% q: sink/source injection and production wells.
% Output:
% P: pressure field (constant pressure value for each grid cell
% VcX, VcY, VcZ: velocity field (velocity components constant for each grid
% cell; already averaged for grid cells, rather than grid edges).
%


%Define Simulation parameters
Nx=Grid.Nx;    Lx=Grid.Lx;
Ny=Grid.Ny;    Ly=Grid.Ly;
Nz=Grid.Nz;    Lz=Grid.Lz;
K = Keff;
hx = Lx/Nx;
hy = Ly/Ny;
hz = Lz/Nz;
%
N=Nx*Ny*Nz; 

Ex=(Nx-1)*Ny*Nz;                            % Number of edges in x-direction
Ey=Nx*(Ny-1)*Nz;                            % Number of edges in y-direction
Ez=Nx*Ny*(Nz-1);                            % Number of edges in z-direction
E=Ex+Ey+Ez;                                 % Total number of edges in grid

RHS=zeros(E+N,1);                             % Right-hand side
RHS(E+1:E+N) = q;                             % For each edge E RHS is 
                                              % gravity term. For each grid
                                              % cell after that, RHS is
                                              % sink/source injection
                                              % wells, as provided as input

B=GenB(Nx, Ny, Nz, Lx, Ly, Lz, K);          % Compute B-block of matrix
C=GenC(Nx, Ny, Nz);                         % Compute C-block of matrix
A=[B,C';-C,sparse(N,N)];                    % Assemble sparse matrix
A(E+1,E+1)=A(E+1,E+1)+1;                    % Add a constant to speficy p

warning off;
x=A\RHS;                                    % Solve linear system for
                                            % both pressure and velocities
warning on;

v=x(1:E);                                   % Extract velocities
P=reshape(x(E+1:E+N),Nx,Ny,Nz);             % Extract pressure
%Pvec = P(:);
Pvec = reshape(P, N,1);                       % Reassign P into Vector for plot

Vx = zeros(Nx+1,Ny,Nz);
Vy = zeros(Nx,Ny+1,Nz);
Vz = zeros(Nx,Ny,Nz+1);
Vx(2:Nx,:,:)  = reshape(v(1:Ex),Nx-1,Ny,Nz); 
Vy(:,2:Ny,:)  = reshape(v(Ex+1:E-Ez),Nx,Ny-1,Nz); 
Vz (:,:,2: Nz) = reshape(v(E-Ez+1:E),Nx,Ny,Nz-1);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function B=GenB(Nx, Ny, Nz, Lx, Ly, Lz, K)
%GenB: Generate B matrix for Mixed Finite Element Method
hx = Lx/Nx;
hy = Ly/Ny;
hz = Lz/Nz;
N=Nx*Ny*Nz;

L = K.^(-1);
ex=N-Ny*Nz; ey=N-Nx*Nz; ez=N-Nx*Ny;                      % Number of edges
tx=hx/(6*hy*hz); ty=hy/(6*hx*hz); tz=hz/(6*hx*hy);       % Transmissibilities

X1=zeros(Nx-1,Ny,Nz); X2=zeros(Nx-1,Ny,Nz);              % Preallocate memory
X0=L(1,1:Nx-1,:,:)+L(1,2:Nx,:,:);   x0=2*tx*X0(:);       % Main diagonal
X1(2:Nx-1,:,:)=L(1,2:Nx-1,:,:);     x1=tx*X1(:);         % Upper diagonal
X2(1:Nx-2,:,:)=L(1,2:Nx-1,:,:);     x2=tx*X2(:);         % Lower diagonal

Y1=zeros(Nx,Ny-1,Nz); Y2=zeros(Nx,Ny-1,Nz);              % Preallocate memory
Y0=L(2,:,1:Ny-1,:)+L(2,:,2:Ny,:); y0=2*ty*Y0(:);         % Main diagonal
Y1(:,2:Ny-1,:)=L(2,:,2:Ny-1,:);   y1=ty*Y1(:);           % Upper diagonal
Y2(:,1:Ny-2,:)=L(2,:,2:Ny-1,:);   y2=ty*Y2(:);           % Lower diagonal

Lz1=L(3,:,:,1:Nz-1); z1=tz*Lz1(:);                       % Upper diagonal
Lz2=L(3,:,:,2:Nz);   z2=tz*Lz2(:);                       % Lower diagonal
z0=2*(z1+z2);                                            % Main diagonal

B=[spdiags([x2,x0,x1],[-1,0,1],ex,ex),sparse(ex,ey+ez);...
   sparse(ey,ex),spdiags([y2,y0,y1],[-Nx,0,Nx],ey,ey),sparse(ey,ez);...
   sparse(ez,ex+ey),spdiags([z2,z0,z1],[-Nx*Ny,0,Nx*Ny],ez,ez)];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function C=GenC(Nx, Ny, Nz)
%GenC: Generate C matrix for Mixed Finite Element Method

C=sparse(0,0);                                             % Empty sparse matrix
Nxy=Nx*Ny; N=Nxy*Nz;                                       % Number of grid-points
vx=ones(Nx,1); vy=ones(Nxy,1); vz=ones(N,1);               % Diagonals

for i=1:Ny*Nz                                              % vx-block of C
  Cx=spdiags([vx,-vx],[-1,0]-(i-1)*Nx,N,Nx-1);             % create bidiagonal block
  C=[C,Cx];                                                % append to C
end

for i=1:Nz                                                 % vy-block of C
  Cy=spdiags([vy,-vy],[-Nx,0]-(i-1)*Nxy,N,Nxy-Nx);         % create bidiagonal block
  C=[C,Cy];                                                % append to C
end

C = [C, spdiags([vz,-vz],[-Nxy,0],N,N-Nxy)];               % vz-block of C
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%