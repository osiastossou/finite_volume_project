%%% 2-phase Flow modeling Scrip
%%% james ross 11/16/2017 
%%%Reservior Modeling ES 5751 

close all;
clear all;
warning('off', 'Octave:possible-matlab-short-circuit-operator');

%Define Simulation parameters

%Grid Properties
Grid.Nx =4;     Grid.Lx = 1;%meters     
Grid.Ny = 4;    Grid.Ly = 1;%meters   
Grid.Nz = 1;     Grid.Lz = 1;%meters    

% %SPE10
% feet = 0.3048;
% Grid.Lx=Grid.Lx * feet; Grid.Ly=Grid.Ly * feet; Grid.Lz=Grid.Lz *feet;
% Grid.Nz = 10; 
% Grid.Nz = 10; 

Grid.dx = Grid.Lx/Grid.Nx; %grid cell sizes
Grid.dy = Grid.Ly/Grid.Ny;
Grid.dz = Grid.Lz/Grid.Nz;
Grid.N = Grid.Nx.*Grid.Ny.*Grid.Nz;

%Rock Properties
Grid.perm = 100*ones(3,Grid.Nx,Grid.Ny,Grid.Nz); %reservior permeability
Grid.por = 0.3*ones(Grid.Nx,Grid.Ny,Grid.Nz); %reservior porosity
Grid.por = max(Grid.por, 0.01); %eliminate small pororsity to keep time-steps viable
Grid.V = (Grid.dx*Grid.dy*Grid.dz);  %reservior volume
Grid.PV = Grid.V*Grid.por(:);        %reservior pore volume
Grid.PVtot = sum(Grid.PV);


% %SPE10
% load spe10;
% Layer = 1;
% Grid.perm=KU(:,1:Grid.Nx,1:Grid.Ny,Layer:Layer+Grid.Nz-1);

%Injection/Procuction Wells
q = zeros(Grid.N,1);     q(1)=.15;     q(Grid.N)=-.15;
% Bore-hole pressure
BHP  =  300 ;%bar

%Fluid Properties
Fluid.vw = 0.4; % Water viscosity
Fluid.vo = 0.1; %Oil Viscosity
Fluid.nw = 1; % power Corey water rel perm
Fluid.no = 1; % power Corey oil rel perm
Fluid.kr0w = 1; %Corey end-point rel perm for water
Fluid.kr0o = 1; %Corey end-point rel perm for oil
Fluid.sor = 0; %Residual oil saturation 
Fluid.swc = 0; % Connate water saturation

nt = 30;        %number of time steps
tmax = 1; %max time(average 70% injected pore volume @which H20 is produced)


% Solver flags
cflx= 10; %number of implicit time steps per pressure time step
timemethod  =  2; %flag for Implicit method
method = 1; %TPFA =1, MFE =2

TwoPhaseSimulatorImplicit( Grid, Fluid, q, BHP, tmax, nt, cflx, timemethod, method);