function [  ] = TwoPhaseSimulatorImplicit( Grid, Fluid, q, BHP, tmax, nt, cflx, timemethod, method) 
%2-phase flow simulator for Reservior Modeling
%   Simulator using Implicit Pressure/ Explicit Saturation calculations to
%   quatify reservior scenarios


%SI Units (because numbers matter)
%Unit Conversions
year = 60*60*24*365;
q= (q *(Grid.PVtot))/year;                                                 %(flow rate m^3/s)
milliDarcy = 9.869233e-16;                                                 %permeability
Grid.perm = Grid.perm * milliDarcy;
Pas = 1e-3;                                                                %Pascal Second (Petroleum industry Poise/[P])
Fluid.vw = Fluid.vw*Pas;
Fluid.vo = Fluid.vo*Pas;
Pascal = 100000;                                                           %Pressure conversion to metric
BHP = BHP*Pascal;
tmax = (tmax *Grid.PVtot)/sum(max(q,0));



Sw = linspace(Fluid.swc,1-Fluid.sor);
[mobw, mobo, mtot] = FracFlow(Fluid, Sw); 
mtot = mobo + mobw;
fracflow = mobw./mtot;

dfdS = (fracflow(2:end)-fracflow(1:(end-1)))./(Sw(2:end)-Sw(1:(end-1)));
maxdfdS = max(dfdS);
dtbig = tmax/nt;
Sw =Fluid.swc*ones(Grid.N,1);                                              %water saturation 

Pc=[0; 1; 0]; Tt=0;                                                        %plotting time sensitive variables
myfig =figure('Color','White','Position',[0 500 1600 500]);                %create plotting figure

for i = 1:nt
        [mobw,mobo,mtot] = FracFlow(Fluid, Sw);
        Keff = reshape([mtot,mtot,mtot]',3,Grid.Nx,Grid.Ny,Grid.Nz).*Grid.perm;
     if method == 1
         tic;
        [P,Flux.VX,Flux.VY,Flux.VZ] = TPFAedge(Grid, Keff, q, BHP);                       %Two point flux approx method (explicit)
     else
         tic
        [P,Flux.VX,Flux.VY,Flux.VZ] = MFEedge(Grid, Keff, q);                             %Mixed Finite element method (explicit)
     end   
        [CFL, A] = TransMatrix(Grid,q,Flux.VX,Flux.VY,Flux.VZ,maxdfdS);                   %Matrix A (Permeability matrix) & CFL minimum time-step 
        Nts = ceil(dtbig/CFL); dfw = dtbig/Nts./Grid.PV;
        if(timemethod~=1)                                                  %implicit method flag 
          Nts =cflx;
          St = dtbig/Nts;
        end
            for tsmall = 1:Nts                                             %small time step for Saturation (not needed for implicit)
                if(timemethod~=1)
                 Sw = ImplicitBackwardEuler(Grid,Sw,Fluid,q,St,A);         %implicit saturation
                else
                [mobw,mobo] = FracFlow(Fluid, Sw);
                fw = mobw./(mobw + mobo);
                Sw = Sw + dfw.*(A*fw + max(q,0));                          %explicit saturation
                end
            end
            toc;
        time =(i)*(dtbig);                                                 %time varaible for plotting
        [Tt, Pc] = TwoPhasePlot( Grid, Fluid, q, Sw, time, P, Flux, Tt, Pc, tmax);  %plot data
end

