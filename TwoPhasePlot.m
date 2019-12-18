function [Tt, Pc] = TwoPhasePlot( Grid, Fluid, q, Sw, time, P, Flux, Tt, Pc, tmax )
%Plotting function for 2-phase flow simulator
%   %James Ross 10/26/17


[max_value, index] = min(q(:));
[mobw,mobo] = FracFlow(Fluid,Sw(index)); mtot = mobw+mobo;

%Recovered Fraction
IOIP = sum((1-Fluid.swc).*Grid.PV);
OIP = sum((1-Sw) .* Grid.PV);
RF = 1-(OIP ./ IOIP);

Tt = [Tt, time];
Pc = [Pc,[mobw/mtot; mobo/mtot; RF]];

%Unit Conversions
milliDarcy = 9.869233e-16 ;
Grid.perm = Grid.perm/milliDarcy;
year = 60*60*24*365;
Timeyear = time./year;

%PLOTTING
%change sliceat variables for different visualization levels
sliceatX = ceil(Grid.Nx/2);
sliceatY = ceil(Grid.Ny/2);
sliceatZ = ceil(Grid.Nz/2);

%define X, Y for plot axis distance measurements
X = linspace(0,Grid.Nx  * ((Grid.Lx/Grid.Nx)),Grid.Nx);
Y = linspace(0,Grid.Ny  * ((Grid.Ly/Grid.Ny)),Grid.Ny);
Z = linspace(0,Grid.Nz  * ((Grid.Lz/Grid.Nz)),Grid.Nz);

%More visualization variables
lxlylz = [Grid.Lx, Grid.Ly, Grid.Lz];
nxnynz = [Grid.Nx, Grid.Ny, Grid.Nz];
%grid3D = cartGrid(nxnynz,lxlylz);



%subplot('position',[0.07 .12 .6 .8]);% Make left subplot
subplot(2,2,2)
plot(Tt/year,Pc(1,:),Tt/year,Pc(2,:),Tt/year,Pc(3,:)); %plot functions
axis([0, tmax/year, -0.05, 1.05]);
sattitle = strcat('Saturation Oil & Water Fluid Recovered with Total Recovered Oil Fraction @ t=',num2str(Timeyear),'yrs');
title(sattitle, 'FontName','Arial', 'FontSize', 14);
xlabel('Time(years)');ylabel('Ratio');
legend('Water','Oil','RF');



%subplot('position',[0.7 .12 .3 .8]);% Make right subplot
subplot(2,2,4)
Sgrid = reshape((1-Sw),Grid.Nx, Grid.Ny, Grid.Nz);
        sliceat = 1;
        contourf(linspace(Grid.dx/2, Grid.Lx-Grid.dx/2,Grid.Nx),linspace(Grid.dy/2,Grid.Ly-Grid.dy/2,Grid.Ny),permute(Sgrid(:,:,sliceat),[2,1,3]),11,'k');
        axis equal; caxis([Fluid.swc (1-Fluid.sor)]);
        title('Oil Pore Space Saturation (@(Depth Total/2)[H20 injection/extraction]');
        xlabel('x-dist from injection well(m)');ylabel('y-dist from injection well(m)');
        h = colorbar; ylabel(h, 'Ratio Oil Remaining')
        
        
%surface plot
subplot(2,2,3);% Make pressure subplot
    surf(X, Y, permute(P(:,:,sliceatZ),[2,1,3]));
    hold on;
    title('Pressure Surfer plot');xlabel('x(m)');ylabel('y(m)');zlabel('Pressure(pascal)');
    hh = colorbar; ylabel(hh, 'Pressure(pascal)')
    %view(3);
    hold off;
    
%P/V Quiver plot
subplot(2,2,1)
contourf(X, Y, permute(P(:,:,sliceatZ), [2,1,3]));
hold on;
title('Pressure gradient in x-y');
xlabel('x(m)');ylabel('y(m)');
hhh = colorbar; ylabel(hhh, 'Pressure(pascal)')
%quiver(X, Y, Flux.VX(:,:,sliceatZ)',Flux.VY(:,:,sliceatZ))';
axis equal
hold off;

% %figure
% %plotCellData(grid3D, Pvec);
% 
% figure
% title('Cell Pressure plot');xlabel('x(m)');ylabel('y(m)');colorbar;
% plotGrid(grid3D, 'FaceAlpha', 0, 'EdgeAlpha', 1);
% plotCellData(grid3D , Pvec , find(Pvec>mean(Pvec)));

drawnow;
end

