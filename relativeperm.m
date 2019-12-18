function [ Kro, Krw ] = relativeperm(Fluid, Sw)
%relativeperm function for calculation of relative permeabilities
%   function intakes water saturation for 2-phase flow model, and ouputs
%   relative permeabilities of water and subsequent oil phase.
%
%   can be used with higher order functions (increase Nw/No respectively) 
%       Nw = order water   No = order oil
%       sor = residual oil saturation     swc connate water saturation
% Kr0 - corey end-point relative permeability


Sweff =(Sw - Fluid.swc)/((1-Fluid.swc) - Fluid.sor); %effective water saturation
Krw = Fluid.kr0w * (Sweff .^ Fluid.nw); %relative permeability water 
Kro = Fluid.kr0o * ((1 - Sweff).^Fluid.no); %relative permeability oil 

end

