function [ mobw, mobo, mtot ] = FracFlow(Fluid, Sw)
%FracFlow function: computed fractional flow for oil/water through reservior 
%   function values vary between 0-1, proportional to which part is
%   water/oil respectively
[ Kro, Krw ] = relativeperm(Fluid, Sw);                                    %need effective permeability Oil & water
mobw = Krw./Fluid.vw;                                                      %mobility water
mobo = Kro./Fluid.vo;                                                      %mobility oil
mtot = mobw  + mobo;                                                       %mobility total (always 1[for 2 species])
end

