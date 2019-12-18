function [Mw,Mo,dfw]=Mobilities(s,Fluid)
%Mobilities This function computes the water and oil mobilities using a
% Brooks-Corey dependence on the saturation for both water and oil. It also
% scales the saturation into an effective saturation, which involves the
% residual saturations: connate water and residual oil.
% Inputs: s = water saturation
% Fluid.swc and Fluid.sor are the connate water and residual oil
% saturations, respectively.
% Fluid.nw and Fluid.no are the water and oil powers
% Fluid.vw and Fluid.vo are the water and oil viscosities

S = (s-Fluid.swc)/(1-Fluid.swc-Fluid.sor);           % Rescale saturations
Mw = Fluid.kr0w * S.^Fluid.nw       /Fluid.vw;       % Water mobility
Mo = Fluid.kr0o * (1-S).^Fluid.no   /Fluid.vo;       % Oil mobility

% Note the Matlab notation below:  if (nargout==3)
% This means that ONLY when you 'ask for it', this function will also
% compute the derivatives of the fractional flow function. 
% In other words, if in a script
% or function you write: [mw,mo]=Mobilities(s,Fluid) with two output
% arguments, it will only compute the mobilities, and not the derivatives.
% If on the other hand you write [mw,mo, df]=Mobilities(s,Fluid)
% it will also compute and return the derivatives. For illustrative
% purposes, I intentially used different notation in the above two function
% calls to emphasize once more that it only matters that you give the
% correct NUMBER of inputs and outputs, but not what NAMES you give it in a
% different function or script. A function just takes a number of
% scalar/vector/matrix inputs, and spits out a number of outputs.
if (nargout==3)
  % Derivative below uses chain rule.  
  dMw = Fluid.nw * Fluid.kr0w * S.^(Fluid.nw-1)       /Fluid.vw;
  dMo = Fluid.no * Fluid.kr0o * (1-S).^(Fluid.no-1)   /Fluid.vo;
  %
  % S itself, though, is effective saturation, while derivative is w.r.t. 
  % saturation itself. Chain rule does derivative of S w.r.t. s:
  dMw = dMw/(1-Fluid.swc-Fluid.sor);
  dMo = dMo/(1-Fluid.swc-Fluid.sor);
  %
  % The derivative of the fractional flow function of water
  % w.r.t. water saturation is calculated using, again, the chain rule.
  % Note that you can also use this to calculate the maximum derivative, to 
  % compute CFL condition in IMPES scheme
  dfw=dMw./(Mw + Mo) - Mw./(Mw+Mo).^2.*(dMw+dMo);              % df_w/ds
end
