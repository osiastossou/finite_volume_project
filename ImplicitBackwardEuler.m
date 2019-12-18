function S=ImplicitBackwardEuler(Grid,S,Fluid,q,T,A);
%NetRaph This function updates the transport equation for saturations
% Using the implicit backward Euler method. The implicit system of 
% equations is solved using the iterative Newton-Raphson method.
% Inputs are the grid and fluid properties, the saturations from the 
% previous timestep, the wells, an initial time-step, and the pre-computed
% matrix A. Output is the new saturation.

N = Grid.Nx*Grid.Ny*Grid.Nz;                                 % number of unknowns

conv=0; IT=0; S00=S;
while conv==0;
 % The time-step (below) is initially provided as input
 % for IT=0, dt=T
 % However, when the Newton method does not converge, because 
 % the initial guess was not good enough, the time-step is split
 % in half (for IT=1). If it still doesn't converge, it is split
 % in half again, until the Newton method converges to an accurate
 % saturation solution.
 dt = T/2^IT;                                                % timestep
 dtx = dt./(Grid.V(:)*Grid.por(:));                          % timestep / pore volume
 % Note, the trick below is not stricktly necessary, 
 % but convenient: the two terms in the transport equation
 % are already multiplied by dt/PV, written in sparse matrix notation.
 fi = max(q,0).*dtx;                                         % injection
 B=spdiags(dtx,0,N,N)*A;

 %testB = full(B)
 I=0;
 % The loop below is only necessary when Newton does not converge the first
 % time, and the time-step needs to be split into two (and maybe into two
 % again).
 while I<2^IT;                                               % loop over sub-timesteps
     % S0 saves the initial input saturation, so we can use it if we have
     % to start over again with a smaller time-step.
  S0=S; dsn=1; it=0; I=I+1;

  % dsn is the TOLERANCE; We require the error in the iterative Newton
  % approximation to be less than 0.001. This is actually not a very
  % 'tight' tolerance (i.e. it is still a rather large error), but we use
  % this to speed up the computation. A smaller error would require more
  % Newton steps.
  while dsn>1e-3 & it<10;                     % Newton-Raphson iteration
      % We need the mobilities, and the derivative of mobilities
      % You may need to modify your relperm/mobility function:
   [Mw,Mo,df]=Mobilities(S,Fluid);    % mobilities and derivatives
   % The line below evaulates the derivative of our transport equation wrt
   % saturations, as given in the lecture notes. Written in sparse matrix
   % nodation: speye is a Matlab function for the unit matrix of size N
   % with all zeros, except ones on the main diagonal. spdiags(df,0,N,N) is
   % notational we used before: it puts the vector df on the main diagonal
   % of a N by N matrix (as given by the zero in the function call) with
   % zeroos everywhere else.
   dG=speye(N)-B*spdiags(df,0,N,N);                          % G'(S)   
   fw = Mw./(Mw+Mo);                                         % fractional flow
   % This is the transport equation itself, where S0 is the saturation from
   % the previous time-step, and S is our current guess/approximation to
   % the new saturation, which will be improved in each iteration of this
   % loop.
   G = S-S0-(B*fw+fi);                                       % G(s)
   ds = -dG\G;                                               % increment ds
   S = S+ds;                                                % update S
   % dsn is the difference between the previous and current
   % guess/approximation to the new saturation. The norm gives the size/'length' of
   % this vector (saturations in each grid cell). If it is very small, then
   % we are very close to the true solution and have a good approximation
   % for the new saturation.
   dsn = norm(ds);                                           % norm of increment
   it = it+1;                                                % number of N-R iterations
  end
    % Below: if we did 10 Newton iterations, but the error in the new
    % saturations is still more than 0.001, this means we didn't converge
    % yet because the initial guess was not good enough. In this case, we
    % need to reduce the time-step and try again.
  if dsn>1e-3; I=2^IT; S=S00; end                            % check for convergence
 end
    % If the error is small, we have an accurate new saturation and we're
    % done.
 if dsn<1e-3; conv=1;                                        % check for convergence
 else
     IT=IT+1;
 end                                           % if not converged, decrease
end                                            %   timestep by factor 2

