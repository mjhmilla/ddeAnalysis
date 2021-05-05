function dx = calcDerivativeYiEqn8(t,x,Z,A,Ad,params, inputFunction, ...
                                  flag_addinU, flag_addInChipForce)
%%
% From the chatter equation example presented in Yi et al.
%
% Eqn 2:
%
% d/dt x1 = x2
% d/dt x2 = -2*zeta*wn*x2 - (wn^2 + kc/m)*x1 + (kc/m)*x1(t-T)
%           + (kc/(8*fo*m))*( (x1-x1(t-T))^2 - (5/(12*fo))*(x1-x1(t-T))^3 )
%
% where
%
%   x1   : tool edge position (+ is into the work piece)
%   x2   : tool edge velocity
%   T    : delay which is 2*pi/Omega, where Omega is the spindle speed.
%   m    : mass of the cutting tool
%   zeta : the damping factor of the spring-mass-damper of the cutter 
%   omega: the natural frequency of the spring-mass-damper of the cutter
%   kc   : the cutting coefficient derived from a stationary cutting force
%          model which is a function of chip speed, chip thickness 
%   f    : chip thickness
%   fo   : chip thickness at steady state
%
% at equilibrium d/dt x1 = d/dt x2 = 0. Assuming no vibration is left from
% the previous cycle then x1 = x2 = x1(t-T) = 0
%
% d/dt x1  = x2
% d/dt x2  =  -2*zeta*wn*x2 - (wn^2 + kc/m)*x1 - (kc/m)*x1(t-T)
%
% which can be put into block form
%
% d/dt [ x1 ] = [             0           1 ](x1) + [ 0     0](x1(t-T))
% d/dt [ x2 ] = [ -(wn^2 + kc/m)  -2*zeta*wn](x2) + [ kc/m  0](x2(t-T)).
%
% Yi et al. replace km = wn^2*m, where km is structural stiffness, to yield
%
% d/dt [ x1 ] = [               0             1](x1) + [ 0           0](x1(t-T))
% d/dt [ x2 ] = [ -(1 + kc/km)*wn^2  -2*zeta*wn](x2) + [ kc*wn^2/km  0](x2(t-T)).
%
%
% @param t time
% @param x state
% @param Z delayed state
% @param params
% @return the state derivative of the time-delayed system.
%%


kckm = params(1,1);
wn   = params(2,1);
zeta = params(3,1);
T    = params(4,1);
fo   = params(5,1);

wn2  = wn*wn;

x1d  = x(1,1)-Z(1,1);    
x1d2 = x1d*x1d;    
x1d3 = x1d2*x1d;

u = zeros(size(Z));
if(flag_addinU==1)
  u = inputFunction(t);
end
f = zeros(size(Z));
if(flag_addInChipForce==1)
  f = [0; (wn2*kckm/(8*fo))*( (x1d2) - (5/(12*fo))*x1d3 )];    
end   
dx = -A*x + -Ad*(Z+u) + f;    

