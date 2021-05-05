function dx = calcDerivativeJunEqn8Exact(t,x,Z,params)
%%
% This function evaluates the exact state derivative of example 3.1 in Sun.
% In this example
%
% dx(t) = a*x(t)+epsilon*x(t)^3 + b*x(t-tau) +f*sin(omega*t)
%
%
% References
%   Sun JQ. A method of continuous time approximation of delayed dynamical 
%   systems. Communications in Nonlinear Science and Numerical Simulation. 
%   2009 Apr 1;14(4):998-1007.
%
%%

a       = params(1);
b       = params(2);
epsilon = params(3);
f       = params(4);
omega   = params(5);
tstep   = params(6);
fftype  = params(7);

xt   = x(1,1);
xt3  = xt*xt*xt;
xdt1 = Z(1,1); %x-delay-no.-1

ffcn = 0;
if(fftype==0)
  ffcn = f*sin(omega*t);  
else
  if(t > tstep)
    ffcn = f;
  end
end

dx = a*xt + epsilon*xt3 +b*xdt1 + ffcn;

