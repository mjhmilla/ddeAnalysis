function dx = calcDerivativeJunEqn8CTA(t,x, params)

a           = params(1,1);
b           = params(1,2);
epsilon     = params(1,3);
f           = params(1,4);
omega       = params(1,5);
tstep       = params(1,6);
fftype      = params(1,7);
tau         = params(1,8);
ctaN        = params(1,9);

dt = tau/ctaN; 
n  = ctaN+1; 

assert((n) == length(x));

%Evalute f
xt = x(1,1);
xdt1= x(n,1);

ffcn = 0;
if(fftype==0)
  ffcn = f*sin(omega*t);  
else
  if(t > tstep)
    ffcn = f;
  end
end

dx = zeros(size(x));
dx(1,1) = a*xt + epsilon*(xt*xt*xt) + b*xdt1 + ffcn;

%Now evaluate the discrete approximations to x(t-tau)

%First difference for the first entry
dx(2,1)                     = (x(1,1)-x(2,1))*(1/dt);

%Central differences for the middle entries
for i= 3:1:(n-1)
  dx(i,1) = (x(i-1,1)-x(i+1,1))./(2*dt);
end

%First difference for the last entry
dx(n,1)                     = (x(n-1,1)-x(n,1))*(1/dt);




