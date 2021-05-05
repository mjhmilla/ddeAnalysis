% calculate the nth deritive of lambert W function
% Shiming Duan, Aug.05, 2002
function z = der_lambertw(branch,x,order)
syms s
y=lambertw(branch,s);
for i=1:order
   y=diff(y);
end
z = subs(y, x);
