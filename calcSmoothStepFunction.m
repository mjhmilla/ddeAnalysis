function y = calcSmoothStepFunction(t, tstart, ton, stepMagnitude)

y = [0;0];
if(t >= tstart && t <= ton)
    pos = stepMagnitude*(0.5-0.5*cos( pi*(t-tstart)/(ton-tstart) ));
    dArg = pi*(1)/(ton-tstart);
    vel = stepMagnitude*(0.5*sin( pi*(t-tstart)/(ton-tstart) ))*dArg;
    y = [pos;vel];
elseif(t > ton)
    y=[stepMagnitude;0];
end