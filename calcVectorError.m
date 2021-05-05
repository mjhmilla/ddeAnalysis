function err = calcVectorError(arg)

err = [sin(arg(1,1))-0.5;...
       sin(arg(2,1))-0.75];