function success = calcUlsoyGitik2019Tables(A,Ad,h,kV, options, flag_0Witer1Qiter)
success = 0;

tol = max(size(A))*eps*max(1,norm(A));


                
for i=1:1:length(kV)
  Qo = expm(h*A);
  Ho = h*Ad*Qo;
  Wo = calcMatrixLambertW(kV(i),Ho,tol);
  %Wo = lambertw_matrix(kV(i),Ho);
  
  fk = 0;
  Wk = [];
  
  switch flag_0Witer1Qiter
    case 0
      errFcnW = @(arg)calcWIterationError2019(arg,A,Ad,h);
      [Wk, fk, exitflag]=fsolve(errFcnW,Wo,options);
    case 1
      errFcnQ = @(arg)calcQIterationError(arg,kV(i),A,Ad,h,tol);      
      [Qk, fk, exitflag]=fsolve(errFcnQ,Qo,options);
      Hk = h*Ad*Qk;
      Wk = calcMatrixLambertW(kV(i),Hk,tol);      
    otherwise assert(0);
  end
  %assert(exitflag==1);
  Sk = (1/h)*Wk+A;
  d=eig(Sk);
  
  disp('__________');
  fprintf('%i\t\t: k\n',kV(i));
  fprintf('%1.3e\t: Err\n',norm(fk));
  fprintf('%i\t\t: Exit Flag\n',exitflag);
  switch flag_0Witer1Qiter
    case 0
      disp('Winit');
      disp(Wo);
    case 1
      disp('Qinit');
      disp(Qo);
  end
  
  disp('Sk');
  disp(Sk);
  disp('eig(Sk)');
  for j=1:1:size(d,1)
    disp(d(j,1));
  end
  
  
end

success = 1;
