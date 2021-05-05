clc;
close all;
clear all;

%%
% From the chatter equation example presented in Yi et al.
%
% Eqn 2:
%
% d/dt x1 = x2
% d/dt x2 = -2*zeta*wn*x2 - (wn^2 + kc/m)*x1 + (kc/m)*x1(t-T)
%           + (kc/(8*fo*m))*( (x1-x1(t-T))^2 - (5/(12*fo))*(x1-x1(t-T))^3 )
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
% After linearizing about the equilbrium point, and substuting km=m*wn^2 
% we are left with:
%
% d/dt [ x1 ] = [               0             1](x1) + [ 0           0](x1(t-T))
% d/dt [ x2 ] = [ -(1 + kc/km)*wn^2  -2*zeta*wn](x2) + [ kc*wn^2/km  0](x2(t-T)).
%
%%

flag_CalcWTransform      = 1;
flag_IterationQ0W1       = 1;

flag_SimRun              = 0;

flag_SimAddInU           = 0;
flag_SimAddInChipForce   = 0;
flag_SimPerturbHistory   = 0;



delta = 1e-4; %Nominal cutting rate per turn.
assert(delta > 0);

%Chatter equation parameters
kckm  =    0.25; %kc/km
wn    =  150.00; 
wn2   =   wn*wn; 
zeta  =    0.05;
T     =    1/50;
fo    =   delta; %Set to the nominal chip thickness at steady state.
                 %(This may not be set correctly.) 
params = [kckm;wn;zeta;T;fo];                  
                                               
A = -[            0,          1;...
     -(1+kckm)*wn2, -2*zeta*wn];
   
Ad = -[       0, 0;...
       kckm*wn2, 0];

if(flag_CalcWTransform==1)     

  tol = eps*100;


  %Test data from Table 1     
  testData(5) = struct('k',zeros(2,1),...
            'S',zeros(size(A,1),size(A,2),2),...
            'lambdaS',zeros(size(A,1),2));

  testData(1).k = [0];
  testData(1).S = [     0,     1;...
                   -33083, -0.24];
                 
  testData(1).lambdaS = [-0.12+181.88i;...
                         -0.12-181.88i];                     

  testData(2).k = [-1];
  testData(2).S(:,:,1) = [             0,          1;...
                           -77988+32093i,  -177-247i];                                
  testData(2).S(:,:,2) = [        0,          1;...
                          -11-1663i,   -92-182i];

  testData(2).lambdaS(:,1) = [ -0.12+181.88i;...
                             -176.73-428.66i];
  testData(2).lambdaS(:,2) = [-91.61;...
                              -0.12-181.88i];

  testData(3).k = [ 1];
  testData(3).S(:,:,1) = [            0,          1;...
                          -77988-32093i,  -177+247i];
  testData(3).S(:,:,2) = [        0,          1;...
                          -11+1663i,   -92+182i];

  testData(3).lambdaS(:,1) = [-0.12-181.88i;...
                            -176.73+428.66i];
  testData(3).lambdaS(:,2) = [       -91.61;...
                              -0.12+181.88i];


  testData(4).k = [-2];
  testData(4).S(:,:,1) = [            0,          1;...
                         -137360+42340i,  -230-570i];
  testData(4).S(:,:,2) = [        0,              1;...
                       77945-31297i,   -177-611i];

  testData(4).lambdaS(:,1) = [ -0.12+181.88i;...
                             -233.30-755.05i];
  testData(4).lambdaS(:,2) = [-0.12-181.88i;...
                            -176.73-428.66i];


  testData(5).k = [ 2];
  testData(5).S(:,:,1) = [            0,          1;...
                         -137360-42340i,  -230+570i];
  testData(5).S(:,:,2) = [        0,              1;...
                       77945+31297i,   -177+611i];

  testData(5).lambdaS(:,1) = [-0.12-181.88i;...
                            -233.30+755.05i];
  testData(5).lambdaS(:,2) = [-0.12+181.88i;...
                            -176.73+428.66i];    


  %%
  % Solve for values of Q and compare it to the data reported in the paper.
  %%

  AdInv = pinv(Ad);
  
  switch flag_IterationQ0W1
    case 0
      disp('Using Q-Iteration');
    case 1
      disp('Using W-Iteration');      
    otherwise assert(0)
      
  end

  fprintf('%s\t%s\t%s\t\t%s\t%s\t%s\n','N:f','N:Q','N:W','P:f','P:Q','P:W');
  for i=1:1:length(testData)

    for j=1:1:size(testData(i).S,3)
      Sk = testData(i).S(:,:,j);
      Wk = (Sk+A).*T;
      Hk = Wk*expm(Wk);
      Qk = -(1/T).*(AdInv*Hk);
      k = testData(i).k;

      w = calcMatrixLambertW(k,Hk,tol);

      errW = w(:,:)-Wk;

      options = optimset('TolX',tol*10,'TolFun',tol*10,...
                        'MaxFunEvals',1e4,'MaxIter',1e4,'Display','none');   
      
      %Outputs from polishing the reported solutions
      fp = [];
      QpErr = [];
      WpErr = [];
      SpErr = [];
      
      %Outputs from solving using a niave initial solution
      fn = [];
      QnErr = [];
      WnErr = [];
      SnErr = [];
                      
      %%
      %Q iteration: polish the solution
      %%
      switch flag_IterationQ0W1
        case 0
          errFcnQ = @(argX)calcQIterationError(argX, k, A, Ad, T, tol);
          [Qp, fp, exitflag] = fsolve(errFcnQ,Qk,options);

          QpErr = norm(Qp-Qk);    
          Hp=-T.*(Ad*Qp);
          Wp = calcMatrixLambertW(k,Hp,tol);    
          WpErr = norm(Wp-Wk);    
          Sp = (1/T).*Wp - A;
          SpErr = Sp-Sk;
          
        case 1
          Wp = Wk;
          errFcnW = @(argX)calcWIterationError(argX,A,Ad,T);
          [Wp, fp, exitflag] = fsolve(errFcnW,Wp,options);

          WpErr = norm(Wp-Wk);    
          Hp = Wp*expm(Wp);
          Qp = (-1/T).*AdInv*Hp;          
          QpErr = norm(Qp-Qk);  
          Sp = (1/T).*Wp - A;
          SpErr = Sp-Sk;
                    
      end
      

      
      %%
      %Q iteration: solve using niave start   
      %%
      switch flag_IterationQ0W1
        case 0
          Qn0   = expm(A.*T);
          errFcnQ = @(argX)calcQIterationError(argX, k, A, Ad, T,tol);       
          [Qn, fn, exitflag] = fsolve(errFcnQ,Qn0,options);
          
          QnErr = norm(Qn-Qk);
          Hn=-T.*(Ad*Qn);
          Wn = calcMatrixLambertW(k,Hn,tol);    
          WnErr = norm(Wn-Wk);       
          Sn = (1/T).*Wn - A;
          SnErr = Sn-Sk;

        case 1
          Qn = expm(A.*T);
          Hn = -(T).*(Ad*Qn);
          Wn = calcMatrixLambertW(k,Hn,tol);
          errFcnW = @(argX)calcWIterationError(argX,A,Ad,T);
          [Wn, fn, exitflag] = fsolve(errFcnW,Wn,options);

          WnErr = norm(Wn-Wk);    
          Hn = Wn*expm(Wn);
          Qn = (-1/T).*AdInv*Hn;          
          QnErr = norm(Qn-Qk);            
          Sn = (1/T).*Wn - A;
          SnErr = Sn-Sk;
      end

      fprintf('%1.1e\t%1.1e\t%1.1e\t\t%1.1e\t%1.1e\t%1.1e\n',...
           norm(fn),QnErr,WnErr,norm(fp),QpErr,WpErr) ;


      here=1;


    end
  end
end

%%
%
% Numerical Simulation
%
%%
  
if(flag_SimRun==1)
  x0 = [0;0];
  if(flag_SimPerturbHistory ==1)
    x0 = [delta;0];
  end

  uFcn = @(argT)calcSmoothStepFunction(argT,0.1,0.2,delta);

  ddeFcn = @(argT,argY,argZ)calcDerivativeYiEqn8(argT,argY,argZ,A,Ad,params, ...
                              uFcn, flag_SimAddInU, flag_SimAddInChipForce);
  ddeHist= @(argT)calcHistoryYiEqn8(argT,params,x0);

  t0    = 0;
  t1    = 5;
  npts  = 100*(t1-t0);

  tspan = [t0:(t1-t0)/(npts-1):t1];

  lags = [T];

  options = ddeset('RelTol',1e-5,'AbsTol',1e-5);
  sol = dde23(ddeFcn,lags,ddeHist,[t0,t1],options);



  %%
  % Sample both solutions
  %%

  dim = size(x0,1);

  x   = zeros(npts,dim);
  xd1 = zeros(npts,dim);
  dx  = zeros(npts,dim);

  for z=1:1:npts
    td = tspan(1,z)-lags(1,1);
    if(td < 0)
      xd1(z,:) = ddeHist(td)';
    else 
      xd1(z,:) = deval(sol,td)';
    end

    x(z,:) = deval(sol,tspan(1,z))';
    dx(z,:)= ddeFcn(tspan(1,z),x(z,:)',xd1(z,:)')';  

  end

  lineTypeExact   = {'-','-'};
  colorExact      = [0.75,0.75,0.75;...
                     0.25,0.25,0.25];
  lineWidthExact  = [1;1];

  set(groot, 'defaultAxesFontSize',10);
  set(groot, 'defaultTextFontSize',10);
  set(groot, 'defaultAxesLabelFontSizeMultiplier',1.2);
  set(groot, 'defaultAxesTitleFontSizeMultiplier',1.2);

  set(groot, 'defaultAxesTickLabelInterpreter','latex');
  set(groot, 'defaultAxesTickLabelInterpreter','latex');
  set(groot, 'defaultLegendInterpreter','latex');
  set(groot, 'defaultTextInterpreter','latex');

  fig_response = figure;

    subplot(4,1,1);
      for i=1:1:dim
        plot(tspan,x(:,i),lineTypeExact{i}, ...
             'Color',colorExact(i,:), 'LineWidth',lineWidthExact(i,1));
        hold on;
      end

      ylabel('$$x(t)$$');
      xlabel('Time (sec)');
      box off;    

      legend('$$x_1$$','$$x_2$$','Location','NorthEast');
      legend boxoff;  

    subplot(4,1,2);
      for i=1:1:dim
        plot(tspan,xd1(:,i),lineTypeExact{i},...
             'Color',colorExact(i,:),'LineWidth',lineWidthExact(i,1));
        hold on;
      end

      xlabel('Time (s)');
      ylabel('$$x(t-T)$$');    
      box off;

      %legend('$$x_1 (t-T)$$','$$x_2 (t-T)$$','Location','NorthEast');
      %legend boxoff;  

    subplot(4,1,3);
      for i=1:1:dim
        plot(tspan,dx(:,i),lineTypeExact{i},'Color',colorExact(i,:),...
            'LineWidth',lineWidthExact(i,1));
        hold on; 
      end
      xlabel('Time (s)');
      ylabel('$$\dot{x}(t)$$');    
      box off; 

      %legend('$$\dot{x}_1 (t-T)$$','$$\dot{x}_2 (t-T)$$','Location','NorthEast');
      %legend boxoff;  

  pageWidth  = 21.;
  pageHeight = 29.7;


  pdfName = 'fig_NumericalDDE_Yi2007.pdf';

  set(fig_response,'Units','centimeters',...
     'PaperUnits','centimeters',...
     'PaperSize',[pageWidth pageHeight],...
     'PaperPositionMode','manual',...
     'PaperPosition',[0 0 pageWidth pageHeight]);     
     %set(findall(figList(i).h,'-property','FontSize'),'FontSize',10);     
     set(fig_response,'renderer','painters');     
     print('-dpdf', pdfName); 
end














