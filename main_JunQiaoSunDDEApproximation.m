clc;
close all;
clear all;

% From Example 3.1 of Sun. The numberical vlaues are listed just after 
% Eqn. 11
%
%   Sun JQ. A method of continuous time approximation of delayed dynamical 
%   systems. Communications in Nonlinear Science and Numerical Simulation. 
%   2009 Apr 1;14(4):998-1007.

%parameters
a           = -0.75;
b           = 0.5;
epsilon     = 0.01;
f           = 0.01;
omega       = 10;
tstep       = 1;
fftype      = 0; % 0: f*sin(omega*t), 1: f*step(t-tstep)
lags        = [0.5];
ctaN        = 10; %number of delay discretizations
params      = [a,b,epsilon,f,omega,tstep,fftype,lags,ctaN];

%%
% Numerically integrate the exact solution
%%

%function handle
ddeFcn = @(argT,argY,argZ)calcDerivativeJunEqn8Exact(argT,argY,argZ,params);
ddeHist= @(argT)calcHistoryJunEqn8(argT,params);

t0   = 0;
t1   = 20;
npts = 100*(t1-t0);

tspan = [t0:((t1-t0)/(npts-1)):t1];

sol = dde23(ddeFcn,lags,ddeHist,[t0,t1]);

%%
% Numerically integrate the CTA solution
%%

x0Cta   = zeros(ctaN+1,1);
dfcnCta = @(argT,argX)calcDerivativeJunEqn8CTA(argT,argX,params);

solnCta = ode45(dfcnCta,[t0,t1],x0Cta);

%%
% Sample both solutions
%%

x   = zeros(npts,1);
xd1 = zeros(npts,1);
dx  = zeros(npts,1);

xCta   = zeros(npts,(ctaN+1));
dxCta   = zeros(npts,(ctaN+1));


for z=1:1:npts
  td = tspan(1,z)-lags(1,1);
  if(td < 0)
    xd1(z,1) = ddeHist(td);
  else 
    xd1(z,1) = deval(sol,td);
  end
    
  x(z,1) = deval(sol,tspan(1,z));
  dx(z,1)= ddeFcn(tspan(1,z),x(z,1),xd1(z,1));  
  
  xCta(z,:) = deval(solnCta,tspan(1,z))';
  
  dxCta(z,:) = dfcnCta(tspan(1,z),xCta(z,:)')';
  
end

lineTypeExact = '-';
colorExact = [1,1,1].*0.75;
lineWidthExact= 2;

lineTypeCta = '-';
colorCta = [0,0,1];
lineWidthCta= 1;


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
    plot(tspan,x,lineTypeExact, 'Color',colorExact,'LineWidth',lineWidthExact);
    hold on;
    plot(tspan,xCta(:,1),lineTypeCta,'Color',colorCta,'LineWidth',lineWidthCta);
    hold on;
    
    ylabel('$$x(t)$$');
    xlabel('Time (sec)');
    box off;    

    titleStr = '';
    eqnStr = '$$ \dot{x}(t) = a\,x(t)+\epsilon\,x^3(t) + b\,x(t-\tau_1)';
    if(fftype == 0)
      titleStr = 'Sin response: ';
      eqnStr = [eqnStr,' + f\,\sin( \omega t) $$'];            
    else
      titleStr = 'Step response: ';
      eqnStr = [eqnStr,' + f\,u(t-\tau_2) $$'];      
    end
    
    title([titleStr,eqnStr]);   
    legend('DDE Numerical Soln','DDE Continuous Time Approx.',...
           'Location','NorthEast');
    legend boxoff;  
    
  subplot(4,1,2);
    plot(tspan,xd1,lineTypeExact,'Color',colorExact,'LineWidth',lineWidthExact);
    hold on;
    plot(tspan,xCta(:,(1+ctaN)),lineTypeCta,'Color',colorCta,'LineWidth',lineWidthCta);
    hold on;
    
    xlabel('Time (s)');
    ylabel('$$x(t-\tau_1)$$');    
    box off;

   

    
  subplot(4,1,3);
    plot(tspan,dx,lineTypeExact,'Color',colorExact,'LineWidth',lineWidthExact);
    hold on; 
    plot(tspan, dxCta(:,1),lineTypeCta,'Color',colorCta,'LineWidth',lineWidthCta);
    xlabel('Time (s)');
    ylabel('$$\dot{x}(t)$$');    
    box off; 
  
    
pageWidth  = 21.;
pageHeight = 29.7;

pdfName = '';

if(fftype == 0)
  pdfName = 'fig_NumericalDDE_vs_ContinousTimeApproxDDE_Sin.pdf';
else
  pdfName = 'fig_NumericalDDE_vs_ContinousTimeApproxDDE_Step.pdf';  
end

set(fig_response,'Units','centimeters',...
   'PaperUnits','centimeters',...
   'PaperSize',[pageWidth pageHeight],...
   'PaperPositionMode','manual',...
   'PaperPosition',[0 0 pageWidth pageHeight]);     
   %set(findall(figList(i).h,'-property','FontSize'),'FontSize',10);     
   set(fig_response,'renderer','painters');     
   print('-dpdf', pdfName); 
    
    