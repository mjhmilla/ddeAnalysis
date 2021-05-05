clc;
close all;
clear all;


H = rand(5,5);

fprintf('k,\tErr:eig\t\tErr:jor\t\tTime:eig/jor\n');
for k=0:1:10

  t0 = tic;
  Weig    = calcMatrixLambertW(k,H,eps*10);
  teig = toc(tic); 
   
  errWeig = H-Weig*expm(Weig);
  
  t0      = tic;
  Wjor    = lambertw_matrix(k,H);
  tjor    = toc(t0);
  
  errWjor = H-Wjor*expm(Wjor);
  fprintf('%i,\t%1.3e\t%1.3e\t%1.3e\n',...
           k,norm(errWeig),norm(errWjor),...
           teig/tjor);
  
end