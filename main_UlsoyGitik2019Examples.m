clc;
close all;
clear all;

flag_table1 = 1;
flag_table2 = 0;
flag_table3 = 0;
flag_table4 = 0;

flag_0Witer1Qiter=1;

%options = [];
options = optimset('TolX',1e-6,'TolFun',1e-6,...
                  'MaxFunEvals',1e5,'MaxIter',1e5,'Display','none');   
%options = optimset('Display','none');   


if(flag_table1==1)
  %Table 1
  disp('==============================');
  disp(' Table 1');
  disp('==============================');
  A  = [0 1; -5 -1];
  Ad = [0 0; -3 -0.6];
  h  = 5;
  kV = [-1;0;1];
  success = calcUlsoyGitik2019Tables(A,Ad,h,kV,options,flag_0Witer1Qiter);
end

if(flag_table2==1)
  %Table 2
  disp('==============================');
  disp(' Table 2');
  disp('==============================');
  
  A  = [0 1; -2.5 2.5];
  Ad = [0 0; 2.5 0];
  h  = 1;
  kV = [-1;0;1];
  success = calcUlsoyGitik2019Tables(A,Ad,h,kV,options,flag_0Witer1Qiter);
end

if(flag_table3==1)
  %Table 3
  disp('==============================');
  disp(' Table 3');
  disp('==============================');
  
  A  = [-27       -0.0097     6;...
          9.5999 -40.2750   -40.6578;...
          0       18.0608     4.1480];
  Ad = [ 0 0 0; ...
        21 0 0;...
         0 0 0];
       
  h  = 0.06;
  kV = [-2;-1;0;1;2];
  success = calcUlsoyGitik2019Tables(A,Ad,h,kV,options,flag_0Witer1Qiter);
end

if(flag_table4==1)
  %Table 3
  disp('==============================');
  disp(' Table 3');
  disp('==============================');
  
  A  = [  0    0     1       0;...
          0    0     0       1;...
        -50    0  -117.8   -67.1;...
          1   -1     0.1    -0.1];
  Ad = [ 0 0   0 0; ...
         0 0   0 0; ...
         0 0 100 0;...
         0 0   0 0];
       
  h  = 0.01;
  kV = [-3;-2;-1;0;1;2;3];
  success = calcUlsoyGitik2019Tables(A,Ad,h,kV,options,flag_0Witer1Qiter);
end

