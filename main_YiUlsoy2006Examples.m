clc;
clear all;
close all;

%%
% Worked examples from Yi & Ulsoy (YU) for delay differential equations (DDE's).
% The text below appears in a few papers, and I'm re-typing it here as
% preparation for implementing the examples.
%
% YU looked at DDE's of the form:
%
%   d/dt (x(t)) + A*x(t) + Ad*x(t-T) = 0      t > 0               [1]
%   x(t) = phi(t)                             t in [-T, 0]        
%
% where x(t) is a vector, T is the constant delay, and the matrices A and
% Ad do not commute. Yi & Ulsoy begin by assuming a solution of the form
%
%   x(t) = e^{St} x_0                                             [2]
%
% subtituting into [3] yields
%
%   S e^{St} x_0 + A e^{St} x_0 + Ad e^{S(t-T)} x_0 = 0           [3] 
% 
%   (S + A + Ad e^{-ST} ) e^{St} x_0 = 0                          [4]
%
% Thus in general
%
%   (S + A + Ad e^{-ST} ) = 0                                     [5]
%
% Multiplying by Te^{ST}e^{AT} and re-arrange
% 
%   T(S + A)e^{ST}e^{AT} = - Ad Te^{AT}                           [6]
%
% In general 
%
%  T(S+A)e^{ST}e^{AT} != T(S+A) e^{ (S+A)T }                      [7]
%
% What YU have is close, but not quite in the form of a Lambert W function
% 
%  W(H)e^{W(H)} = H                                               [8]
% 
% YU introduce a matrix Q to bring 6 into a Lambert W equation form
%
%  T(S+A) e^{(S+A)T} = -Ad T Q                                    [9]
%
% Comparing 9 and 8
%
%  (S+A)T = W(-Ad T Q)                                            [10]
%
% Isolating for S
%
%  S = W( -Ad T Q )/T - A                                         [11]
% 
% Now subsituting S into 6
%
%  W( -Ad T Q )e^{W( -Ad T Q )-AT} = - Ad T                       [12]
%
%
% Yi S, Ulsoy AG. Solution of a system of linear delay differential 
% equations using the matrix Lambert function. In 2006 American Control 
% Conference 2006 Jun 14 (pp. 6-pp). IEEE.
%%

tolZeroEigenValue = 1e-12;

%Test: Lambert W transform of a matrix
hTest = rand(2,2);
wH = calcMatrixLambertW(0,hTest,tolZeroEigenValue);
err = norm(wH*expm(wH) - hTest);
assert(err < 1e-6);


%%
%Examples from Yi & Ulsoy 2006
%%

disp('Example: Yi & Ulsoy 2006 Equation 17');
disp('  Warning: I cannot match the results in Table 2');
disp('         : Why? It is not clear. There are several');
disp('           typos in this example, T is not reported');
disp('           and though the example comes from ');
disp('           Lee & Dianat there are sign differences. ');
A = -[ -1, -3; ...
        2, -5];
   
Ad =  -[1.66, -0.697;...
        0.93, -0.330];

T = 1;    
    
%From Table II
testData(3) = struct('k',0,'Q',zeros(size(A)),'S',zeros(size(A)));

testData(3).k =       1;               
testData(3).Q =       [-18.8024 - 10.2243i , 6.0782 - 2.2661i;...
                      -61.1342 - 23.6812i, 1.0161 - 0.2653i];  
testData(3).S =       [-0.3499 + 4.9801i, -1.6253 - 0.1459i;...
                        2.4174 - 0.1308i, -5.1048 + 4.5592i];
                      
testData(3).lambda =  [-1.3990 + 5.0935i;...
                      -4.0558 + 4.4458i ];                   


testData(1).k =       -1;
testData(1).Q =       [ -18.8024 + 10.2243i, 6.0782 + 2.2661i;...
                        -61.1342 + 23.6812i, 1.0161 + 0.2653i];
testData(1).S =       [-0.3499 - 4.9801i, -1.6253 + 0.1459i;...
                        2.4174 + 0.1308i, -5.1048 - 4.5592i];
testData(1).lambda =  [-1.3990 + 5.0935i;...
                      -4.0558 + 4.4458i];            
              

testData(2).k =       0;
testData(2).Q =       [  -9.9183, 14.2985;...
                        -32.7746,  6.5735];                      
testData(2).S =       [ 0.3055,  -0.4150;...
                        2.1317,  -3.3015];
                      
testData(2).lambda =  [-1.0119;...
                      -1.9841];                         
                    
                    
                    
                                                        
%%    
% Now we have to numerically solve for the matrix Q
%
% Q = [q11 q12;
%      q21 q22]
%
% W(-Ad*T*Q)e^{W(-Ad*T*Q)-AT} = -Ad*T
%  
% Where W(-Ad*T*Q) is found by first solving for the eigen values and
% vectors of -Ad*T*Q = V S V^-1, and then applying the Lambert W transform
% to the eigen values 
%
% W(-Ad*T*Q) = V [W_k(S_00) ... 
%                 .
%                 ...       W_k(S_nn)] V^{-1}
%
% The first thing to test is the Lambert W transform of a matrix
%%

disp('Comparing Si and lambdai from Table II of Yi & Ulsoy 2006');

for i=1:1:length(testData)
  
  k = testData(i).k;
  hi  = T.*(-Ad*testData(i).Q);
  wHi = calcMatrixLambertW(k,hi,tolZeroEigenValue);
  Si  = (1/T).*wHi(:,:,1) - A;
  
  serr        = Si        - testData(i).S;
  lambdaErr   = [wLambdaHi(1,1,1);wLambdaHi(2,2,1)]  - testData(i).lambda;
  
  fprintf('  %1.2e\t%1.2e : |Serr|\t|lambdaErr|\n',norm(serr),norm(lambdaErr));
end




%%
% Next we evaluate this constraint equation to see if it is satisfied:
%
% W(-Ad*T*Q)e^{W(-Ad*T*Q)-AT} = -AdT 
%
%%

%Test constraint function

% q0 = [9.9183, 14.2985,-32.7746, 6.5735];
% k  = 0;
% errFcn = @(argQ)calcLambertWTransformError(argQ, k, A, Ad, T);
% err0 = errFcn(q0);

%Solve for the Q matrix for the k^th branch

% options = optimoptions('fminunc','FunctionTolerance',1e-6,'StepTolerance',1e-8);
% [qsol,fval,exitflag] = fsolve(errFcn,q0,options);
% 
% disp('Error');
% disp(fval');

%Display the Q matrix

% Qi  = reshape(qsol,size(A,1),size(A,2)); 
% fprintf('Q-%i\n',k);
% disp(Qi);

%Solve for the S matrix for the k^th branch

% [v,d] = eig( -T.*(Ad*Qi) );
% Wdi = zeros(size(d));
% for i=1:1:size(d,1)
%   Wdi(i,i) = lambertw(k,d(i,i));  
% end
% Whi = v*Wdi/v;
% Si = (1/T)*Whi-A;
% fprintf('S-%i\n',k);
% disp(Si)

%Display the eigenvalues (just for completness)

% fprintf('lambda-%i\n',k);
% disp(d);





    
   