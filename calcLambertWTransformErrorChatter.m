function err = calcLambertWTransformErrorChatter(args, k,A, Ad, T, params, htpy, scaling, tol, qMap)


assert(length(args)==length(qMap)*2);

if(size(args,1)==1)
  args = args';
end

kckm = params(1,1);
wn   = params(2,1);
zeta = params(3,1);
T    = params(4,1);

Q = assembleQ(args,qMap,size(A,1),size(A,2));

% nq = length(qMap);
% 
% n = size(A,1)*size(A,2);
% qReal = zeros(n,1);
% qImag = zeros(n,1);
% 
% qReal(qMap,:) = args(1:1:nq,:);
% qImag(qMap,:) = args((nq+1):1:end,:);
% 
% Qr = reshape(qReal,size(A,1),size(A,2));
% Qi = reshape(qImag,size(A,1),size(A,2));
% 
% Q = Qr+Qi;

%[vn,dn] = eig( -T.*(Ad*Q) );

d = [Q(1,2)*kckm*wn*wn*T 0;...
     0                   0];
v = [0 -Q(1,2)/Q(1,1);...
     1              1];
vinv = pinv(v);
   
wd = zeros(size(d));
for i=1:1:size(d,1)
  if( abs(d(i,i)) > tol)
    wd(i,i) = lambertw(k,d(i,i));  
  else
    wd(i,i) = lambertw(0,d(i,i));
  end
end

wh = v*wd*vinv;

errM =  wh*expm(wh-A.*T) + Ad.*T;
err = (norm(errM)./scaling).*(1-htpy);
%err = reshape(errM,size(args,1),size(args,2));