function [wH] = calcMatrixLambertW(k,H,toleranceBranching)

[v,d] = eig( H );
for j=1:1:size(d,1)
  if(abs(d(j,j))>toleranceBranching)  
    d(j,j) = lambertw( k, d(j,j));    
  else
    d(j,j) = lambertw( 0, d(j,j));        
  end
end
wH = v*d*pinv(v);

