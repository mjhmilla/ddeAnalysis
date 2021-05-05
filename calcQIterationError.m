function err = calcQIterationError(Q, k,A, Ad, T, tol)

H = -T.*(Ad*Q);
W = calcMatrixLambertW(k,H, tol);

err =  W*expm(W-A.*T) + Ad.*T;
