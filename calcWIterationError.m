function err = calcWIterationError(W,A,Ad,T)

err = W*expm(W-A.*T)+Ad.*T;