function err = calcWIterationError2019(W,A,Ad,h)

err = W*expm(W+h*A)-h*Ad;