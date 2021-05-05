function Q = assembleQ(qVec,qMap,rows,cols)

n = rows*cols;
nq = length(qMap);

qReal = zeros(n,1);
qImag = zeros(n,1);

qReal(qMap,:) = qVec(1:1:nq,:);
qImag(qMap,:) = qVec((nq+1):1:end,:).*1i;

Qr = reshape(qReal,rows,cols);
Qi = reshape(qImag,rows,cols);

Q = Qr+Qi;