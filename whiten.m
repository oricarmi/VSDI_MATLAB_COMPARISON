function [Xw] = whiten(X)
% whiten data
Xc = center(X);
[U,S,V] = svd(cov(Xc'));
D = diag(1./diag(sqrt(S)));
whiteM = U*D*U';
Xw = whiteM*Xc;
end

function [Xc] = center(X)
% center data
Xc = X - mean(X);
end

