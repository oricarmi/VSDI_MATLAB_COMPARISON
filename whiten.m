function [Xw] = whiten(X)
% whiten data
Xc = center(X);% center the columns
[U,S,V] = svd(cov(Xc')); % perform svd on the covariance of centered data
D = diag(1./diag(sqrt(S))); % make the diagonal matrix of singular vlues
whiteM = U*D*U'; % whitening matrix
Xw = whiteM*Xc;
end

function [Xc] = center(X)
% center data
Xc = X - mean(X);
end

