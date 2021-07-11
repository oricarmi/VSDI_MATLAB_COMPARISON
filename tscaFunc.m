function [projected,components,D,Alpha,output] = tscaFunc( Z,X,Y,gamma,toProject,reduceComp )
% implementation of tsca algorithm 
% X - signal, Y - noise. Gamma - vector of gammas (first one is signal)
% <---- make all inputs row vectors
if isrow(X.time)
    X.time = X.time';
end
% ---->
% <--- Check that # of gammas is correct
if length(gamma)<length(Y)+length(X)
    disp('not enough gamma inputs');
    return;
elseif length(gamma)>length(Y)+length(X)
    disp('too many gamma inputs');
    return;
end
% --->
T = length(X.time); % P = numel(X.space);
if isvector(X.time)
    Cx = X.time*X.time'/T; %autocorrelation matrix of signal
else
    Cx = X.time;
end
Cy = zeros(size(Cx,1),size(Cx,2),length(Y));
for i=1:length(Y) % autocorrelation matrices of noise
    if isvector(Y(i).time) % if it is a time vector, estimate autocorrelation matrix
        if isrow(Y(i).time)
            Y(i).time = Y(i).time';
        end
        Cy(:,:,i) = Y(i).time*Y(i).time'/T;
    else % it is already autocorrelation matrix
        Cy(:,:,i) = Y(i).time;
    end
end
C = cat(3,Cx,Cy); N = size(C,3);
MAT = zeros(N); Vec = zeros(1,N);
for i=1:N % calculate matrix inner products between the autocorrelation matrices and weights vector  
    for j=1:N
        MAT(i,j) = MatIP(C(:,:,j),C(:,:,i));
    end
    Vec(i) = gamma(i)*trace(C(:,:,i));
end
Alpha = MAT^(-1)*Vec'; % calculate alpha values
Q = zeros(T);
for i=1:length(Alpha) % construct Q matrix
    Q = Q + Alpha(i)*C(:,:,i);
end
if nargin<6 || reduceComp==0 % if no procedure is done to reduce computation time
    M = Z*Q*Z';
    [components,D] = eig(M);
    [~,Ind]=sort(diag(D),'desc');
    components = components(:,Ind);
else % perform svd to reduce computation time
    [U,S,V] = svd(Z,'econ');
    M = S*V'*Q*V*S;
    [W,D] = eig(M);
    [~,Ind]=sort(diag(D),'desc');
	W=W(:,Ind);
    components = U*W;
end
components=components.*repmat(sign(mean(components)),size(components,1),1);
if nargin>4 && toProject>0 % project the data on components
    projected = components(:,1:toProject)'*Z;
    if isvector(X.time)
        components(:,1) = components(:,1).*sign(projected(1,:)*X.time); % fix inversion
    end
else
    projected = [];
end
output = struct('projected',projected,'components',components,'D',D,'Alpha',Alpha,'gammas',gamma,'C',C);
end

