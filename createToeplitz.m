function [TopMat] = createToeplitz(f0,bandWidth,numHarmonics,omega,T)
global fs
% Create Toeplitz matrix for TSCA function of a stationary quasi-sinusoid 
%   Detailed explanation goes here
dt = 1/fs;
if isempty(numHarmonics)
    numHarmonics = 5;
end
if isempty(omega)   
    omega = ones(1,numHarmonics);
end
if length(omega) ~= numHarmonics
    disp('omega and number of harmonics specified not the same');
    return
end
% [j,k] = meshgrid(0:T-1,0:T-1);
% [j2,k2] = meshgrid(cos(omega*cos(2*pi*
Den = T*sum(omega);
% Nom = sinc(2*pi*bandWidth*(j-k)*dt)*omega*cos(2*pi*[1:numHarmonics]'*f0*(j-k)*dt));
TopMat = zeros(T);
for j=1:T % rows
    for k=1:T % cols
        temp = 0;
        for h = 1:numHarmonics
            temp = temp + omega(h)*cos(2*pi*h*f0*(j-k)*dt);
        end
        TopMat(j,k) = sinc(2*pi*bandWidth*(j-k)*dt)*temp/Den;
    end
end

end

