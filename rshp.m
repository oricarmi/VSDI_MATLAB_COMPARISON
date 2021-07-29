function [reshapedZ] = rshp(Z) 
% reshape input from 3d matrix to 2d (stack columns vertically, so 3d matrix of size NxMxD
% becomes 2d matrix of (NxM)xD 
global brn brn0
%     if size(Z,1) == size(brn,1)*size(brn,2) && size(Z,2)>1
    try
        reshapedZ = reshape(Z,size(brn,1),size(brn,2),[]); % experimental data
    catch
        reshapedZ = reshape(Z,40,40,[]); % for the simulation
    end
%     elseif
end

