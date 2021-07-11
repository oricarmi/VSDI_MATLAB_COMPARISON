function [reshapedZ] = rshp(Z)
global brn brn0
%     if size(Z,1) == size(brn,1)*size(brn,2) && size(Z,2)>1
    try
        reshapedZ = reshape(Z,size(brn,1),size(brn,2),[]); % experimental data
    catch
        reshapedZ = reshape(Z,40,40,[]); % for the simulation
    end
%     elseif
end

