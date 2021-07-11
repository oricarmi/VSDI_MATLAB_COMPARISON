function [normData] = MinMaxNorm(data)
% perform min max normalzation
if length(size(data))<3
    temp = data(:);
else
    temp = cat(3,data(:));
end
normData = (temp-min(temp))./(max(temp)-min(temp));
if length(size(data))<3 
    normData = reshape(normData,size(data,1),[]);
else
    normData = reshape(normData,size(data,1),[],size(data,3));
end