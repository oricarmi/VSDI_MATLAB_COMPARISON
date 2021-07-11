function x = CorrelationParam(orig,rec)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    H = fspecial('laplacian',0);    
    orig2 = conv2(orig,-H,'same');
    rec2 = conv2(rec,-H,'same');
    x = corr2(orig2,rec2);
end

