function x = CorrelationParam(orig,rec)
% calculate correlation parameter 
    H = fspecial('laplacian',0);    
    orig2 = conv2(orig,-H,'same');
    rec2 = conv2(rec,-H,'same');
    x = corr2(orig2,rec2);
end

