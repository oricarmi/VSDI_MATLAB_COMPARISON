function [CNR] = calcCNR(rec,ROI,offROI)
% Calculate CNR according to pixels of roi and background (offroi)
    t = rec(logical(ROI)); r = rec(logical(offROI));
    mu_t = mean(t(:)); mu_r = mean(r(:)); 
    if mu_t<mu_r
        CNR = -100;
        return
    end
    var_t = var(t(:)); var_r = var(r(:));
    CNR = 10*log10((mu_t-mu_r)/(sqrt(var_t+var_r)));
end

