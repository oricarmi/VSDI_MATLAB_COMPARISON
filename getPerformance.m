function [MSE,PSNR,CNR,mssim,Corr,CP] = getPerformance(rec,orig,ROI,offROI)
% Calculate 6 performance measures of image reconstruction
% input should be an image m by n (not vector)
% ROI and offROI are same size as image (m by n), logical matrices)
    if nargin<4
        [orig,ROI,offROI] = createOrig(rec);
    end
    MSE = sum((orig(:)-rec(:)).^2);
    PSNR = psnr(rec,orig);
    mssim = MSSIM(orig,rec,4);
    CNR = calcCNR(rec,ROI,offROI);
    Corr = corr2(orig,rec);
    CP = CorrelationParam(orig,rec);
end

function [orig,ROI,offROI] = createOrig(rec)
    global brn
        img2fit = MinMaxNorm(rec); img2fit(img2fit<prctile(img2fit(:),85)) = 0;
        [fitresult, zfit, fiterr1, zerr, resnorm, rr] = fmgaussfit(1:size(brn,2),1:size(brn,1),img2fit); % fit gaussian to image
        hFig = figure; subplot(1,2,1); imagesc(rec); subplot(1,2,2);imagesc(zfit); % plot image and fitted gaussian
        isFitGood = input('is fit good? [0/1]'); % prompt user if fit is good 
        close(hFig);
        while ~isFitGood % while fit isn't good, keep asking for new ROI until it is 
            R = roipoly(MinMaxNorm(rec)); % select ROI
            img2fit = double(R); 
            [fitresult, zfit, fiterr1, zerr, resnorm, rr] = fmgaussfit(1:size(brn,2),1:size(brn,1),img2fit); % fit image after ROI selction
            hFig = figure; subplot(1,2,1); imagesc(rec); subplot(1,2,2);imagesc(zfit);
            isFitGood = input('is fit good? [0/1]');% prompt user if fit is good 
            close(hFig);
        end
        orig = zfit; ROI = zfit>=prctile(zfit(:),99); offROI = ~ROI;
end