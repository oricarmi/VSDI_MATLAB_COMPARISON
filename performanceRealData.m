function [TSCAwGLM,TSCAnoGLM,Tmax,AOF,Corr,GLM,Nadav] = performanceRealData(result)
% Get the performance measures of each individual map and combined
% retinotopic map
global params
    fn = fieldnames(result);
    performances = cell(1,length(fn)); 
    for i=1:length(fn) % iterate the methods 
        if contains(fn{i},'cluster')
            continue
        end
        thisPerformance = zeros(params.experiment.N,6); % number of maps x 6 performance measures
        for j=1:params.experiment.N % iterate the maps in each method
            [MSE,PSNR,CNR,mssim,Corr,CP] = getPerformance(result.(fn{i}).maps(:,:,j),params.experiment.optimalMaps.orig(:,:,j),params.experiment.optimalMaps.ROI(:,:,j),params.experiment.optimalMaps.offROI(:,:,j));
            thisPerformance(j,:) = [MSE,PSNR,CNR,mssim,Corr,CP]; 
        end
        performances{i} = thisPerformance;
    end
    TSCAwGLM = performances{1};
    TSCAnoGLM = performances{2};
    Tmax = performances{3};
    AOF = performances{4};
    Corr = performances{5};
    GLM = performances{6};
    Nadav = performances{7};
end

function [orig,ROI,offROI] = createOrig(rec)
    global brn
    img2fit = MinMaxNorm(rec); img2fit(img2fit<prctile(img2fit(:),85)) = 0;
    [fitresult, zfit, fiterr1, zerr, resnorm, rr] = fmgaussfit(1:size(brn,2),1:size(brn,1),img2fit);
    hFig = figure; subplot(1,2,1); imagesc(rec); subplot(1,2,2);imagesc(zfit);
    isFitGood = input('is fit good? [0/1]');
    close(hFig);
    while ~isFitGood % while fit isn't good, keep asking for new ROI until it is
        again = 1; R = zeros(size(brn));
        while again
            thisR = roipoly(MinMaxNorm(rec));
            R = R | thisR;
            again = input('draw another ROI? [0/1]');
        end
        img2fit = double(R);
        [fitresult, zfit, fiterr1, zerr, resnorm, rr] = fmgaussfit(1:size(brn,2),1:size(brn,1),img2fit);
        hFig = figure; subplot(1,2,1); imagesc(rec); subplot(1,2,2);imagesc(zfit);
        isFitGood = input('is fit good? [0/1]');
        close(hFig);
    end
    orig = zfit; ROI = zfit>=prctile(zfit(:),99); offROI = ~ROI;
end