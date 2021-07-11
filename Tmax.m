function [mapT] = Tmax(Z)
% Apply Tmax method. r_thresh: 1. mean + r_thresh*std, gaussfltSTD = std of
% gaussian filter (if zero no filtering)
    global bsln fs sz ump rot fgn brn brn0 frq cmap lgn scl lgn00 fnm c_f vms plt_on pt pc vc xs prms cfn cfn0 basis params
    refff = [params.experiment.responseSig zeros(1,params.experiment.T-params.experiment.T1)];
    r = zeros(size(Z,1),1); tmax = zeros(size(r));
    for i=1:size(Z,1)
        [rtemp,lags] = xcorr(Z(i,:)-mean(Z(i,:)),refff-mean(refff),'normalized');
        [r(i),I] = max(rtemp);
        tmax(i) = lags(I);
    end
    tmaxx = tmax;
    tmaxx(r<(mean(r)+params.Tmax.Thresh*std(r))) = -100;
    A = tmaxx;
    [~,responseMax] = max(params.experiment.responseSig);
    B = [-100, responseMax:params.experiment.T1:params.experiment.T];
    [~,I] = min(abs(bsxfun(@minus,A',B')));
    Anew = B(I)';
    mapT = zeros(size(brn,1),size(brn,2),params.experiment.N);
    for i=1:size(mapT,3)
        temp = zeros(size(Anew))';
        temp(Anew == B(i+1)) = 1;
        temp = temp.*r';
        mapT(:,:,i) = postProcess(reshape(temp,size(brn,1),[]));
    end
%     plotMaps(mapT,what,0);
end

