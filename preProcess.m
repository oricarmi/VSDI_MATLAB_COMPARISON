function [ZZZ] = preProcess(Z)
% preprocess vsdi video.
% input Z should be in numpixels x numframes
    global brn fs params
%     ZZ = (Z-min(Z(:)))/(max(Z(:)-min(Z(:))))*255; % make all values [0 255]
%     ZZ = Z; ZZ(ZZ<0) = 0; ZZ(ZZ>1) = 1;
%     ZZ = filtfilt(ones(1,3),3,ZZ'); ZZ = ZZ'; % temporal low pass filtering
%     t = linspace(0,(T-1)./fs,T); 
%     ZZZ = ZZ;
%     for i=1:size(ZZ,1)
%         ZZZ(i) = max(xcorr(Z(40430,:),sin(2*pi*2.*t),'coeff')).*ZZ(i);
%     end  
%     ZZZ = reshape(ZZ,size(brn,1),[],T); % reshape to image frames
%     ZZZ = Z;
%     for i=1:size(ZZZ,3)
%         ZZZ(:,:,i) = medfilt2(ZZZ(:,:,i));
% %         [gradThresh,numIter] = imdiffuseest(ZZZ(:,:,i),'ConductionMethod','quadratic');
% %         ZZZ(:,:,i) = imdiffusefilt(ZZZ(:,:,i),'ConductionMethod','quadratic', ...
% %             'GradientThreshold',gradThresh,'NumberOfIterations',numIter);
%         ZZZ(:,:,i) = imopen(ZZZ(:,:,i),strel('disk',3));
%         ZZZ(:,:,i) = imgaussfilt(ZZZ(:,:,i),3.5);
% %         PSF = fspecial('gaussian',7,2.5);
% %         WT = zeros(size(brn));
% %         WT(5:end-4,5:end-4) = 1;
% %         INITPSF = ones(size(PSF));
% %         ZZZ(:,:,i) = deconvblind(ZZZ(:,:,i),INITPSF);
%     end
    switch params.pre.filter
        case 'high'
            [ZZZ,d] = highpass(Z',params.pre.cutoff,fs);
        case 'low'
            [ZZZ,d] = lowpass(Z',params.pre.cutoff,fs);
        case 'band'
            [ZZZ,d] = bandpass(Z',[params.pre.cutoff(1) params.pre.cutoff(2)],fs);
        otherwise 
            ZZZ = Z';
    end
    switch params.pre.normalization
        case 'z'
            ZZZ = (ZZZ-mean(ZZZ))./std(ZZZ);
        case 'minmax'
            ZZZ = MinMaxNorm(ZZZ);
        otherwise
            ;
    end
    params.pre.filterType = d;
    if params.pre.whiten % whiten data
        ZZZ = whiten(ZZZ);
    end
    ZZZ = ZZZ';
end

