function [ZZZ] = preProcess(Z)
% preprocess vsdi video.
% input Z should be in numpixels x numframes
    global brn fs params

    switch params.pre.filter % filter according to parameters set by user in csv file
        case 'high'
            [ZZZ,d] = highpass(Z',params.pre.cutoff,fs);
        case 'low'
            [ZZZ,d] = lowpass(Z',params.pre.cutoff,fs);
        case 'band'
            [ZZZ,d] = bandpass(Z',[params.pre.cutoff(1) params.pre.cutoff(2)],fs);
        otherwise 
            ZZZ = Z';
    end
    switch params.pre.normalization % normalize according to parameters set by user in csv file
        case 'z'
            ZZZ = (ZZZ-mean(ZZZ))./std(ZZZ);
        case 'minmax'
            ZZZ = MinMaxNorm(ZZZ);
        otherwise
            ;
    end
    params.pre.filterType = d;
    if params.pre.whiten % whiten data if user specified to do so
        ZZZ = whiten(ZZZ);
    end
    ZZZ = ZZZ';
end

