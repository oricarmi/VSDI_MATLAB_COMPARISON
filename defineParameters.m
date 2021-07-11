function defineParameters(paramsFilePath,what,Z)
% Create parameters of this experiment
    global params
    params = struct('experiment',struct,'pre',struct,'post',struct,'AOF',struct,'GLM',struct,'Tmax',struct,'TSCA',struct,'Nadav',struct);
    params.experiment.T = size(Z,3);
    params.experiment.what = what;
    params.experiment.Z = Z;
    prms = readtable(paramsFilePath);
    if what <10
        params.experiment.N = what; % loc 8 or loc 9 or loc 4
    else % what >10
        params.experiment.N = (what-2)/10; % it is #2 (the 2 indicating it is 2 hz)
    end
    if what/10<1 % 1 hz experiments
        params.experiment.responseSig = load('responseSig.mat').responseSig;
        params.experiment.T1 = 100; % one repetition of stimulus samples
    else % 2 hz experiments
        params.experiment.responseSig = load('responseSig2Hz.mat').responseSig;
        params.experiment.T1 = 50;% one repetition of stimulus asmples
    end
    if ~isempty(prms.numval{strcmp(prms.method,'experiment') & contains(prms.parameter,'circshift')}) || str2num(prms.numval{strcmp(prms.method,'experiment') & contains(prms.parameter,'circshift')})
        params.experiment.responseSig = circshift(params.experiment.responseSig,str2num(prms.numval{strcmp(prms.method,'experiment') & contains(prms.parameter,'circshift')}));
    end
    params.experiment.N = params.experiment.T/params.experiment.T1;
    params.experiment.theoreticalSigs = buildTheoreticalSigs();
    % read parameters from csv file
    params.pre.filter = prms.textval{strcmp(prms.method,'preprocess params') & strcmp(prms.parameter,'filter')};
    params.pre.normalization = prms.textval{strcmp(prms.method,'preprocess params') & strcmp(prms.parameter,'normalization')};
    params.post.normalization = prms.textval{strcmp(prms.method,'postprocess params') & strcmp(prms.parameter,'normalization')};
    params.pre.cutoff = str2num(prms.numval{strcmp(prms.method,'preprocess params') & strcmp(prms.parameter,'cutoff')});
    params.pre.whiten = str2num(prms.numval{strcmp(prms.method,'preprocess params') & strcmp(prms.parameter,'whiten')});
    params.post.gaussfltSTD = str2num(prms.numval{strcmp(prms.method,'postprocess params') & strcmp(prms.parameter,'gaussfltSTD')});
    params.post.medFiltSize = str2num(prms.numval{strcmp(prms.method,'postprocess params') & strcmp(prms.parameter,'medFiltSize')});
    params.AOF.numFramesFrom = str2num(prms.numval{strcmp(prms.method,'AOF params') & contains(prms.parameter,'From')});
    params.AOF.numFramesUntil = str2num(prms.numval{strcmp(prms.method,'AOF params') & contains(prms.parameter,'Until')});
    params.GLM.Noise = str2num(prms.numval{strcmp(prms.method,'GLM params') & strcmp(prms.parameter,'oscillatory noise freqs')});
    params.Tmax.Thresh = str2num(prms.numval{strcmp(prms.method,'Tmax params') & strcmp(prms.parameter,'threshold')});
    if ~isempty(prms.numval{contains(prms.method,'TSCA') & contains(prms.parameter,'freqs')})
        params.TSCA.Noise.freqs = str2num(prms.numval{contains(prms.method,'TSCA') & contains(prms.parameter,'freqs')});
        params.TSCA.Noise.bw = str2num(prms.numval{contains(prms.method,'TSCA') & contains(prms.parameter,'bandwidths')});
        params.TSCA.Noise.numharmonics = str2num(prms.numval{contains(prms.method,'TSCA') & contains(prms.parameter,'harmonics') & ~contains(prms.parameter,'weights')});
        params.TSCA.Noise.harmWeights = str2num(prms.numval{contains(prms.method,'TSCA') & contains(prms.parameter,'weights')});
    end
    params.TSCA.gammas = str2num(prms.numval{contains(prms.method,'TSCA') & contains(prms.parameter,'gammas')});
    params.TSCA.numProj = str2num(prms.numval{contains(prms.method,'TSCA') & strcmp(prms.parameter,'numProjections')});
    params.TSCA.reduceComp = str2num(prms.numval{contains(prms.method,'TSCA') & contains(prms.parameter,'reduceComp')});
    params.Nadav.p = str2num(prms.numval{contains(prms.method,'Nadav') & contains(prms.parameter,'p')});
    params.Nadav.x = str2num(prms.numval{contains(prms.method,'Nadav') & contains(prms.parameter,'x')});
    params.Nadav.t_lmts = str2num(prms.numval{contains(prms.method,'Nadav') & contains(prms.parameter,'limit')});
    params.Nadav.settle = str2num(prms.numval{contains(prms.method,'Nadav') & contains(prms.parameter,'settle')});
end

