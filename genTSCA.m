function [mapTSCA] = genTSCA(Z,ZZ)
% Generate TSCA
    Z = MinMaxNorm(Z); % perform minmax normalization prior (no GLM)
    ZZ = MinMaxNorm(ZZ);% perform minmax normalization prior (after GLM)
global fs basis params brn cfn
    noise1.time = eye(params.experiment.T)/params.experiment.T; % autocorrelation matrix of white noise
    if ~isempty(params.TSCA.Noise.freqs) % build oscillatory noises matrices
        for i=1:length(params.TSCA.Noise.freqs)
            noise2(i).time = createToeplitz(params.TSCA.Noise.freqs(i),params.TSCA.Noise.bw,params.TSCA.Noise.numharmonics,params.TSCA.Noise.harmWeights,params.experiment.T);
        end
    else
        noise2 = [];
    end
    mapTSCA1 = cell(params.experiment.N,1); % toepliz denoising
    for i=1:params.experiment.N
        signal.time = params.experiment.theoreticalSigs(i,:); % this stimulus theoretical signal
        [projected,components,~,~,output] = tscaFunc(Z,signal,[noise1 noise2],[params.TSCA.gammas(1:2) params.TSCA.gammas(3).*ones(1,length(noise2))],params.TSCA.numProj,params.TSCA.reduceComp); % call TSCA with this signal's parameters
        [~,I] = max(corr(abs(projected(1:9,:)'),signal.time')); % get the index of component with highest correlation to original signal
        mapTSCA1{i} = postProcess(rshp(components(:,I)));      
    end
    mapTSCA1 = cat(3,mapTSCA1{:});
    if nargin>1
        mapTSCA2 = cell(params.experiment.N,1); % with GLM denoising
        for i=1:params.experiment.N
            signal.time = params.experiment.theoreticalSigs(i,:); % this stimulus theoretical signal
            [projected,components,~,~,output] = tscaFunc(ZZ,signal,[noise1],params.TSCA.gammas(1:2),params.TSCA.numProj,params.TSCA.reduceComp); % call TSCA with this signal's parameters
            [~,I] = max(corr(abs(projected(1:9,:)'),signal.time')); % get the index of component with highest correlation to original signal
            mapTSCA2{i} = postProcess(rshp(components(:,I)));      
        end
        mapTSCA2 = cat(3,mapTSCA2{:});
    else
        mapTSCA2 = [];
    end    
    mapTSCA = struct('noGLM',mapTSCA1,'withGLM',mapTSCA2);
end

