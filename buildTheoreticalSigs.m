function [theoreticalSig] = buildTheoreticalSigs()
% Build theoretical signals
    global params
    theoreticalSig = zeros(params.experiment.N,params.experiment.T);
    for i=1:params.experiment.N
        theoreticalSig(i,:) = [zeros(1,(i-1)*params.experiment.T1) params.experiment.responseSig zeros(1,(params.experiment.N-(i))*params.experiment.T1)];
    end
end

