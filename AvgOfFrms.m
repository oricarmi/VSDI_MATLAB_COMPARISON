function [mapAOF] = AvgOfFrms(ZZ)
global params brn
% Average of frames
mapAOF2 = zeros(size(brn,1),size(brn,2),4);
for i=1:params.experiment.N % iterate the conditions
    mapAOF(:,:,i) = postProcess(mean(ZZ(:,:,(i-1)*params.experiment.T1+params.AOF.numFramesFrom:(i-1)*params.experiment.T1+params.AOF.numFramesUntil),3));
end
