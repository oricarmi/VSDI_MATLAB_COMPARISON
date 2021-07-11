function [maps] = Correlation2All(Z)
% Calculate the correlation to all. Z - numel(brn) x time
    global basis brn params
%     Z = (Z-mean(Z,2))./std(Z,[],2);
    maps = Z*params.experiment.theoreticalSigs'; maps = rshp(maps);
    for i=1:size(maps,3)
        maps(:,:,i) = postProcess(maps(:,:,i));
    end
%     plotMaps(maps,what);
end

