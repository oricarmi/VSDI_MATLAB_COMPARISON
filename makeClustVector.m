function [clustVec] = makeClustVector(x)
% input: receive vector of repetitions of each digit, # of digits
% (clusters) is the length of x
clustNums = 1:length(x);
clustVec = repelem(clustNums,x);
end

