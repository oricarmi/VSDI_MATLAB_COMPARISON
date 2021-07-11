function [res] = ClusterEvaluationSim(map)
% Perform cluster similiary analysis from simulated data
if length(size(map))<3 % reshape so that 3rd dimension is map of different signal
    map = reshape(map,40,40,[]);
end
X = cell(1,size(map,3));
for i=1:size(map,3) % iterate maps 
    thisMap = map(:,:,i);
    [row,col] = find(thisMap>prctile(reshape(thisMap,[],1),97));
%     mapp(:,:,i) = thisMap>prctile(reshape(thisMap,[],1),97);
    X{i} = [row,col];
end
[R,db,dbsI,dbns,dbnsI,DBI] = ClusterSimilarity(X);
res = struct('R',R,'db',db,'dbsI',dbsI,'dbns',dbns,'dbnsI',dbnsI,'DBI',DBI);
end

