function [result] = clusterEvaluation(result)
global brn params ump
% perform cluster evaluation and add to result 
    optimalMaps = params.experiment.optimalMaps.orig;
    N = params.experiment.N;
    X = cell(N,7); % cell array of clusters (rows - maps, cols - methods + optimal maps (7 total) )
    fn = fieldnames(result);
    RR = mean(params.experiment.optimalMaps.orig,3)>prctile(reshape(mean(params.experiment.optimalMaps.orig,3),[],1),80);
    for j=1:length(fn)+1 % iterate the methods and the optimal maps and perform cluster analysis between the maps
        if j~=length(fn)+1 % if not in optimal maps
            for k=1:N % iterate the maps of this method
                thisMap = result.(fn{j}).maps(:,:,k).*double(RR);
                [row,col] = find(thisMap>prctile(reshape(thisMap,[],1),99));
%                 mapp(:,:,k) = thisMap>prctile(reshape(thisMap,[],1),99);
                X{k,j} = [row,col];
            end
            [R,dbs,dbsI,dbns,dbnsI,DBI] = ClusterSimilarity(X(:,j));
            result.(fn{j}).clusterEval.R = R;
            result.(fn{j}).clusterEval.dbs = dbs;
            result.(fn{j}).clusterEval.dbsI = dbsI;
            result.(fn{j}).clusterEval.dbns = dbns;
            result.(fn{j}).clusterEval.dbnsI = dbnsI;
            result.(fn{j}).clusterEval.DBI = DBI;
        else % it is optimal maps
            for k=1:N % iterate the optimal maps
                [row,col] = find(optimalMaps(:,:,k)>prctile(reshape(optimalMaps(:,:,k),[],1),99));
%                 mapp(:,:,k) = optimalMaps(:,:,k)>prctile(reshape(optimalMaps(:,:,k),[],1),99);
                X{k,j} = [row,col];
            end
            [R,dbs,dbsI,dbns,dbnsI,DBI] = ClusterSimilarity(X(:,j));
            params.experiment.optimalMaps.clusterEval.R = R;
            params.experiment.optimalMaps.clusterEval.dbs = dbs;
            params.experiment.optimalMaps.clusterEval.dbsI = dbsI;
            params.experiment.optimalMaps.clusterEval.dbns = dbns;
            params.experiment.optimalMaps.clusterEval.dbnsI = dbnsI;
            params.experiment.optimalMaps.clusterEval.DBI = DBI;
        end
    end
    for k=1:N % iterate the maps and perform cluster analysis between the methods
        [R,dbs,dbsI,dbns,dbnsI,DBI] = ClusterSimilarity(X(k,:));
        result.clusterEvalAll{k}.R = R;
        result.clusterEvalAll{k}.dbs = dbs;
        result.clusterEvalAll{k}.dbsI = dbsI;
        result.clusterEvalAll{k}.dbns = dbns;
        result.clusterEvalAll{k}.dbnsI = dbnsI;
        result.clusterEvalAll{k}.DBI = DBI;
    end 

end

