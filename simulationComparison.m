for i=1:5
    thisMap = signals(i).space;
    [row,col] = find(thisMap>prctile(reshape(thisMap,[],1),97));
    mapp(:,:,i) = thisMap>prctile(reshape(thisMap,[],1),97);
    X{i} = [row,col];
end
[R,dbs,dbsI,dbns,dbnsI,DBI] = ClusterSimilarity(X);
for i=1:5
    thisMap = mapTSCA(:,:,i);
    [row,col] = find(thisMap>prctile(reshape(thisMap,[],1),90));
    mapp(:,:,i) = thisMap>prctile(reshape(thisMap,[],1),90);
    X{i} = [row,col];
end
[Rtsca,dbstsca,dbsItsca,dbnstsca,dbnsItsca,DBItsca] = ClusterSimilarity(X);
for i=1:5
    thisMap = mapORIG(:,:,i);
    [row,col] = find(thisMap>prctile(reshape(thisMap,[],1),97));
    mapp(:,:,i) = thisMap>prctile(reshape(thisMap,[],1),97);
    X{i} = [row,col];
end
[RAOF,dbsAOF,dbsIAOF,dbnsAOF,dbnsIAOF,DBIAOF] = ClusterSimilarity(X);