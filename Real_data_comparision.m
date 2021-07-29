%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GENERAL %%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;
fname = "E:\181218\m181218.mat"; n=2; % path to .mat file
addpath("C:\Users\orica\Dropbox\fcns_and_decript");
addpath('C:\Users\orica\Dropbox\master degree\codes');
addpath("C:\Users\Ori\Dropbox\fcns_and_decript")
addpath("C:\Users\Ori\Dropbox\master degree\codes")
global bsln fs sz ump rot fgn brn brn0 frq cmap lgn scl lgn00 fnm c_f vms plt_on pt pc vc xs prms cfn cfn0 signal2_smooth basis params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
result = struct('TSCA',struct,'Tmax',struct,'AOF',struct,'Corr',struct,'GLM',struct,'Nadav',struct);
[result.TSCA.maps,result.Tmax.maps,result.AOF.maps,result.Corr.maps,result.GLM.maps,result.Nadav.maps,tscaBoth] = GenerateMaps(fname,n,8); % 3rd parameter: what, 8=loc8,9=loc9,92=mvngbars2hz
% <----- generate retinotopic maps from the individual maps
fn = fieldnames(result);
for i=1:length(fn) % iterate the methods 
    [~,result.(fn{i}).retinotopicMap] = retinotopicMapFromIndividualMaps(result.(fn{i}).maps,1,fn{i},90);
end
[~,resultBoth] = retinotopicMapFromIndividualMaps(tscaBoth.withGLM,1,'GLM+TSCA',90);
% ---- RUN UNTIL HERE -----
% ---->
[result.TSCA.performance,result.Tmax.performance,result.AOF.performance,result.Corr.performance,result.GLM.performance,result.Nadav.performance] = performanceRealData(result);
result = clusterEvaluation(result);
Summary = struct('params',params,'result',result','description','loc 4,191119 n=2'); 
%% silhouette + DB index, cluster analysis with maxind
% path = 'E:\comparision results 2';
path = 'D:\dataForComparison\comparision results 2';
files = dir(path);
AllRMats = cell(1,7);
AllDBI = cell(1,7);
AllS = cell(1,7);
global params brn
for i=3:length(files)-1 % iterate files (except last one which is NIR)  
    load(fullfile(files(i).folder,files(i).name)); % load summary
    result = Summary2.result;
    fn = fieldnames(result);
    params = Summary2.params;
    RR = mean(params.experiment.optimalMaps.orig,3)>prctile(reshape(mean(Summary2.params.experiment.optimalMaps.orig,3),[],1),85);
    for j=1:length(fn)-1 % iterate the methods and the optimal maps and perform cluster analysis between the maps
        [~,~,maxind{j}] = retinotopicMapFromIndividualMaps(result.(fn{j}).maps,0,'',93);
        maxind{j} = maxind{j}.*RR;
        for k=1:Summary2.params.experiment.N % iterate the maps of this method
            [row,col] = find(maxind{j}==k);
            X{k,j} = [row,col];
        end
    end
    for j=1:length(fn)-1
        figure;
        [s{j},h] = silhouette(cat(1,X{:,j}),makeClustVector(cellfun(@(x) size(x,1),X(:,j)))');
        AllS{j} = [AllS{j};s{j}];
        [R{j},~,~,~,~,DBI{j}] = ClusterSimilarity(X(:,j));
        if contains(files(i).name,'loc 8')
            AllRMats{j} = cat(3,AllRMats{j},R{j});
        end
        AllDBI{j} = [AllDBI{j};DBI{j}];
    end
    close all
    clear X maxind R DBI s
%     figure;
%     for i=1:7
%         subplot(3,3,i)
%         imagesc(maxind{i});
%     end
end
set(0,'DefaultTextInterpreter','tex')
figure;
boxplot([AllDBI{2},AllDBI{3},AllDBI{4},AllDBI{5},AllDBI{6},AllDBI{7},AllDBI{1}],char({'TSCA';'Tmax';'AOF';'Corr';'GLM';'MPT';'TSCA+GLM'}));
title('Davies-Bouldin Index');
ylabel('DB Index');
figure;
AllS = cellfun(@(x) max(x, -ones(size(x))),AllS,'uniformoutput',false);
boxplot([AllS{2};AllS{3};AllS{4};AllS{5};AllS{6};AllS{7};AllS{1}],[makeClustVector(cellfun(@(x) size(x,1),AllS))]','labels',char({'TSCA';'Tmax';'AOF';'Corr';'GLM';'MPT';'TSCA+GLM'}));
title('Silhouette index');
ylabel('Si');
meanR = cell(1,7); stdR = cell(1,7);
for i=1:7
    meanR{i} = mean(AllRMats{i},3);
    stdR{i} = std(AllRMats{i},0,3);
end
meanDBI = cellfun(@mean,AllDBI);
meanS = cellfun(@nanmean,AllS);
for i=1:7
    writematrix(meanR{i},'C:\Users\Ori\Desktop\Ori\2nd degree\matlab codez\vsdi - matlab\comparision results 2\meanR.xls','Sheet',i)
    writematrix(stdR{i},'C:\Users\Ori\Desktop\Ori\2nd degree\matlab codez\vsdi - matlab\comparision results 2\stdR.xls','Sheet',i)
end
%% figure 7 main
figure;
subplot 121
boxplot([AllmedS{1};AllmedS{2};AllmedS{3};AllmedS{4};AllmedS{5};AllmedS{6};AllmedS{7}],[makeClustVector(cellfun(@(x) size(x,1),AllmedS))]','labels',char({'T&G';'TSCA';'Tmax';'AOF';'Corr';'GLM';'MPT'}),'symbol', '');
% title('Median Silhouette index');
ylabel('median Adjusted SI');set(gca,'XTickLabelRotation',45);
box off
subplot 122
boxplot([AllDBI{1},AllDBI{2},AllDBI{3},AllDBI{4},AllDBI{5},AllDBI{6},AllDBI{7}],char({'T&G';'TSCA';'Tmax';'AOF';'Corr';'GLM';'MPT'}),'symbol', '');
ylabel('Adjusted DBI');set(gca,'XTickLabelRotation',45);
box off
