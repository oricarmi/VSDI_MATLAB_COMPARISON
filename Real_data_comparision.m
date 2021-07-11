%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GENERAL %%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;
% fname = "G:\2020.01.21\m200121.mat"; n=4;
% fname = "E:\2018.12.18\MAT\m181218.mat"; n=2;
% fname = "E:\191119\m191119.mat"; n=2;
% fname = "E:\200121\m200121.mat"; n=4;
% fname = "E:\180904\m180904.mat"; n=4;
% fname ="G:\comparision results 2\loc 8,181218 n=2.mat"; n=2;
fname = "H:\2021.03.30\m210330.mat";n=4;
% fname = "G:\181218\m181218.mat"; n=2;
% fname = "E:\180801\m180801.mat"; n=3;
% fname = "D:\2019.07.10\m190710.mat"; n=2;
% fname = "C:\Users\orica\OneDrive\Desktop\2nd degree\matlab codez\matlab - vsdi\2021.03.03 - taget flankers\m210303.mat"; n=3;
% fname = "C:\Users\orica\OneDrive\Desktop\2nd degree\matlab codez\matlab - vsdi\2021.03.15 - target flankers\m210315.mat";n=3;
% fname = "G:\2019.12.11\m191211.mat"; n=3;
% fname = "H:\2020.11.30\m201130.mat"; n=2; 
% cd('C:\Users\Ori\Desktop\Ori\2nd degree\mtdt');
addpath("C:\Users\orica\Dropbox\fcns_and_decript");
addpath('C:\Users\orica\Dropbox\master degree\codes');
addpath("C:\Users\Ori\Dropbox\fcns_and_decript")
addpath("C:\Users\Ori\Dropbox\master degree\codes")
global bsln fs sz ump rot fgn brn brn0 frq cmap lgn scl lgn00 fnm c_f vms plt_on pt pc vc xs prms cfn cfn0 signal2_smooth basis params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
result = struct('TSCA',struct,'Tmax',struct,'AOF',struct,'Corr',struct,'GLM',struct,'Nadav',struct);
[result.TSCA.maps,result.Tmax.maps,result.AOF.maps,result.Corr.maps,result.GLM.maps,result.Nadav.maps,tscaBoth] = GenerateMaps(fname,n,9); % 3rd parameter: what, 8=loc8,9=loc9,92=mvngbars2hz
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
%% DB index for cluster similarity
X = cell(8,1); Xraw = cell(8,1); mapp = zeros(size(brn,1),size(brn,2),8);
for i=1:8
    [row,col] = find(result.TSCA.maps(:,:,i)>prctile(reshape(result.TSCA.maps(:,:,i),[],1),99));
    mapp(:,:,i) = result.TSCA.maps(:,:,i)>prctile(reshape(result.TSCA.maps(:,:,i),[],1),99);
    X{i} = [row,col];
    Xraw{i} = result.TSCA.maps(row,col,i);
end
[R,dbs,dbsI,dbns,dbnsI,DBI] = ClusterSimilarity(X);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%% Show Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath("C:\Users\orica\Dropbox\fcns_and_decript");
addpath('C:\Users\orica\Dropbox\master degree\codes');
addpath("C:\Users\Ori\Dropbox\fcns_and_decript")
addpath("C:\Users\Ori\Dropbox\master degree\codes")
global bsln fs sz ump rot fgn brn brn0 frq cmap lgn scl lgn00 fnm c_f vms plt_on pt pc vc xs prms cfn cfn0 signal2_smooth basis params
% description = Summary.description
fname = "G:\180905\m180905.mat";n=2; % take from description
% fname = "G:\181218\m181218.mat"; n=2;
% fname = "G:\191119\m191119.mat"; n=2;
% fname = "E:\200121\m200121.mat"; n=5;
[cf1 cfn trsh0]=strt_up(fname, n);  
result = Summary2.result;
params = Summary2.params;
fn = fieldnames(result);
expectedResponse = imread("C:\Users\orica\OneDrive\Desktop\2nd degree\comparison article\ExpectedResponsesImages\loc8_expected.bmp");
expectedResponse = imread("C:\Users\orica\OneDrive\Desktop\2nd degree\comparison article\ExpectedResponsesImages\loc9_expected3_leftright.bmp");
expectedResponse = imread("C:\Users\orica\OneDrive\Desktop\2nd degree\comparison article\ExpectedResponsesImages\loc4_expected.bmp");
%% Plot retmaps of a loaded Summary2
result = Summary2.result;
params = Summary2.params;
fn = fieldnames(result);
retmapsAll = figure("name","retmapsAll",'Position', [100 100 1400 600]);
% indmapsAll = figure("name","indmapsAll");
% Titles = {'(1)','(2)','(3)','(4)','(5)','(6)','(7)','(8)'};
Titles ={'TSCA and GLM','TSCA','Tmax','AOF','Corr','GLM','MPT','Expected Retinotopy'};
flag=0;
for i=1:length(fn)-1 % iterate the methods  
%     if i==1
%         for iii=1:size(result.(fn{i}).maps,3)
%             result.(fn{i}).maps(:,:,iii) = result.(fn{i}).maps(:,:,iii).*R;
%         end   
%     end
    [~,r] = retinotopicMapFromIndividualMaps(result.(fn{i}).maps,0,fn{i},95);
    figure(retmapsAll); subplot(2,4,i)
    imf2(r);xlabel('');set(gca,'xticklabel',[]);set(gca,'xtick',[]);ylabel('');set(gca,'yticklabel',[]);set(gca,'ytick',[]);
%     imagesc(rshp(result.(fn{i}).retinotopicMap));
    title(Titles{i}); hold on;
    if ~flag
        rectangle('Position',[4.5,5,1,0.1],'FaceColor',[1 1 1],'EdgeColor',[1 1 1]);
        text(4.4,4.8,'1mm','color',[1 1 1],'fontsize',12);
        flag=1;
    end
%     annotation(retmapsAll,'doublearrow',[0.1 0.1],[.1 0]);
%     if contains(fn{i},'nadav','ignorecase',true)
%         title('M.P.T');
%     elseif contains(fn{i},'tscano','ignorecase',true)
%         title('TSCA');
%     elseif contains(fn{i},'tscaw','ignorecase',true)
%         title('TSCA & GLM');
%     else
%         title(fn{i});
%     end
%     r = plotMaps(result.(fn{i}).maps,fn{i},1);
end
subplot(2,4,8);
imshow(expectedResponse);
title(Titles{end});
% with and without GLM
% figure;
% subplot 121
% [~,r] = retinotopicMapFromIndividualMaps(result.TSCAnoGLM.maps,0,'TSCA without GLM',92);
% imf2(r); title('TSCA without GLM');
% subplot 122
% [~,r] = retinotopicMapFromIndividualMaps(result.TSCAwGLM.maps,0,'TSCA with GLM',92);
% imf2(r);title('TSCA after GLM denoising');
%% average response of 1 condition
t = linspace(0,1,100);
% figure;
% plot(t,params.experiment.responseSig)
% ylabel('amp [au]');xlabel('time [sec]');
avgRes = zeros(1,params.experiment.T1);
zz = params.experiment.ZZ;
zz3d = rshp(zz);
avgAll = MinMaxNorm(max(zz3d,[],3));
[row,col] = find(avgAll>prctile(avgAll(:),90));
ROI = zeros(size(brn));
for ind=1:length(row)
    ROI(row(ind),col(ind)) = 1;
end
ROI = medfilt2(ROI,[5,5]);
AvgAll2 = avgAll.*ROI;
brn3d = reshape(repmat(MinMaxNorm(brn),1,1,3),[],3);
% avgAll(avgAll(:,1)<prctile(avgAll(:,1),95),:) = brn3d(avgAll(:,1)<prctile(avgAll(:,1),95),:);
%% plot for main 
figure;
subplot 221;
imf2(rshp(avgAll)); colorbar;
xlabel('');set(gca,'xticklabel',[]);set(gca,'xtick',[]);ylabel('');set(gca,'yticklabel',[]);set(gca,'ytick',[]);
rectangle('Position',[4.5,5,1,0.1],'FaceColor',[1 1 1],'EdgeColor',[1 1 1]);
text(4.4,4.8,'1mm','color',[1 1 1],'fontsize',12);
subplot 222;
imf2(brn3d);hold on;imf2(rshp(AvgAll2),prctile(avgAll(:),90)); colorbar;
xlabel('');set(gca,'xticklabel',[]);set(gca,'xtick',[]);ylabel('');set(gca,'yticklabel',[]);set(gca,'ytick',[]);
zzROI2 = [];
for i =1:params.experiment.N
    thisMap = result.TSCAwGLM.maps(:,:,i);
    [row,col] = find(thisMap>0.7);
%     [row,col] = find(ROI);
    thisZZ = zz(sub2ind(size(brn),row,col),(i-1)*params.experiment.T1+1:i*params.experiment.T1);
    zzROI2 = [zzROI2,mean(thisZZ)];
    avgRes = avgRes + mean(thisZZ);
end
avgRes = avgRes/(params.experiment.N); avgRes = avgRes-avgRes(1);
zzROI2 = 5*(zzROI2-min(zzROI2))/(max(zzROI2)-min(zzROI2));zzROI2=zzROI2-zzROI2(1);
subplot 223;plot( 0:0.01:8-0.01, zzROI2,'linewidth',1.5 );ylabel('Z score');xlabel('time [sec]');box off;
subplot 224;plot(t,avgRes,'linewidth',1.5);ylabel('Z score');xlabel('time [sec]'); box off;
xticks([0:0.1:1]);
%% plot for supplmentary materil
figure;
yyaxis left
plot(t,avgRes/params.experiment.N);ylabel('Z score');xlabel('time [sec]'); box off;
yyaxis right
plot(t,params.experiment.responseSig)
ylabel('amp [au]');xlabel('time [sec]');xticks([0:0.1:1]);
legend('Experimental Average Response','Theoretical Response');
% end for supplmentary
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
