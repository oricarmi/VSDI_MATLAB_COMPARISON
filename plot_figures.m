%% show simulation results - Figure 3
path = 'C:\Users\Ori\Desktop\Ori\2nd degree\matlab codez\vsdi - matlab\simulation results';
files = dir(path);
for m=3:length(files)
    load(fullfile(files(m).folder,files(m).name)); 
    TSCA = thisSNR_Summary{1};
    Tmax = thisSNR_Summary{2};
    ORIG = thisSNR_Summary{3};
    Corr = thisSNR_Summary{4};
    GLM = thisSNR_Summary{5};
    NADAV = thisSNR_Summary{6};
    Title = {'MSE','PSNR','CNR','MSSIM','Pearson''s Corr','CP'};
    Title2 = {'TSCA','Tmax','AOF','Corr','GLM','MPT'};
    retMaps = thisSNR_Summary{7};
    clusterEvalAll = thisSNR_Summary{8};
    figure(m*1000);
    figure(m*10000);
    for i=1:6 % iterate the 6 performance measures
        figure(m*1000);
        subplot(2,3,i)
        boxplot([squeeze(mean(TSCA(i,:,:),2)) squeeze(mean(Tmax(i,:,:),2)) squeeze(mean(ORIG(i,:,:),2))...
            squeeze(mean(Corr(i,:,:),2)) squeeze(mean(GLM(i,:,:),2)) squeeze(mean(NADAV(i,:,:),2))]...
            ,char({'TSCA';'Tmax';'AOF';'Corr';'GLM';'MPT'}));
        title(Title{i});
        figure(m*10000);
        subplot(2,3,i)
        imagesc(reshape(retMaps{i},40,40,3)); title([Title2{i} ' Retinotopic Map']);
    end
end

%% Figure 4 main
path = 'C:\Users\orica\OneDrive\Desktop\2nd degree\matlab codez\matlab - vsdi\simulation results 2';
path = 'C:\Users\Ori\Desktop\Ori\2nd degree\matlab codez\vsdi - matlab\simulation results 2';
path = 'C:\Users\orica\OneDrive\Desktop\2nd degree\matlab codez\matlab - vsdi\simulation results 2';
path = 'C:\Users\Ori\Desktop\Ori\2nd degree\matlab codez\vsdi - matlab\simulation results 2';
files = dir(path);
SNRorder = [5,3,4,7,6];
clusterEvalDBI = zeros(100,5,7);
clusterEvalS = zeros(100,5,7);
for m=1:length(SNRorder)
    load(fullfile(files(SNRorder(m)).folder,files(SNRorder(m)).name)); 
    thisclusterEvalDBI = thisSNR_Summary{9};
    thisclusterEvalS = thisSNR_Summary{10};
    % fix silhouette indexes
    for j=1:7 % iterate method
        for k=1:100 % iterate repitition
            if length(thisclusterEvalS{j}(k).s)~=240
                thisclusterEvalS{j}(k).s = [thisclusterEvalS{j}(k).s;-1*ones(240-length(thisclusterEvalS{j}(k).s),1)];
            end
        end
        clusterEvalDBI(:,m,j) = [thisclusterEvalDBI{j}.DBI]';
        clusterEvalS(:,m,j) = [mean([thisclusterEvalS{j}.s])]';
    end
    % end fix
end
Title2 = {'TSCA','Tmax','AOF','Corr','GLM','MPT','TSCA and GLM'};
figure("name","All",'Position', [50 100 1450 620]);
for i=1:7
    subplot 121;
    plot(1:5,nanmean(clusterEvalS(:,:,i)),'-o');hold on;
    subplot 122;
    plot(1:5,nanmean(clusterEvalDBI(:,:,i)),'-o'); hold on;

end
subplot 121;
xlim([0,6]);ylabel('adjusted Silhouette Index');
xticks(1:5);
xticklabels({'SNR=inf','SNR=5[dB]','SNR=0[dB]','SNR=-5[dB]','SNR=-10[dB]'});
set(gca,'XTickLabelRotation',45);
legend(Title2,'FontSize',12,'location','southwest');box off;
subplot 122;
xlim([0,6]);ylabel('adjusted DBI');
xticks(1:5);
xticklabels({'SNR=inf','SNR=5[dB]','SNR=0[dB]','SNR=-5[dB]','SNR=-10[dB]'});
set(gca,'XTickLabelRotation',45);
legend(Title2,'FontSize',12,'location','northwest');box off;
%%  Figure 5 main 
global bsln fs sz ump rot fgn brn brn0 frq cmap lgn scl lgn00 fnm c_f vms plt_on pt pc vc xs prms cfn cfn0 signal2_smooth basis params
[cf1 cfn trsh0]=strt_up(fname, n);  
result = Summary2.result;
params = Summary2.params;
fn = fieldnames(result);
t = linspace(0,1,100);
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
%% plot for supplmentary materil (supp. figure 1)
figure;
yyaxis left
plot(t,avgRes/params.experiment.N);ylabel('Z score');xlabel('time [sec]'); box off;
yyaxis right
plot(t,params.experiment.responseSig)
ylabel('amp [au]');xlabel('time [sec]');xticks([0:0.1:1]);
legend('Experimental Average Response','Theoretical Response');
% end for supplmentary
%% Figure 6,7 and supp 3,4,5 (need to load appropriate .mat results file)
result = Summary2.result;
params = Summary2.params;
fn = fieldnames(result);
retmapsAll = figure("name","retmapsAll",'Position', [100 100 1400 600]);
Titles ={'TSCA and GLM','TSCA','Tmax','AOF','Corr','GLM','MPT','Expected Retinotopy'};
flag=0;
for i=1:length(fn)-1 % iterate the methods  
    [~,r] = retinotopicMapFromIndividualMaps(result.(fn{i}).maps,0,fn{i},95);
    figure(retmapsAll); subplot(2,4,i)
    imf2(r);xlabel('');set(gca,'xticklabel',[]);set(gca,'xtick',[]);ylabel('');set(gca,'yticklabel',[]);set(gca,'ytick',[]);
    title(Titles{i}); hold on;
    if ~flag
        rectangle('Position',[4.5,5,1,0.1],'FaceColor',[1 1 1],'EdgeColor',[1 1 1]);
        text(4.4,4.8,'1mm','color',[1 1 1],'fontsize',12);
        flag=1;
    end
end
subplot(2,4,8);
imshow(expectedResponse);
title(Titles{end});

%% figure 8 main
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

%% figure 9 main - silhouette between moving bars and loc 8 
AllS_movingBars = [];
AllS_loc8 = [];
AllS_loc4 = [];
AllDBI_movingBars = [];
AllDBI_loc8 = [];
AllDBI_loc4 = [];
AllR_movingBars = [];
AllR_loc8 = [];
AllR_loc4 = [];
path = 'C:\Users\Ori\Desktop\Ori\2nd degree\matlab codez\vsdi - matlab\comparision results 2';
files = dir(path);
for i=[6,7,10,11,17,18] % iterate files (except last one which is NIR)  
    load(fullfile(files(i).folder,files(i).name)); % load summary
    result = Summary2.result;
    fn = fieldnames(result);
    params = Summary2.params;
    RR = mean(params.experiment.optimalMaps.orig,3)>prctile(reshape(mean(Summary2.params.experiment.optimalMaps.orig,3),[],1),85);
    [~,~,maxind2] = retinotopicMapFromIndividualMaps(result.TSCAwGLM.maps,1,'',93);
    maxind2 = maxind2.*RR;
        for k=1:Summary2.params.experiment.N % iterate the maps of this method
            [row,col] = find(maxind2==k);
            XX{k} = [row,col];
        end
    [ss,h] = silhouette(cat(1,XX{:}),makeClustVector(cellfun(@(x) size(x,1),XX))');
    [R,~,~,~,~,DBI] = ClusterSimilarity(XX);
    R = unique(R(R~=1));
    if contains(files(i).name,'horz')
        AllS_movingBars = [AllS_movingBars;ss];
        AllDBI_movingBars = [AllDBI_movingBars;DBI];
        AllR_movingBars = [AllR_movingBars;R(:)];
    elseif contains(files(i).name,'loc 8')
        AllS_loc8 = [AllS_loc8;ss];
        AllDBI_loc8 = [AllDBI_loc8;DBI];
        AllR_loc8 = [AllR_loc8;R(:)];
    else % loc 4
        AllS_loc4 = [AllS_loc4;ss];
        AllDBI_loc4 = [AllDBI_loc4;DBI];
        AllR_loc4 = [AllR_loc4;R(:)];
    end
    clear XX maxind R ss
end
figure;boxplot([AllS_movingBars;AllS_loc8],makeClustVector([length(AllS_movingBars);length(AllS_loc8)]),'labels',{'9 moving bars','8 location grid'});
ylabel('Si');title('silhouette values - 9 moving bars and 8 location grid');
figure;boxplot([AllDBI_movingBars,;AllDBI_loc8],makeClustVector([length(AllDBI_movingBars);length(AllDBI_loc8)]),'labels',{'9 moving bars','8 location grid'});
title('Davies-Bouldin Index- 9 moving bars and 8 location grid');
ylabel('DB Index');
figure;boxplot([AllS_movingBars;AllS_loc8],makeClustVector([length(AllS_movingBars);length(AllS_loc8)]),'labels',{'9 moving bars','8 location grid'});
ylabel('Si');title('silhouette values - 9 moving bars and 8 location grid');
figure;boxplot([AllR_movingBars,;AllR_loc8],makeClustVector([length(AllR_movingBars);length(AllR_loc8)]),'labels',{'9 moving bars','8 location grid'});
title('Davies-Bouldin Index- 9 moving bars and 8 location grid');
ylabel('DB Index');

%% like figure 4 but for the statistical measures (supp figure 2)
path = 'C:\Users\orica\OneDrive\Desktop\2nd degree\matlab codez\matlab - vsdi\simulation results 2';
path = 'C:\Users\Ori\Desktop\Ori\2nd degree\matlab codez\vsdi - matlab\simulation results 2';
files = dir(path);
SNRorder = [3,4,7,6];
MSE_all = zeros(100,4,7);PSNR_all = zeros(100,4,7);CNR_all = zeros(100,4,7);
MSSIM_all = zeros(100,4,7); CC_all = zeros(100,4,7); CP_all = zeros(100,4,7);
for m=1:length(SNRorder)
    load(fullfile(files(SNRorder(m)).folder,files(SNRorder(m)).name)); 
    for j=1:7 % iterate method
        MSE_all(:,m,j) = squeeze(mean(thisSNR_Summary{j}(1,:,:)));
        PSNR_all(:,m,j) = squeeze(mean(thisSNR_Summary{j}(2,:,:)));
        CNR_all(:,m,j) = squeeze(mean(thisSNR_Summary{j}(3,:,:)));
        MSSIM_all(:,m,j) = squeeze(mean(thisSNR_Summary{j}(4,:,:)));
        CC_all(:,m,j) = squeeze(mean(thisSNR_Summary{j}(5,:,:)));
        CP_all(:,m,j) = squeeze(mean(thisSNR_Summary{j}(6,:,:)));
    end
end
Title2 = {'TSCA','Tmax','AOF','Corr','GLM','MPT','TSCA and GLM'};
figure("name","All",'Position', [50 100 1450 620]);
for i=1:7
    subplot 231;
    plot(1:4,nanmean(MSE_all(:,:,i)),'-o');hold on;
    subplot 232;
    plot(1:4,nanmean(PSNR_all(:,:,i)),'-o'); hold on;
    subplot 233;
    plot(1:4,nanmean(CNR_all(:,:,i)),'-o'); hold on;    
    subplot 234;
    plot(1:4,nanmean(MSSIM_all(:,:,i)),'-o'); hold on;    
    subplot 235;
    plot(1:4,nanmean(CC_all(:,:,i)),'-o'); hold on;
    subplot 236;
    plot(1:4,nanmean(CP_all(:,:,i)),'-o'); hold on;
end
subplot 231;
xlim([0,5]);ylabel('MSE');
xticks(1:4);
xticklabels({'SNR=5[dB]','SNR=0[dB]','SNR=-5[dB]','SNR=-10[dB]'});
set(gca,'XTickLabelRotation',45);
legend(Title2,'FontSize',8,'location','northwest');box off;

subplot 232;
xlim([0,5]);ylabel('PSNR');
xticks(1:4);
xticklabels({'SNR=5[dB]','SNR=0[dB]','SNR=-5[dB]','SNR=-10[dB]'});
set(gca,'XTickLabelRotation',45);
% legend(Title2,'FontSize',12,'location','northwest');box off;

subplot 233;
xlim([0,5]);ylabel('CNR');
xticks(1:4);
xticklabels({'SNR=5[dB]','SNR=0[dB]','SNR=-5[dB]','SNR=-10[dB]'});
set(gca,'XTickLabelRotation',45);
% legend(Title2,'FontSize',12,'location','northwest');box off;

subplot 234;
xlim([0,5]);ylabel('MSSIM');
xticks(1:4);
xticklabels({'SNR=5[dB]','SNR=0[dB]','SNR=-5[dB]','SNR=-10[dB]'});
set(gca,'XTickLabelRotation',45);
% legend(Title2,'FontSize',12,'location','northwest');box off;

subplot 235;
xlim([0,5]);ylabel('CC');
xticks(1:4);
xticklabels({'SNR=5[dB]','SNR=0[dB]','SNR=-5[dB]','SNR=-10[dB]'});
set(gca,'XTickLabelRotation',45);
% legend(Title2,'FontSize',12,'location','northwest');box off;

subplot 236;
xlim([0,5]);ylabel('CP');
xticks(1:5);
xticklabels({'SNR=5[dB]','SNR=0[dB]','SNR=-5[dB]','SNR=-10[dB]'});
set(gca,'XTickLabelRotation',45);
% legend(Title2,'FontSize',12,'location','northwest');box off;


