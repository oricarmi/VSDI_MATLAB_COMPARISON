%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% SIMULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath("C:\Users\orica\Dropbox\fcns_and_decript");
addpath('C:\Users\orica\Dropbox\master degree\codes');
global fs 
T = 1000; fs = 100;
P = 1600;
m = sqrt(P);
signals = struct;
noises = struct;
colors = [0 0 1; 0 1 0; 0 1 1; 1 0 0; 1 0 1];
stim_time = 10:20:90;
noiseSig = [0 0.1 0.8 2 10];
t = linspace(0,(T-1)/fs,T);
flag=0; % flag for calculating signal cluster evaluation
for kk = 1:length(noiseSig) % iterate different noise sig
    runSummary = cell(7,100); % 7 methods, 100 repetitions
    retinotopicMapTSCA = zeros(m*m,3); % preallocate 40x40x3(hsv)
    retinotopicMapTSCA2 = zeros(m*m,3); % preallocate 40x40x3(hsv)
    retinotopicMapTmax = zeros(m*m,3); % preallocate 40x40x3(hsv)
    retinotopicMapAOF = zeros(m*m,3); % preallocate 40x40x3(hsv)
    retinotopicMapCorr = zeros(m*m,3); % preallocate 40x40x3(hsv)
    retinotopicMapGLM = zeros(m*m,3); % preallocate 40x40x3(hsv)
    retinotopicMapNADAV = zeros(m*m,3); % preallocate 40x40x3(hsv)
    for k=1:size(runSummary,2)
        %% construct signals
        [I,J] = ndgrid(1:m,1:m);
        % locs = [0.25 0.25; 0.25 0.75; 0.5 0.5; 0.75 0.25; 0.75 0.75];
        locs = [0.4 0.4; 0.4 0.6; 0.5 0.5; 0.6 0.4; 0.6 0.6]; ind2plot = [1,5,9,13,17];
        DC = 0.02; r = 4;
        amp = (1.5-0.5).*rand(1,5)+0.5; % amplitude distributed U~[0.5 1.5]
        for i=1:5
            signals(i).time = amp(i).*(2.5.*normpdf(0:0.1:(T-1)/10,stim_time(i),1)); % normpdf with peak amplitude 1, times random amplitude
            signals(i).space = MinMaxNorm(imgaussfilt(double((I-m*locs(i,1)).^2+(J-m*locs(i,2)).^2<r^2),sqrt(r)))+DC; % what is larger than r^2 is 1 (white), what is smaller is 0 (black)'
%             signals(i).space = double((I-m*locs(i,1)).^2+(J-m*locs(i,2)).^2<r^2); % what is larger than r^2 is 1 (white), what is smaller is 0 (black)
        end
        % plot
        figure(1122);ind2plot = [1,3,5,7,9];%ind2plot = [1,5,9,13,17];
        figure(2211);
        mapSigs = cat(3,signals(:).space);
        colorz = squeeze(hsv2rgb([1:5]./5,ones(1,5),ones(1,5)));
        [retinotopicMap,~,maxind] = retinotopicMapFromIndividualMaps(mapSigs,0,'sigs',85);
        figure;
        imagesc(rshp(retinotopicMap));xlabel('pixels'); ylabel('pixels');
        for i=1:5
            figure(1122);
            [retinotopicMap,~,maxind] = retinotopicMapFromIndividualMaps(mapSigs(:,:,i),0,['simulation_' num2str(i)],85);
            subplot(3,3,ind2plot(i));
            imagesc(rshp(retinotopicMap));xlabel('pixels'); ylabel('pixels');
            title(['sig. ' num2str(i) ' - spatial']);
            figure(2211);
            plot(t,signals(i).time,'color',colorz(i,:));hold on;
%             title(['sig. ' num2str(i) ' - temporal']); 
            xlabel('time [sec]'); ylabel('Amp. [au]');box off;
            xticks(1:2:9);
        end
        if ~flag % do this only once
            mapSigs = cat(3,signals(:).space);
            [retinotopicMap,~,maxind] = retinotopicMapFromIndividualMaps(mapSigs,0,'sigs',85);
            [clusterEvalSig,clusterEval2Sig] = ClusterEvaluationSim2(maxind);
            flag=1;
        end

        %% construct noise
        [I1,I2] = ndgrid([repmat(linspace(0,2*pi,m/2),1,2)]',[repmat(linspace(0,2*pi,m/2),1,2)]);
        noises(1).time = normrnd(0,noiseSig(kk),T,1);
        noises(1).space = MinMaxNorm(normrnd(0,1,m,m))+DC;%cos(I1);
        freqs = 2*pi.*[2.8:0.05:3.2]'; % bandwidth around 3 hz
        phase = 2*pi.*rand(length(freqs),1); % random phase
        amp = noiseSig(kk)*4.*rand(size(phase)); % random amplitude
        noises(2).time = mean(repmat(amp,1,T).*cos(freqs.*t+phase))+normrnd(0,noiseSig(kk)*0.05,1,T); % create almost periodic signal around 3hz
        noises(2).space = MinMaxNorm(cos(3.*I1+2*pi*rand(1)))+DC;
        freqs = 2*pi.*[0.47:0.05:0.87]'; % bandwidth around 0.67 hz
        phase = 2*pi.*rand(length(freqs),1);% random phase
        amp = noiseSig(kk)*4.*rand(size(phase));% random amplitude
        noises(3).time = mean(repmat(amp,1,T).*cos(freqs.*t+phase))+normrnd(0,noiseSig(kk)*0.05,1,T); % create almost periodic signal around 0.67 hz %noiseSig(kk)*cos(2*pi*0.67.*t'+2*pi*rand(1))+normrnd(0,noiseSig(kk)/2,T,1);
        noises(3).space = MinMaxNorm(normrnd(0,1,m,m))+DC;
        figure; ind2plot = [1,3,5];
        for i=1:length(noises)
            subplot(3,2,ind2plot(i))
            imagesc(noises(i).space); colormap(gray);
            xlabel('pixels'); ylabel('pixels');
%             title(['noise ' num2str(i) ' - spatial']);
            subplot(3,2,ind2plot(i)+1)
            plot(t,noises(i).time,'k');
%             title(['noise ' num2str(i) ' - temporal']); 
            xlabel('time [sec]'); ylabel('Amp. [au]');
        end
        %% construct Z and show 9 frames
        Z = zeros(m*m,T); % preallocate memory
        for i = 1:T
            Z(:,i) = reshape( ...
                signals(1).time(i)*signals(1).space+signals(2).time(i)*signals(2).space+...
                signals(3).time(i)*signals(3).space+signals(4).time(i)*signals(4).space+...
                signals(5).time(i)*signals(5).space+...
                noises(1).time(i)*noises(1).space+...
                noises(2).time(i)*noises(2).space+...
                noises(3).time(i)*noises(3).space...
                ,[],1);
        end
        % for i=1:size(Z,1) % iterate pixels and add exponential trend)
        %     ZZ(i,:) = Z(i,:)+(0.5.*exp(-t./4));
        % end
%         figure; %suptitle('the recorded signal (Z) at random frames');
%         ind2show = sort([100;501;899;randi(1000,7,1)]);
%         for i=1:length(ind2show)
%             subplot(2,5,i)
%             imshow(reshape(Z(:,ind2show(i)),m,[])); colormap(gray);
%             title(['Frame ' num2str(ind2show(i))]);
%         end
        theoreticalSig = zeros(length(signals(1).time),length(signals));
        for i=1:length(signals)
            theoreticalSig(:,i) = normpdf(0:0.1:(T-1)/10,stim_time(i),1);
        end
        [ZZ,ZZZ,ZZZZ,betas] = GLM_VSDI(Z,[0.67 3],theoreticalSig);
%         ind2show = sort([100;501;899;randi(1000,7,1)]);
%         for i=1:length(ind2show)
%             subplot(2,5,i)
%             imshow(reshape(ZZZZ(:,ind2show(i)),m,[])); colormap(gray);
%             title(['Frame ' num2str(ind2show(i))]);
%         end
        figure;
        subplot(2,2,1);
        imagesc(signals(2).space);colormap(gray);xlabel('pixels'); ylabel('pixels');
        subplot(2,2,2);
        plot(t(1:200),signals(1).time(1:200)); xlabel('time [sec]'); ylabel('amp [au]');
        subplot(2,2,3);
        imagesc(rshp(Z(:,300)));colormap(gray);xlabel('pixels'); ylabel('pixels');
        subplot(2,2,4);
        plot(t(1:200),(Z(sub2ind([40,40],repmat([15:17]',3,1),repmat([23:25]',3,1)),201:400))); xlabel('time [sec]'); ylabel('amp [au]');
        %% ============================= TSCA no GLM ===================================
        noiseNew.time = eye(T)/T; noiseNew.space = []; mapTSCA = zeros(m,m,length(signals));
        noise2New.time = createToeplitz(3,0.1,1,1,T); noise2New.space = [];
        noise3New.time = createToeplitz(0.67,0.1,1,1,T); noise3New.space = [];
        for i=1:length(signals)
            sig.time = theoreticalSig(:,i)';
            [projected,components,D,Alpha,output] = tscaFunc(Z,sig,[noiseNew noise2New noise3New],[1 -0.2*ones(1,3)],100,1);
            %     tscaAnalyze(output,3,[],0,T);
            [~,I] = max(corr(abs(projected(1:4,:)'),theoreticalSig(:,i))); % get the index of component with highest correlation to original signal
            mapTSCA(:,:,i) = MinMaxNorm(abs(reshape(components(:,I),m,m)));
            [tscaMSE(i),tscaPSNR(i),tscaCNR(i),tscaMSSIM(i),tscaCorr(i),tscaCP(i)] = getPerformance(mapTSCA(:,:,i),signals(i).space,signals(i).space>10*DC,signals(i).space<=10*DC);
        end
        [r,~,maxind] = retinotopicMapFromIndividualMaps(mapTSCA,0,'TSCA',85);
        retinotopicMapTSCA = retinotopicMapTSCA + r;
        [clusterEvalTSCA(k),clusterEval2TSCA(k)] = ClusterEvaluationSim2(maxind);

        %% ============================= TSCA with GLM ===================================
        for i=1:length(signals)
            sig.time = theoreticalSig(:,i)';
            [projected,components,D,Alpha,output] = tscaFunc(ZZ,sig,[noiseNew],[1 -0.2*ones(1,1)],100,1);
            %     tscaAnalyze(output,3,[],0,T);
            [~,I] = max(corr(abs(projected(1:4,:)'),theoreticalSig(:,i))); % get the index of component with highest correlation to original signal
            mapTSCA2(:,:,i) = MinMaxNorm(abs(reshape(components(:,I),m,m)));
            [tsca2MSE(i),tsca2PSNR(i),tsca2CNR(i),tsca2MSSIM(i),tsca2Corr(i),tsca2CP(i)] = getPerformance(mapTSCA2(:,:,i),signals(i).space,signals(i).space>10*DC,signals(i).space<=10*DC);
        end
        [r,~,maxind] = retinotopicMapFromIndividualMaps(mapTSCA2,0,'TSCA+GLM',85);
        retinotopicMapTSCA2 = retinotopicMapTSCA2 + r;
        [clusterEvalTSCA2(k),clusterEval2TSCA2(k)] = ClusterEvaluationSim2(maxind);        
        % [maxmap,maxind] = max(map,[],3); maxmap = maxmap(:); maxind = maxind(:);
        % map22 = colors(maxind,:); map22(maxmap(:)< prctile(maxmap(:),70),:) = repmat([0 0 0],length(find(maxmap(:)< prctile(maxmap(:),70))),1);
        % retinotopicMap = hsv2rgb(maxind/length(unique(maxind)),ones(size(maxind)),(maxmap-min(maxmap))./(max(maxmap)-min(maxmap)));
        % figure;imagesc(reshape(retinotopicMap,40,40,3));
        
        %% calc performance measures between theoretical signal and Z
        for i=1:length(signals)
            [~,I] = max(signals(i).time);
            mapAOF(:,:,i) = MinMaxNorm(reshape(mean(Z(:,I-25:I+25),2),40,40));
            [origMSE(i),origPSNR(i),origCNR(i),origMSSIM(i),origCorr(i),origCP(i)] = getPerformance( mapAOF(:,:,i),signals(i).space,signals(i).space>10*DC,signals(i).space<=10*DC);
            %     figure;imagesc(reshape(mean(Z(:,I-25:I+25),2),40,40));
        end
        [r,~,maxind] = retinotopicMapFromIndividualMaps(mapAOF,0,'AOF',85);
        retinotopicMapAOF = retinotopicMapAOF + r;
        [clusterEvalAOF(k),clusterEval2AOF(k)] = ClusterEvaluationSim2(maxind);        
        %% ========================== T_max method =============================
        refff = normpdf(0:0.1:(T-1)/10,10,1);
        r = zeros(1,size(Z,1)); tmax = zeros(size(r));
        for i=1:size(Z,1)
            [rtemp,lags] = xcorr(Z(i,:),refff);
            [r(i),I] = max(rtemp);
            tmax(i) = lags(I);
        end
        tmax(r<mean(r)+0.2*std(r)) = -100;
        A = tmax;
        B = [-100, 0:200:800];
        [~,I] = min(abs(bsxfun(@minus,A,B')));
        Anew = B(I);
        for i=1:length(signals)
            temp = zeros(size(Anew))';
            temp(Anew == B(i+1)) = 1;
            temp(temp==1) = r(temp==1);
            if max(temp(:))~=min(temp(:)) % perform maxminnorm if possible
                mapTmax(:,:,i) = MinMaxNorm(reshape(temp,40,40));
            else
                mapTmax(:,:,i) = reshape(temp,40,40);
            end
            [T_MSE(i),T_PSNR(i),T_CNR(i),T_MSSIM(i),T_corr(i),T_CP(i)] = getPerformance(mapTmax(:,:,i),signals(i).space,signals(i).space>10*DC,signals(i).space<=10*DC);
        end
        [r,~,maxind] = retinotopicMapFromIndividualMaps(mapTmax,0,'Tmax',85);
        retinotopicMapTmax = retinotopicMapTmax + r;
        [clusterEvalTmax(k),clusterEval2Tmax(k)] = ClusterEvaluationSim2(maxind);  
        % tmax2 = Anew';
        % tmax2(tmax2>=0) = (tmax2(tmax2>=0)./100+2)./2; tmax22 = zeros(length(tmax2),3);% transform 0,200,400,600,800 to 1,2,3,4,5
        % tmax22(tmax2>0,:) = colors(tmax2(tmax2>0),:);
        % figure;imagesc(reshape(tmax22,40,40,3));
        
        %% ======================== Correlation method =========================
        normalizedZ=Z;
        mapCorr = normalizedZ*theoreticalSig; mapCorr = rshp(mapCorr);
        for i=1:length(signals)
            mapCorr(:,:,i) = MinMaxNorm(mapCorr(:,:,i));
            [Corr_MSE(i),Corr_PSNR(i),Corr_CNR(i),Corr_MSSIM(i),Corr_corr(i),Corr_CP(i)] = getPerformance(mapCorr(:,:,i),signals(i).space,signals(i).space>10*DC,signals(i).space<=10*DC);
        end
        [r,~,maxind] = retinotopicMapFromIndividualMaps(mapCorr,0,'Corr',85);
        retinotopicMapCorr = retinotopicMapCorr + r;
        [clusterEvalCorr(k),clusterEval2Corr(k)] = ClusterEvaluationSim2(maxind);
        
        %%  =========================== GLM method =============================
        mapGLM = rshp(betas(6:end,:)');
        for i=1:length(signals)
            mapGLM(:,:,i) = MinMaxNorm(mapGLM(:,:,i));
            [GLM_MSE(i),GLM_PSNR(i),GLM_CNR(i),GLM_MSSIM(i),GLM_corr(i),GLM_CP(i)] = getPerformance(mapGLM(:,:,i),signals(i).space,signals(i).space>10*DC,signals(i).space<=10*DC);
        end
        [r,~,maxind] = retinotopicMapFromIndividualMaps(mapGLM,0,'GLM',85);
        retinotopicMapGLM = retinotopicMapGLM + r;
        [clusterEvalGLM(k),clusterEval2GLM(k)] = ClusterEvaluationSim2(maxind);
        %% ========================== Nadav's method ===========================
        for i=1:length(signals)
            thisSig = Z(:,stim_time(i)*10-99:stim_time(i)*10+100);
            [ind_mx1 w_mx i_rx w_rx]=xplmts(thisSig,[],mean(thisSig(:))+0.2*std(thisSig(:)),[0.9 1.1],[mean(thisSig(:))+1*std(thisSig(:)) 30], 10);
            if max(ind_mx1(:))~=min(ind_mx1(:)) % if can perform min max norm
                mapMPT(:,:,i) = MinMaxNorm(reshape(ind_mx1,40,40));
            else
                mapMPT(:,:,i) = reshape(ind_mx1,40,40);
            end
            %     figure;imagesc(reshape(mapN(:,:,i),40,40));colormap(gray);
            [nadavMSE(i),nadavPSNR(i),nadavCNR(i),nadavMSSIM(i),nadavCorr(i),nadavCP(i)] = getPerformance(mapMPT(:,:,i),signals(i).space,signals(i).space>10*DC,signals(i).space<=10*DC);
        end
        [r,~,maxind] = retinotopicMapFromIndividualMaps(mapMPT,0,'MPT',85);
        retinotopicMapNADAV = retinotopicMapNADAV + r;
        [clusterEvalNadav(k),clusterEval2Nadav(k)] = ClusterEvaluationSim2(maxind);
        %% finish up this iteration
        % close all;
        runSummary{1,k} = [origMSE;origPSNR;origCNR;origMSSIM;origCorr;origCP];
        runSummary{2,k} = [nadavMSE;nadavPSNR;nadavCNR;nadavMSSIM;nadavCorr;nadavCP];
        runSummary{3,k} = [tscaMSE;tscaPSNR;tscaCNR;tscaMSSIM;tscaCorr;tscaCP];
        runSummary{4,k} = [T_MSE;T_PSNR;T_CNR;T_MSSIM;T_corr;T_CP];
        runSummary{5,k} = [Corr_MSE;Corr_PSNR;Corr_CNR;Corr_MSSIM;Corr_corr;Corr_CP];
        runSummary{6,k} = [GLM_MSE;GLM_PSNR;GLM_CNR;GLM_MSSIM;GLM_corr;GLM_CP];  
        runSummary{7,k} = [tsca2MSE;tsca2PSNR;tsca2CNR;tsca2MSSIM;tsca2Corr;tsca2CP];
    end
    %% Compare results (all runs)
    ORIG = cat(3,runSummary{1,:});
    NADAV = cat(3,runSummary{2,:}); NADAV(isnan(NADAV)) = 0;
    TSCA = cat(3,runSummary{3,:});
    Tmax = cat(3,runSummary{4,:}); Tmax(isinf(Tmax)) = 100; Tmax(isnan(Tmax)) = 0;
    Corr = cat(3,runSummary{5,:});
    GLM = cat(3,runSummary{6,:});
    TSCAwGLM = cat(3,runSummary{7,:});
    clusterEvalAll = {clusterEvalTSCA,clusterEvalAOF,clusterEvalTmax,clusterEvalCorr,clusterEvalGLM,clusterEvalNadav,clusterEvalTSCA2};
    clusterEvalAll2 = {clusterEval2TSCA,clusterEval2AOF,clusterEval2Tmax,clusterEval2Corr,clusterEval2GLM,clusterEval2Nadav,clusterEval2TSCA2};
    Title = {'MSE','PSNR','CNR','MSSIM','Corr','CP'};
    Title2 = {'TSCA','T_{max}','AOF','Corr','GLM','MPT'};
    retMaps = {retinotopicMapTSCA./k,retinotopicMapTmax./k,retinotopicMapAOF./k,retinotopicMapCorr./k,retinotopicMapGLM./k,retinotopicMapNADAV./k,retinotopicMapTSCA2./k};
    figure(kk*10);
    figure(kk*100);
    for i=1:6 % iterate the 6 performance measures
        figure(kk*10);
        subplot(2,3,i)
        boxplot([squeeze(mean(TSCA(i,:,:),2)) squeeze(mean(Tmax(i,:,:),2)) squeeze(mean(ORIG(i,:,:),2))...
            squeeze(mean(Corr(i,:,:),2)) squeeze(mean(GLM(i,:,:),2)) squeeze(mean(NADAV(i,:,:),2))...
            squeeze(mean(TSCAwGLM(i,:,:),2))]...
            ,char({'TSCA';'Tmax';'AOF';'Corr';'GLM';'MPT';['GLM+','TSCA']}));
        title(Title{i});
        figure(kk*100);
        subplot(2,3,i)
        imagesc(reshape(retMaps{i},40,40,3)); title([Title2{i} ' Retinotopic Map']);
    end
    figure(kk*1000);
    boxplot([cat(1,clusterEvalAll{1}(:).DBI),cat(1,clusterEvalAll{2}(:).DBI),cat(1,clusterEvalAll{3}(:).DBI),cat(1,clusterEvalAll{4}(:).DBI),cat(1,clusterEvalAll{5}(:).DBI),cat(1,clusterEvalAll{6}(:).DBI),cat(1,clusterEvalAll{7}(:).DBI)],char({'TSCA';'Tmax';'AOF';'Corr';'GLM';'MPT';'TSCA+GLM'}));
    title('boxplot of Davies-Bouldin Index of all methods - simulation');
    ylabel('DB Index');
    figure(kk*10000);
    boxplot([cat(1,clusterEvalAll2{1}(:).s);cat(1,clusterEvalAll2{2}(:).s);cat(1,clusterEvalAll2{3}(:).s);cat(1,clusterEvalAll2{4}(:).s);cat(1,clusterEvalAll2{5}(:).s);cat(1,clusterEvalAll2{6}(:).s);cat(1,clusterEvalAll2{7}(:).s)]...
        ,[makeClustVector([length(cat(1,clusterEvalAll2{1}(:).s));length(cat(1,clusterEvalAll2{2}(:).s));length(cat(1,clusterEvalAll2{3}(:).s));length(cat(1,clusterEvalAll2{4}(:).s));length(cat(1,clusterEvalAll2{5}(:).s));length(cat(1,clusterEvalAll2{6}(:).s));length(cat(1,clusterEvalAll2{7}(:).s))])]','labels',char({'TSCA';'Tmax';'AOF';'Corr';'GLM';'MPT';'TSCA+GLM'}));
    title('boxplot of silhouette index of all methods - simulation');
    ylabel('Si');
    thisSNR_Summary = {TSCA,Tmax,ORIG,Corr,GLM,NADAV,TSCAwGLM,retMaps,clusterEvalAll,clusterEvalAll2};
    save(fullfile('C:\Users\Ori\Desktop\Ori\2nd degree\matlab codez\vsdi - matlab\simulation results 2',...
        ['SimulationSummary_NoiseSig=' num2str(noiseSig(kk)) '.mat']),'thisSNR_Summary');
end