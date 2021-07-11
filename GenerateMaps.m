function [mapTSCA,mapTmax,mapAOF,mapCorr,mapGLM,mapNadav,mapTSCAboth] = GenerateMaps(fname,n,what)
% A function that generates maps from all methods
global bsln fs sz ump rot fgn brn brn0 frq cmap lgn scl lgn00 fnm c_f vms plt_on pt pc vc xs prms cfn cfn0 basis params
[cf1 cfn trsh0]=strt_up(fname, n);  
% allBrn = input('Whole Brain? [0,1]');
% if allBrn == 0
%    cfn=pk_zn2(cfn,2);
% end
 % <---------------- generate signals
%  try 
%      if what/10<1 %8locs/9locs/... (not 2hz)
%          load('responseSig.mat');
%      else
%          load('responseSig2Hz.mat');
%      end
%  catch
%     if isempty(basis)
%         if what/10<1 %8locs/9locs/... (not 2hz)
%             t = 0:1/fs:(1-1/fs); % 1 hz
%         else %2hz
%             t = 0:1/fs:(0.5-1/fs); % 2 hz
%         end
%         alphas = [50 200; 10 100; 10 50;50 150]./1000;
%         randomAlphas = [...
%             alphas(1,1)+rand(15000,1)*(alphas(1,2)-alphas(1,1))...
%             alphas(2,1)+rand(15000,1)*(alphas(2,2)-alphas(2,1))...
%             alphas(3,1)+rand(15000,1)*(alphas(3,2)-alphas(3,1))...
%             alphas(4,1)+rand(15000,1)*(alphas(4,2)-alphas(4,1))];
%         r = zeros(size(randomAlphas,1),length(t));
%         for i=1:size(r,1) % iterate random signals
%             for j=1:size(r,2) % iterate time
%                 if t(j)<=randomAlphas(i,1) || t(j)>=(randomAlphas(i,1)+randomAlphas(i,2)+randomAlphas(i,3)+randomAlphas(i,4)) % if before latency or after decay, it is zero
%                     continue
%                 elseif randomAlphas(i,1)<=t(j) && t(j)<=(randomAlphas(i,1)+randomAlphas(i,2)) % if in rise time
%                     r(i,j) = 0.5*(1-cos(pi*(t(j)-randomAlphas(i,1))/randomAlphas(i,2)));
%                 elseif (randomAlphas(i,1)+randomAlphas(i,2))<=t(j) && t(j)<=(randomAlphas(i,1)+randomAlphas(i,2)+randomAlphas(i,3)) % if in plateau time
%                     r(i,j) = 1;
%                 else %if (randomAlphas(i,1)+randomAlphas(i,2)+randomAlphas(i,3))<=t(j) && t(j)<=(randomAlphas(i,1)+randomAlphas(i,2)+randomAlphas(i,3)+randomAlphas(i,4)) % if in decaying time
%                     r(i,j) = 0.5*(1+cos(pi*(t(j)-randomAlphas(i,1)-randomAlphas(i,2)-randomAlphas(i,3))/randomAlphas(i,4)));
%                 end
%             end
%         end
%         [~,~,V] = svd(r,'econ');
%         basis = V(:,1:20);
%     end
%     responseSig = [basis(:,1)-basis(:,2)]'; % figure;plot(responseSig);
%  end
% -------------------- > 
% < --- concatenate all 10 seconds to run in video
Z = [];
switch what
    case 92  % moving bars
        for i=4:length(cfn)
            Z = [Z cfn{i}];
        end
    case 9 % loc 9
        for i=2:length(cfn)
            Z = [Z cfn{i}];
        end
    otherwise % loc 8 or other
        for i=3:length(cfn)
            Z = [Z cfn{i}];
        end
end 
% if input('Choose noise freqs? [0/1]')
%     ChooseNoiseFreqs(Z);
% end
% <--- define parameters 
try
    defineParameters("C:\Users\orica\OneDrive\Desktop\2nd degree\matlab codez\matlab - vsdi\VSDI-MATLAB\paramsori.csv",what,rshp(Z)); % personal pc
catch
    defineParameters("C:\Users\Ori\Desktop\Ori\2nd degree\matlab codez\vsdi - matlab\VSDI-MATLAB\paramsori.csv",what,rshp(Z)); % lab pc
end
% ---->
ZZ = preProcess(Z);
implay(rshp(ZZ),20); % play video
[ZZZ,ZZZZ,ZZZZZ,beta] = GLM_VSDI(ZZ,[0.78 3.3 6.6],params.experiment.theoreticalSigs');
params.experiment.ZZ = ZZ; params.experiment.ZZZ = ZZZ;
mapGLM = postProcess(rshp([beta(end-(params.experiment.N-1):end,:)]'));
mapAOF = AvgOfFrms(rshp(ZZ));
mapTmax = Tmax(ZZ); 
mapCorr = Correlation2All(ZZ);
[mapTSCAboth] = genTSCA(ZZ,ZZZ);
mapTSCA = mapTSCAboth.withGLM;
mapNadav = MPT(ZZ);