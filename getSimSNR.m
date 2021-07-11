function [SNR] = getSimSNR()
    stim_time = 10:20:90;
    T = 1000; fs = 100;
    P = 1600;
    m = sqrt(P);
    noiseSig = [0.1 0.8 2 10];
    t = linspace(0,(T-1)/fs,T);
    signals = struct;
    noises = struct;
    SNR = zeros(1,length(noiseSig));
    for kk = 2:length(SNR) % iterate different noise sig
        thisNoiseSigSNR = zeros(1,100);
        for k=1:100
            %% construct signals
            [I,J] = ndgrid(1:m,1:m);
            locs = [0.4 0.4; 0.4 0.6; 0.5 0.5; 0.6 0.4; 0.6 0.6]; ind2plot = [1,5,9,13,17];
            DC = 0.02; r = 4;
            amp = (1.5-0.5).*rand(1,5)+0.5; % amplitude distributed U~[0.5 1.5]
            for i=1:5
                signals(i).time = amp(i).*(2.5.*normpdf(0:0.1:(T-1)/10,stim_time(i),1)); % normpdf with peak amplitude 1, times random amplitude
                signals(i).space = MinMaxNorm(imgaussfilt(double((I-m*locs(i,1)).^2+(J-m*locs(i,2)).^2<r^2),sqrt(r)))+DC; % what is larger than r^2 is 1 (white), what is smaller is 0 (black)
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
            thisNoiseSigSNR(k) = thisIterSNR(signals,noises,m,T);
        end
        SNR(kk) = mean(thisNoiseSigSNR);
    end
end
function SNR = thisIterSNR(signals,noises,m,T)
    ZSig = zeros(m*m,T); % preallocate memory
    ZNoise = zeros(m*m,T); % preallocate memory
    for i = 1:T
        ZSig(:,i) = reshape( ...
            signals(1).time(i)*signals(1).space+signals(2).time(i)*signals(2).space+...
            signals(3).time(i)*signals(3).space+signals(4).time(i)*signals(4).space+...
            signals(5).time(i)*signals(5).space,[],1);
        ZNoise(:,i) = reshape( ...
            noises(1).time(i)*noises(1).space+...
            noises(2).time(i)*noises(2).space+...
            noises(3).time(i)*noises(3).space,[],1);
    end
    ZSig = reshape(ZSig,40,40,[]);
    ZNoise = reshape(ZNoise,40,40,[]);
    ZSig = ZSig(0.3*m:0.7*m,0.3*m:0.7*m,[50:150 250:350 450:550 650:750 850:950]);
    ZNoise = ZNoise(0.3*m:0.7*m,0.3*m:0.7*m,[50:150 250:350 450:550 650:750 850:950]);
    Zsig = reshape(ZSig,[],size(ZSig,3))'; Znoise = reshape(ZNoise,[],size(ZSig,3))';
    SNR = 10*log10((mean(Zsig,2)'*mean(Zsig,2))/(mean(Znoise,2)'*mean(Znoise,2)));
%     SNR = 10*log10((ZSig(:)'*ZSig(:))/(ZNoise(:)'*ZNoise(:)));
end

