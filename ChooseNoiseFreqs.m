function [] = ChooseNoiseFreqs(Z)
global fs cfn cfn0 brn brn0 sz sz0
more = 1;
cfn0=cfn; sz0=sz; brn0=brn; % save original condition
while more
    cfn=pk_zn2(cfn,2);
    Z = Z - mean(Z,2);
    % [ppx,f] = periodogram(mean(Z),[],[],fs);
    [ppx,f] = pmtm(mean(Z),2,pow2(nextpow2(length(mean(Z)))),fs); % calculate spectrum
    figure; % show time course and spectrum
    subplot(1,2,1);plot(mean(Z)); xlabel('frame'); ylabel('amp'); title('time course of average signal');
    subplot(1,2,2); plot(f,ppx); xlabel('frequency [Hz]');ylabel('amp');title('spectrum of average signal');
    disp('insert to csv file the noise frequencies');
    cfn=cfn0; sz=sz0; brn=brn0; %restore to original condition
    more = input('choose more noise freqs (from differenet region)? [0/1]');
end

