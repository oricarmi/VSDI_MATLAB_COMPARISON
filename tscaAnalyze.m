function [] = tscaAnalyze(tscaStruct,numComp,Title,largeComp,T)
% Show results analysis of algorithm
global brn0 brn
if nargin<2 % default is 25 components (5x5 subplot)
    numComp = 5;
    Title = [];
    largeComp = [];
    T=100;
elseif nargin<3
    Title = [];
    largeComp = [];
    T=100;
elseif nargin<4
    largeComp = [];
    T=100;
end
if isempty(numComp)
    numComp = 5;
end
figure;suptitle([Title ' spatial components, gammas: ' num2str(tscaStruct.gammas)]);
for i=1:numComp^2
    subplot(numComp,numComp,i);
    try
        imagesc(reshape(tscaStruct.components(:,i),size(brn0,1),[]));colormap(jet);
    catch
        try
            imagesc(reshape(tscaStruct.components(:,i),size(brn,1),[]));colormap(jet);
        catch
            imagesc(abs(reshape(tscaStruct.components(:,i),sqrt(length(tscaStruct.components(:,i))),[])));colormap(gray);
        end
    end
    title(['\phi = ' num2str(i)]);
end
t= linspace(0,10,T); 
figure;suptitle([Title ' time course, gammas: ' num2str(tscaStruct.gammas)]);
for i=1:numComp^2
    subplot(numComp,numComp,i); hold on;
    plot(abs(tscaStruct.projected(i,:))); xlabel('time [sec]');ylabel('amplitude');
    title(['$\phi$ = ' num2str(i)]);
    for j=1:10
        xline(i);
    end
end
figure;suptitle([Title ' - eigen-values']); eigenVals = sort(abs(diag(tscaStruct.D)),'desc');
plot(eigenVals(1:numComp^2),'*-');xlabel('component #');ylabel('eigen value');
if ~isempty(largeComp) && largeComp>0 % show enlarged component of interest
    figure;suptitle([Title ', component:' num2str(largeComp)]);
    try
        imagesc(reshape(tscaStruct.components(:,largeComp),size(brn0,1),[]));colormap(jet);colorbar('southoutside');
    catch
        imagesc(reshape(tscaStruct.components(:,largeComp),size(brn,1),[]));colormap(jet);colorbar('southoutside');
    end
end
figure;suptitle('autocorrelation matrices - 1st is signal')
for i=1:size(tscaStruct.C,3)
    subplot(1,size(tscaStruct.C,3),i);
    imagesc(tscaStruct.C(:,:,i));colorbar('southoutside');
end
