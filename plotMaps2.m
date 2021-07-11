function plotMaps2(map,Ttle,CA)
% plot maps (3d matrix, 3rd dimension is differnet condition)
global params lgn brn ump
    if nargin<2
        CA = 0;
        Ttle = 'unknown';
    elseif nargin<3
        CA = 0;
    end
    if iscell(map)
        map = cat(3,map{:});
    end 
    for i=1:size(params.experiment.optimalMaps.orig,3) % calc Center Of Mass of each image
        [~,r(i),c(i)] = calcCOM(map(:,:,i));
    end
    figure("name",sprintf('indMaps %s', Ttle)); 
%     suptitle([Ttle]);% ", caxis: " num2str(CA)]);
    cAxis = [prctile(map(:),1) prctile(map(:),100)];
    tempbrn = repmat(reshape((brn-min(brn))./(max(brn)-min(brn)),[],1),1,1,3);
    for i=1:params.experiment.N
        switch params.experiment.what
            case 8
                index2plot = [11,3,15,7,5,13,1,9]; % 8 locs
                subplot(3,5,index2plot(i));
            case 9 % 9 locs
                index2plot = [5,7,3,8,2,9,1,6,4];
                subplot(3,3,index2plot(i));
            case 4 % 4 locs 
                index2plot = [3,2,1,4];
                subplot(2,2,index2plot(i));
            case 52 % 5 locs 2[Hz]
                index2plot = [7,3,1,9,5];
                subplot(3,3,index2plot(i));
            case 92 % 9 moving bars 2[Hz]
                subplot(3,3,i)
        end
        imf2(rshp(tempbrn)); hold on;
        imf2(map(:,:,i),prctile(reshape(map(:,:,i),[],1),90));
        imf2(params.experiment.optimalMaps.orig(:,:,i),prctile(reshape(params.experiment.optimalMaps.orig(:,:,i),[],1),98));
%         plot(ump/1000*(c(i)-1),ump/1000*(r(i)-1),'*r','MarkerSize',2);
        colormap('parula'); title(lgn(i+2,:));
        if CA
            caxis(cAxis);
        end
    end
    
end

