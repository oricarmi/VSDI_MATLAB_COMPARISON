function [D] = calcDistances(maps)
% input: maps. output: (symmetric) matrix of distnaces D between conditions
r = zeros(1,size(maps,3));
c = zeros(1,size(maps,3));
for i=1:size(maps,3)
    [~,r(i),c(i)] = calcCOM(maps(:,:,i));
end
D = zeros(size(maps,3),size(maps,3));
for i=1:length(D)
    for j=i:length(D)
        if i~=j % if not in same index (diagonal of D is zero)
            D(i,j) = sqrt((r(j)-r(i))^2+(c(j)-c(i))^2);
            D(j,i) = D(i,j);
        end
    end
end
    
end

