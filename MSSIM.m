function mssim = MSSIM(orig,rec,win)
% org: Original image
% sup: Super-resolution image
% win: Size of the window
[m,n,p]= size(orig);
c1= 0.0001;
c2= 0.0001;
sz= floor(win/2);
totalErr= [];
for k= 1:p % iterate rgb
    for i= sz+1:m-sz % iterate rows
        for j= sz+1:n-sz % iterate cols
            six= orig(i-sz:i+sz,j-sz:j+sz,k); % define Sx (original image)
            six= six(:);
            siy= rec(i-sz:i+sz,j-sz:j+sz,k); % define Sy (reconstructed image)
            siy= siy(:);
            mux= mean(six); % mu of orig
            muy= mean(siy); % mu of recon
            sigx= std(six); % std of orig
            sigy= std(siy); % std of recon
            sigxy= (six-mux)'*(siy-muy)/(numel(six)-1); % covariance of orig and recon
            err= ((2*mux*muy+c1)*(2*sigxy+c2))/((mux^2+muy^2+c1)*(sigx^2+sigy^2+c2)); % this window error
            totalErr = [totalErr err]; % append to vector
        end
    end
end
mssim = sum(totalErr)/(numel(totalErr)); % mean of all windows