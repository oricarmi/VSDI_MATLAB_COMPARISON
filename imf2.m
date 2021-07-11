function im_out=imf2(im,alphaData)
%% show scaled image
global sz rot ump

if size(im,2)==1
    im_out=rot90(rshp(im),rot);
else
    im_out=rot90(rshp(im) ,rot);
end
if nargin<2
    imagesc(ump.*[0:sz(2)-1]./1000, ump.*[0:sz(1)-1]./1000, im_out);
else
    imagesc(ump.*[0:sz(2)-1]./1000, ump.*[0:sz(1)-1]./1000, im_out,'alphadata',im_out>alphaData);
end
xlabel('mm');ylabel('mm');
end
