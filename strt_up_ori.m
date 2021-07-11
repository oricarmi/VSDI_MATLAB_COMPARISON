function [cf1 cfn trsh0]=strt_up(fname, n)

global bsln fs sz ump rot fgn brn frq cmap lgn scl lgn00 fnm c_f vms lgnd vc prms

vc(1)=0;
whitebg([1 1 1]);
if strcmp(fname, fnm)==0; 
    load(fname); fnm=fname; 
    vms.mltpl=v_mltpl;
    vms.sngl=v_sngl;
    vms.ft=v_ft;
end 

% lgn=lgnd{n};

brn0=brn; scl=[]; rot=0; fs=100; ump=20;

if exist('msz')==1; sz=msz(n,:); else; sz=[270 327]; end

if round(mean(size(brn)==sz))~=1; brn=zeros(sz); else; brn=brn0; end; figure; img(brn); colormap Gray;

if exist('frqm')==1; frq=frqm(n); 
elseif isempty(prms)==0; frq=prms(n).frq; 
else; frq=1; end

trsh0=std(c_f{n},0,2);

for k0=1:size(c_f,1);
    if isempty(c_f{k0,n})==0; cf1(k0,1)=c_f(k0,n); l(k0)=size(cf1{k0},2); cfn{k0,1}=c_f{k0,n}./repmat(trsh0,1, size(c_f{k0,n},2)); end;
end

if std(l)~=0
    lmx=max(l);
   for k0=1:length(cf1);
       if l(k0)<lmx; cf1{k0}(:,l(k0)+1: lmx)=0; cfn{k0}(:,l(k0)+1: lmx)=0; end
   end
    
end
    



