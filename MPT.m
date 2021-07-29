function [mapNadav] = MPT(ZZ)
    global bsln fs sz ump rot fgn brn brn0 frq cmap lgn scl lgn00 fnm c_f vms plt_on pt pc vc xs prms cfn cfn0 basis params
    % set all the user defined thresholds
    vc(2)=0; xs=[0.3 1]; %% remove limits and title from plots if vc(2)=1 and show scale xs=[x y]
    vc(3)=1; %show reject/anti wave if vc(3)=1;
    plt_on=0; cnd=[];  s_flt=[params.post.gaussfltSTD 2]; t_flt = [];
    fgn=100; scl=[]; cmap=colormap(jet);
    if params.experiment.what<10
        Until = params.experiment.what;
    else
        Until = floor(params.experiment.what/10);
    end
    cfn2 = cell(Until,1);
    for i=1:Until
        cfn2{i} = ZZ(:,(i-1)*params.experiment.T1+1:i*params.experiment.T1);
    end
    [mxc rxc]=xp_mp(cfn, cnd, params.Nadav.p, params.Nadav.x, params.Nadav.t_lmts, params.Nadav.settle,  t_flt, s_flt); % perform the thresholding
    mapNadav = postProcess(rshp(mxc(:,end-Until+1:end))); % post process the resulting map 
end 

