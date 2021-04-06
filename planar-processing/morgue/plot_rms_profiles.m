%% Plot RMS spatial profiles for error and uncertainty
clear;
% basedir='Z:\Planar_Uncertainty_work\Results\final_plots\07_31_2018\';
basedir='/scratch/shannon/a/aether/Projects/PIV_PTV_Uncertainty_Quantification/Planar_uncertainty_work/';
% casename='WS64';
% casename='WS32';
% casename='WS32\New_gradient\';
% casename='WS32\SOC_gradient\';
casename='WS64\SOC_gradient\';
% outdir=[basedir,'spatial_RMS_profiles\',casename,filesep];
outdir=fullfile(basedir,'spatial_RMS_profiles',casename);
if ~exist(outdir,'dir');
    mkdir(outdir);
end
savemat=0;
saveplot=0;
% keyboard;
%% 2003B
A1=load(fullfile(basedir,casename,'matfiles', '\'PivChal03B.mat'));
% Upmean=mean(A1.Up,3);
% Vpmean=mean(A1.Vp,3);
% figure;imagesc(A1.Xp(:),A1.Yp(:),Upmean);
Xp=A1.Xp/max(A1.Xp(:));
Yp=A1.Yp/max(A1.Yp(:));
velmagmean=mean((A1.Up.^2+A1.Vp.^2).^0.5,3);
velmagmean=velmagmean(end:-1:1,:);
figure;imagesc(Xp(:),Yp(:),velmagmean);

% error_prana=[A1.err_up;A1.err_vp];
% error_davis=[A1.err_ud;A1.err_vd];
% UMC=[A1.UMCx;A1.UMCy];
% UIM=[A1.UIMx;A1.UIMy];
% UCS=[A1.UCSx;A1.UCSy];
error_prana=A1.err_up;
error_davis=A1.err_ud;
UMC=A1.UMCx;
UIM=A1.UIMx;
UCS=A1.UCSx;

[yig,xig]=meshgrid(linspace(12,500,62),linspace(12,1524,190));

locx=192/2;% 768
locx1=locx-1;%764

EP=squeeze(error_prana(end:-1:1,locx,:));
ED=squeeze(error_davis(end:-1:1,locx1,:));
MC=squeeze(UMC(end:-1:1,locx,:));
IM=squeeze(UIM(end:-1:1,locx,:));
CS=squeeze(UCS(end:-1:1,locx1,:));

thresh=1;
[indx1,indy1]=find(abs(EP)>thresh);
[indx2,indy2]=find(abs(ED)>thresh);

for i=1:62
    k=1;
    l=1;
    for j=1:100
        if ismember(i,indx1) && ismember(j,indy1)
            continue;
        else
            tempstoreEP(k)=EP(i,j);
            tempstoreMC(k)=MC(i,j);
            tempstoreIM(k)=IM(i,j);
            k=k+1;
        end
        if ismember(i,indx2) && ismember(j,indy2)
            continue;
        else
            tempstoreED(l)=ED(i,j);
            tempstoreCS(l)=CS(i,j);
            l=l+1;
        end
    end
    RMSerrorprana(i)=rms(tempstoreEP);
    RMSerrordavis(i)=rms(tempstoreED);
    RMSUMC(i)=rms(tempstoreMC);
    RMSUIM(i)=rms(tempstoreIM);
    RMSUCS(i)=rms(tempstoreCS);
end
% yax=[8:8:504 8:8:504]';
% yax1=[8:8:502 8:8:502]';
yax=(8:8:503)';
yax1=(12:8:500)';
yax=yax/max(yax(:));
yax1=yax1/max(yax1(:));

xax=A1.Xp(1,:)/max(A1.Xp(1,:));
xpoint=xax(locx)*ones(length(yax),1);
hold on;plot(xpoint,yax,'k--');hold off;

lw=2;
figure;hold on;set(gcf,'DefaultLineLineWidth',lw);
% plot(RMSerrorprana,yax,'k-');
% plot(RMSerrordavis,yax1,'k--');
% plot(RMSUMC,yax,'r-');
% plot(RMSUIM,yax,'c-');
% plot(RMSUCS,yax1,'m-');
plot(yax,RMSerrorprana,'k-');
plot(yax1,RMSerrordavis,'k--');
plot(yax,RMSUMC,'r-');
plot(yax,RMSUIM,'c-');
plot(yax1,RMSUCS,'m-');
hold off;
axis([0 1 0 0.2]);

return;
if saveplot==1
    print(gcf,'-dpng',fullfile(outdir,'2003B_spatial_profile.png'),'-r360');
end
if savemat==1
    save(fullfile(outdir,'2003B_spatial_profilemat.mat'),'yax','yax1','RMSerrorprana','RMSerrordavis','RMSUMC','RMSUIM','RMSUCS','xpoint','yax','velmagmean','Xp','Yp');
end

%% 2005B
clear RMSerrorprana RMSerrordavis RMSUMC RMSUIM RMSUCS;
clear tempstoreEP tempstoreMC tempstoreIM tempstoreED tempstoreCS;
A1=load([basedir,casename,'\matfiles\PivChal05B.mat']);
% Upmean=mean(A1.Up,3);
% figure;imagesc(A1.Xp(:),A1.Yp(:),Upmean);
Xp=A1.Xp/max(A1.Xp(:));
Yp=A1.Yp/max(A1.Yp(:));
velmagmean=mean((A1.Up(:,:,1:5).^2+A1.Vp(:,:,1:5).^2).^0.5,3);
% velmagmean=velmagmean(end:-1:1,:);
figure;imagesc(Xp(:),Yp(:),velmagmean(:,:));

error_prana=A1.err_up;
error_davis=A1.err_ud;
UMC=A1.UMCx;
UIM=A1.UIMx;
UCS=A1.UCSx;

locx=44;% 712
locx1=44;%712

EP=squeeze(error_prana(:,locx,:));
ED=squeeze(error_davis(:,locx1,:));
MC=squeeze(UMC(:,locx,:));
IM=squeeze(UIM(:,locx,:));
CS=squeeze(UCS(:,locx1,:));

thresh=1;
[indx1,indy1]=find(abs(EP)>thresh);
[indx2,indy2]=find(abs(ED)>thresh);

for i=1:41
    k=1;
    l=1;
    for j=1:5
        if ismember(i,indx1) && ismember(j,indy1)
            continue;
        else
            tempstoreEP(k)=EP(i,j);
            tempstoreMC(k)=MC(i,j);
            tempstoreIM(k)=IM(i,j);
            k=k+1;
        end
        if ismember(i,indx2) && ismember(j,indy2)
            continue;
        else
            tempstoreED(l)=ED(i,j);
            tempstoreCS(l)=CS(i,j);
            l=l+1;
        end
    end
    RMSerrorprana(i)=rms(tempstoreEP);
    RMSerrordavis(i)=rms(tempstoreED);
    RMSUMC(i)=rms(tempstoreMC);
    RMSUIM(i)=rms(tempstoreIM);
    RMSUCS(i)=rms(tempstoreCS);
end
% yax=[8:8:504 8:8:504]';
% yax1=[8:8:502 8:8:502]';
yax=(24:16:664)';
yax1=(24:16:664)';
yax=yax/max(yax(:));
yax1=yax1/max(yax1(:));
xax=A1.Xp(1,:)/max(A1.Xp(1,:));
xpoint=xax(locx)*ones(length(yax),1);
hold on;plot(xpoint,yax,'k--');hold off;

lw=2;
figure;hold on;set(gcf,'DefaultLineLineWidth',lw);
plot(yax,RMSerrorprana,'k-');
plot(yax1,RMSerrordavis,'k--');
plot(yax,RMSUMC,'r-');
plot(yax,RMSUIM,'c-');
plot(yax1,RMSUCS,'m-');
hold off;
if saveplot==1
    print(gcf,'-dpng',fullfile(outdir,'2005B_spatial_profile.png'),'-r360');
end
if savemat==1
    save(fullfile(outdir,'2005B_spatial_profilemat.mat'),'yax','yax1','RMSerrorprana','RMSerrordavis','RMSUMC','RMSUIM','RMSUCS','xpoint','yax','velmagmean','Xp','Yp');
end

%% stagnation flow
clear RMSerrorprana RMSerrordavis RMSUMC RMSUIM RMSUCS;
clear tempstoreEP tempstoreMC tempstoreIM tempstoreED tempstoreCS;
A1=load([basedir,casename,'\matfiles\stagnation_flow.mat']);
% Upmean=mean(A1.Up,3);
% Vpmean=mean(A1.Vp,3);
% figure;imagesc(A1.Xp(:),A1.Yp(:),Upmean);
% figure;quiver(A1.Xp,A1.Yp,Upmean,Vpmean,5);
% figure;quiver(A1.Xp,A1.Yp,Upmean,Vpmean,5);
Xp=A1.Xp/max(A1.Xp(:));
Yp=A1.Yp/max(A1.Yp(:));
% Yp=Yp(end:-1:1,:);
velmagmean=mean((A1.Up.^2+A1.Vp.^2).^0.5,3);
% velmagmean=velmagmean(end:-1:1,:);
figure;imagesc(Xp(:),Yp(:),velmagmean);

error_prana=A1.err_up;
error_davis=A1.err_ud;
UMC=A1.UMCx;
UIM=A1.UIMx;
UCS=A1.UCSx;

locx=40;% 632
locx1=40;%632

EP=squeeze(error_prana(:,locx,:));
ED=squeeze(error_davis(:,locx1,:));
MC=squeeze(UMC(:,locx,:));
IM=squeeze(UIM(:,locx,:));
CS=squeeze(UCS(:,locx1,:));

thresh=1;
[indx1,indy1]=find(abs(EP)>thresh);
[indx2,indy2]=find(abs(ED)>thresh);

for i=1:64
    k=1;
    l=1;
    for j=1:99
        if ismember(i,indx1) && ismember(j,indy1)
            continue;
        else
            tempstoreEP(k)=EP(i,j);
            tempstoreMC(k)=MC(i,j);
            tempstoreIM(k)=IM(i,j);
            k=k+1;
        end
        if ismember(i,indx2) && ismember(j,indy2)
            continue;
        else
            tempstoreED(l)=ED(i,j);
            tempstoreCS(l)=CS(i,j);
            l=l+1;
        end
    end
    RMSerrorprana(i)=rms(tempstoreEP);
    RMSerrordavis(i)=rms(tempstoreED);
    RMSUMC(i)=rms(tempstoreMC);
    RMSUIM(i)=rms(tempstoreIM);
    RMSUCS(i)=rms(tempstoreCS);
end
% yax=[8:8:504 8:8:504]';
% yax1=[8:8:502 8:8:502]';
yax=(8:16:1016)';
yax1=(8:16:1016)';
yax=yax/max(yax(:));
yax1=yax1/max(yax1(:));
xax=A1.Xp(1,:)/max(A1.Xp(1,:));
xpoint=xax(locx)*ones(length(yax),1);
hold on;plot(xpoint,yax,'k--');hold off;

lw=2;
figure;hold on;set(gcf,'DefaultLineLineWidth',lw);
plot(yax,RMSerrorprana,'k-');
plot(yax1,RMSerrordavis,'k--');
plot(yax,RMSUMC,'r-');
plot(yax,RMSUIM,'c-');
plot(yax1,RMSUCS,'m-');
hold off;
axis([0 1 0 0.2])
if saveplot==1
    print(gcf,'-dpng',fullfile(outdir,'stagnation_spatial_profile.png'),'-r360');
end
if savemat==1
    save(fullfile(outdir,'stagnation_spatial_profilemat.mat'),'yax','yax1','RMSerrorprana','RMSerrordavis','RMSUMC','RMSUIM','RMSUCS','xpoint','yax','velmagmean','Xp','Yp');
end

%% Vortex Ring
clear RMSerrorprana RMSerrordavis RMSUMC RMSUIM RMSUCS;
clear tempstoreEP tempstoreMC tempstoreIM tempstoreED tempstoreCS;
A1=load([basedir,casename,'\matfiles\Vortex_Ring.mat']);
% Upmean=mean(A1.Up,3);
% figure;imagesc(A1.Xp(:),A1.Yp(:),Upmean);
Xp=(A1.Xp-min(A1.Xp(:)))/range(A1.Xp(:));
Yp=(A1.Yp-min(A1.Yp(:)))/range(A1.Yp(:));
velmagmean=mean((A1.Up.^2+A1.Vp.^2).^0.5,3);
velmagmean=velmagmean(end:-1:1,:);
figure;imagesc(Xp(:),Yp(:),velmagmean);
% figure;imagesc(Xp(:),Yp(:),Upmean);

error_prana=A1.err_up;
error_davis=A1.err_ud;
UMC=A1.UMCx;
UIM=A1.UIMx;
UCS=A1.UCSx;

locy=17;% 632
locy1=17;%632

EP=squeeze(error_prana(locy,:,:));
ED=squeeze(error_davis(locy,:,:));
MC=squeeze(UMC(locy,:,:));
IM=squeeze(UIM(locy,:,:));
CS=squeeze(UCS(locy,:,:));

thresh=1;
[indx1,indy1]=find(abs(EP)>thresh);
[indx2,indy2]=find(abs(ED)>thresh);

for i=1:80
    k=1;
    l=1;
    for j=1:50
        if ismember(i,indx1) && ismember(j,indy1)
            continue;
        else
            tempstoreEP(k)=EP(i,j);
            tempstoreMC(k)=MC(i,j);
            tempstoreIM(k)=IM(i,j);
            k=k+1;
        end
        if ismember(i,indx2) && ismember(j,indy2)
            continue;
        else
            tempstoreED(l)=ED(i,j);
            tempstoreCS(l)=CS(i,j);
            l=l+1;
        end
    end
    RMSerrorprana(i)=rms(tempstoreEP);
    RMSerrordavis(i)=rms(tempstoreED);
    RMSUMC(i)=rms(tempstoreMC);
    RMSUIM(i)=rms(tempstoreIM);
    RMSUCS(i)=rms(tempstoreCS);
end
% yax=[8:8:504 8:8:504]';
% yax1=[8:8:502 8:8:502]';
yax=(1:1:80)';
yax1=(1:1:80)';
yax=yax/max(yax(:));
yax1=yax1/max(yax1(:));

xax=Yp(:,1);
xpoint=xax(locy)*ones(length(yax),1);
hold on;plot(yax,xpoint,'k--');hold off;

lw=2;
figure;hold on;set(gcf,'DefaultLineLineWidth',lw);
plot(yax,RMSerrorprana,'k-');
plot(yax1,RMSerrordavis,'k--');
plot(yax,RMSUMC,'r-');
plot(yax,RMSUIM,'c-');
plot(yax1,RMSUCS,'m-');
hold off;
%axis([0 0.2 0 1016])
if saveplot==1
    print(gcf,'-dpng',fullfile(outdir,'vortexring_spatial_profile.png'),'-r360');
end
if savemat==1
    save(fullfile(outdir,'vortexring_spatial_profilemat.mat'),'yax','yax1','RMSerrorprana','RMSerrordavis','RMSUMC','RMSUIM','RMSUCS','xpoint','yax','velmagmean','Xp','Yp');
end

%% Jetdata
clear RMSerrorprana RMSerrordavis RMSUMC RMSUIM RMSUCS;
clear tempstoreEP tempstoreMC tempstoreIM tempstoreED tempstoreCS;
A1=load([basedir,casename,'\matfiles\Jetdata.mat']);
% Upmean=mean(A1.Up,3);
% figure;imagesc(A1.Xp(:),A1.Yp(:),Upmean);
Xp=(A1.Xp-min(A1.Xp(:)))/range(A1.Xp(:));
Yp=-(A1.Yp-254)/(10.2*7.2);%max(A1.Yp(:));
% Yp=max(A1.Yp(:));
velmagmean=mean((A1.Up.^2+A1.Vp.^2).^0.5,3);
% velmagmean=velmagmean(end:-1:1,:);
figure;imagesc(Xp(:),Yp(:),velmagmean);

error_prana=A1.err_up;
error_davis=A1.err_ud;
UMC=A1.UMCx;
UIM=A1.UIMx;
UCS=A1.UCSx;

locx=16;% 632
locx1=locx;%632

EP=squeeze(error_prana(:,locx,:));
ED=squeeze(error_davis(:,locx1,:));
MC=squeeze(UMC(:,locx,:));
IM=squeeze(UIM(:,locx,:));
CS=squeeze(UCS(:,locx1,:));

% RMSerrorprana=rms(EP,2);
% RMSerrordavis=rms(ED,2);
% RMSUMC=rms(MC,2);
% RMSUIM=rms(IM,2);
% RMSUCS=rms(CS,2);

thresh=1;
[indx1,indy1]=find(abs(EP)>thresh);
[indx2,indy2]=find(abs(ED)>thresh);

for i=1:31
    k=1;
    l=1;
    for j=1:495
        if ismember(i,indx1) && ismember(j,indy1)
            continue;
        else
            tempstoreEP(k)=EP(i,j);
            tempstoreMC(k)=MC(i,j);
            tempstoreIM(k)=IM(i,j);
            k=k+1;
        end
        if ismember(i,indx2) && ismember(j,indy2)
            continue;
        else
            tempstoreED(l)=ED(i,j);
            tempstoreCS(l)=CS(i,j);
            l=l+1;
        end
    end
    RMSerrorprana(i)=rms(tempstoreEP);
    RMSerrordavis(i)=rms(tempstoreED);
    RMSUMC(i)=rms(tempstoreMC);
    RMSUIM(i)=rms(tempstoreIM);
    RMSUCS(i)=rms(tempstoreCS);
end

% yax=[8:8:504 8:8:504]';
% yax1=[8:8:502 8:8:502]';
% yax=(1:1:31)';
% yax1=(1:1:31)';
Yd=A1.Yd(:,1);
yax=-((Yd-254)./(10.2*7.2));
yax1=yax;
xax=Xp;%A1.Xp(1,:)/max(A1.Xp(1,:));
xpoint=xax(1,locx)*ones(length(yax),1);
hold on;plot(xpoint,yax,'k--');hold off;

lw=2;
figure;hold on;set(gcf,'DefaultLineLineWidth',lw);
plot(yax,RMSerrorprana,'k-');
plot(yax1,RMSerrordavis,'k--');
plot(yax,RMSUMC,'r-');
plot(yax,RMSUIM,'c-');
plot(yax1,RMSUCS,'m-');
hold off;
 axis([-1 1 0 0.2])
if saveplot==1
    print(gcf,'-dpng',fullfile(outdir,'jet_spatial_profile.png'),'-r360');
end
if savemat==1
    save(fullfile(outdir,'jet_spatial_profilemat.mat'),'yax','yax1','RMSerrorprana','RMSerrordavis','RMSUMC','RMSUIM','RMSUCS','xpoint','yax','velmagmean','Xp','Yp');
end
