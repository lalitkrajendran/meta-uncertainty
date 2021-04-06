function compute_error_uncertainty(casename,inputdir,outputdir,savemat)
%diameter index
% dn=6;%autod
dn=26;%crossd
% addpath Y:\Projects\TurbulentFlame\analysis\src\readimx_v2.0_win64\readimx;
% addpath Y:\Projects\turbulent_flame\analysis\src\readimx_v2.0_win64\readimx
addpath Z:\Planar_Uncertainty_work\codes\readimx_v2.0_win64\readimx;
if strcmp(casename{1},'PivChal03B')
    fprintf('\n PIV Challenge 2003B case uncertainty calculation \n');
    %% Load True Solution
    truesol=load(fullfile(inputdir.Pivchal03B.truesoldir,[inputdir.Pivchal03B.truesolbase,'.mat']));

    %% Set Parameters
    N=100;
    fstart=1;
%     foffset=0;
    sxp=63;syp=191;
    sxd=62;syd=190;
    Up=zeros(sxp,syp,N);
    Vp=zeros(sxp,syp,N);
    Ud=zeros(sxd,syd,N);
    Vd=zeros(sxd,syd,N);
    
    Utrue=zeros(sxp,syp,N);
    Vtrue=zeros(sxp,syp,N);
    Utrued=zeros(sxd,syd,N);
    Vtrued=zeros(sxd,syd,N);
    
    err_up=zeros(sxp,syp,N);
    err_vp=zeros(sxp,syp,N);
    err_ud=zeros(sxd,syd,N);
    err_vd=zeros(sxd,syd,N);
    
    UMCx=zeros(sxp,syp,N);
    UMCy=zeros(sxp,syp,N);
    UIMx=zeros(sxp,syp,N);
    UIMy=zeros(sxp,syp,N);
    UCSx=zeros(sxd,syd,N);
    UCSy=zeros(sxd,syd,N);
    Udiffstore=zeros(sxp,syp,N);
    Vdiffstore=zeros(sxp,syp,N);
    %% Loop through N fields
    for i=fstart:N
        
        %% Load prana solution
        pranasol=load(fullfile(inputdir.Pivchal03B.pranasol,[inputdir.Pivchal03B.pranabaseMC, num2str(2*i-1,'%03.0f'),'.mat']));
        imsol=load(fullfile(inputdir.Pivchal03B.IMdir,[inputdir.Pivchal03B.pranabaseIM, num2str(2*i-1,'%03.0f'),'.mat']));
%         keyboard;
        Up(:,:,i)=pranasol.U(end:-1:1,:,1);
        Vp(:,:,i)=pranasol.V(end:-1:1,:,1);
        if i==fstart
            Xp=pranasol.X;
            Yp=pranasol.Y;
        end
        %% Load Davis Solution
        davissol=readimx(fullfile(inputdir.Pivchal03B.davisol,[inputdir.Pivchal03B.davisbase,num2str(i,'%05.0f'),'.vc7']));
        Utemp=cell2mat(davissol.Frames{1,1}.Components{1,1}.Planes)';
        Vtemp=cell2mat(davissol.Frames{1,1}.Components{2,1}.Planes)';
        Ud(:,:,i)=Utemp(2:end-1,2:end-1);
        Vd(:,:,i)=-Vtemp(2:end-1,2:end-1);
        %Davis X and Y grid 
        if i==fstart
            [D] = create2DVec(davissol.Frames{1,1});
            Xd=D.X(2:end-1,2:end-1);
            Yd=D.Y(2:end-1,2:end-1);
            clear D;
            Xt=truesol.x;
            Yt=truesol.y(:,end:-1:1);
            
            [yig,xig]=meshgrid(linspace(12,500,62),linspace(12,1524,190));
            
        end
        
        %% Get True Solution
        Utrue(:,:,i)=truesol.u(:,:,i)';
        Vtrue(:,:,i)=truesol.v(:,:,i)';
        %Interpolated true solution onto Davis Grid
        Utrued(:,:,i)=interp2(Xt',Yt',Utrue(:,:,i),xig',yig','linear',0);
        Vtrued(:,:,i)=interp2(Xt',Yt',Vtrue(:,:,i),xig',yig','linear',0);
        %% Calculate Error
        err_up(:,:,i)=(Up(:,:,i)-Utrue(:,:,i));
        err_vp(:,:,i)=(Vp(:,:,i)-Vtrue(:,:,i));
        err_ud(:,:,i)=(Ud(:,:,i)-Utrued(:,:,i));
        err_vd(:,:,i)=(Vd(:,:,i)-Vtrued(:,:,i));
        
        %% Get MC uncertainty estimates
        Ixx=pranasol.ixx(end:-1:1,:);
        Iyy=pranasol.iyy(end:-1:1,:);
        scaling=pranasol.mi(end:-1:1,:);
        biasx=reshape(pranasol.uncertainty(:,15),sxp,syp);
        biasy=reshape(pranasol.uncertainty(:,16),sxp,syp);
        Autod=reshape(pranasol.uncertainty(:,dn),sxp,syp);
        biasx=biasx(end:-1:1,:);
        biasy=biasy(end:-1:1,:);
        Autod=Autod(end:-1:1,:);

        %Gradient correction
        Udiff=socdiff(Up(:,:,i),8,1);
        Vdiff=socdiff(Vp(:,:,i),8,2);
%         [dudx ]=gradient_compact_rich(N,dx,dir)
%         Udiff=gradient_compact_rich(Up(:,:,i),8,1);
%         Vdiff=gradient_compact_rich(Vp(:,:,i),8,2);
        Udiffstore(:,:,i)=Udiff;
        Vdiffstore(:,:,i)=Vdiff;
        
        Ixxt= real(sqrt(Ixx.^2 - (Autod.^2/16).*(Udiff).^2 ));
        Iyyt= real(sqrt(Iyy.^2 - (Autod.^2/16).*(Vdiff).^2 ));
        %scaling and bias correction
        UMCx(:,:,i)=sqrt(biasx.^2+(Ixxt.^2)./scaling);
        UMCy(:,:,i)=sqrt(biasy.^2+(Iyyt.^2)./scaling);
        
        
        %% Get IM uncertainty estimates
        UIMx(:,:,i)=imsol.imx(end:-1:1,:);
        UIMy(:,:,i)=imsol.imy(end:-1:1,:);
        
        %% Get CS uncertainty estimates
        Unrxtemp=cell2mat(davissol.Frames{1,1}.Components{22,1}.Planes)';
        Unrytemp=cell2mat(davissol.Frames{1,1}.Components{23,1}.Planes)';
        Unbxtemp=cell2mat(davissol.Frames{1,1}.Components{25,1}.Planes)';
        Unbytemp=cell2mat(davissol.Frames{1,1}.Components{26,1}.Planes)';
                
        Unrx=Unrxtemp(2:end-1,2:end-1);
        Unry=Unrytemp(2:end-1,2:end-1);
        Unbx=Unbxtemp(2:end-1,2:end-1);
        Unby=Unbytemp(2:end-1,2:end-1);
        
        UCSx(:,:,i)=sqrt(Unrx.^2+Unbx.^2);
        UCSy(:,:,i)=sqrt(Unry.^2+Unby.^2);
        
   
    end
    Xd=xig';
    Yd=yig';
    %% Save error and uncertainty values
    if savemat==1
        save(fullfile(outputdir{1},'PivChal03B.mat'),'err_up','err_vp','err_ud','err_vd','UMCx','UMCy','UIMx','UIMy','UCSx','UCSy','Up','Vp','Ud','Vd','Xp','Yp','Xt','Yt','Xd','Yd','Udiff','Vdiff');
    end
    clearvars -except casename inputdir outputdir savemat dn;
    fprintf('\nDone...\n');
    
end
%%
if strcmp(casename{2},'PivChal05B')
    fprintf('\n PIV Challenge 2005B case uncertainty calculation \n');
    %% Load True Solution
    truesol=load(fullfile(inputdir.Pivchal05B.truesoldir,[inputdir.Pivchal05B.truesolbase,'.mat']));
    
    %% Set Parameters
    N=6;
    fstart=1;
%     foffset=0;
    f=[9 29 49 69 89 109];

    sxp=41;syp=88;
    sxd=41;syd=88;
    Up=zeros(sxp,syp,N);
    Vp=zeros(sxp,syp,N);
    Ud=zeros(sxd,syd,N);
    Vd=zeros(sxd,syd,N);
    
    Utrue=zeros(sxp,syp,N);
    Vtrue=zeros(sxp,syp,N);
%     Utrued=zeros(sxd,syd,N);
%     Vtrued=zeros(sxd,syd,N);
    
    err_up=zeros(sxp,syp,N);
    err_vp=zeros(sxp,syp,N);
    err_ud=zeros(sxd,syd,N);
    err_vd=zeros(sxd,syd,N);
    
    UMCx=zeros(sxp,syp,N);
    UMCy=zeros(sxp,syp,N);
    UIMx=zeros(sxp,syp,N);
    UIMy=zeros(sxp,syp,N);
    UCSx=zeros(sxd,syd,N);
    UCSy=zeros(sxd,syd,N);
    Udiffstore=zeros(sxp,syp,N);
    Vdiffstore=zeros(sxp,syp,N);
    
    %% Loop through N fields
    for i=fstart:N
        
        %% Load prana solution
        pranasol=load(fullfile(inputdir.Pivchal05B.pranasol,[inputdir.Pivchal05B.pranabaseMC, num2str(f(i),'%03.0f'),'.mat']));
        imsol=load(fullfile(inputdir.Pivchal05B.IMdir,[inputdir.Pivchal05B.pranabaseIM, num2str(f(i),'%03.0f'),'.mat']));
%         keyboard;
        Up(:,:,i)=pranasol.U./2;
        Vp(:,:,i)=pranasol.V./2;
        if i==fstart
            Xp=pranasol.X;
            Yp=pranasol.Y;
        end
        %% Load Davis Solution
        davissol=readimx(fullfile(inputdir.Pivchal05B.davisol,[inputdir.Pivchal05B.davisbase,num2str(i,'%05.0f'),'.vc7']));
%         keyboard;
        Utemp=cell2mat(davissol.Frames{1,1}.Components{1,1}.Planes)';
        Vtemp=cell2mat(davissol.Frames{1,1}.Components{2,1}.Planes)';
        Ud(:,:,i)=Utemp(end-1:-1:2,2:end-1)./2;
        Vd(:,:,i)=-Vtemp(end-1:-1:2,2:end-1)./2;
        %Davis X and Y grid 
        if i==fstart
            [D] = create2DVec(davissol.Frames{1,1});
            Xd=D.X(end-1:-1:2,2:end-1);
            Yd=D.Y(end-1:-1:2,2:end-1);
            clear D;
            Xt=truesol.Xt;
            Yt=truesol.Yt;
        end
        
        %% Get True Solution
        Utrue(:,:,i)=truesol.Ut(:,:,i);
        Vtrue(:,:,i)=truesol.Vt(:,:,i);
                
        %% Calculate Error
        err_up(:,:,i)=(Up(:,:,i)-Utrue(:,:,i));
        err_vp(:,:,i)=(Vp(:,:,i)-Vtrue(:,:,i));
        err_ud(:,:,i)=(Ud(:,:,i)-Utrue(:,:,i));
        err_vd(:,:,i)=(Vd(:,:,i)-Vtrue(:,:,i));
        
        %% Get MC uncertainty estimates
        Ixx=pranasol.ixx./2;
        Iyy=pranasol.iyy./2;
        scaling=pranasol.mi;
        biasx=reshape(pranasol.uncertainty(:,15),sxp,syp)./2;
        biasy=reshape(pranasol.uncertainty(:,16),sxp,syp)./2;
        Autod=reshape(pranasol.uncertainty(:,dn),sxp,syp);
       
        %Gradient correction
        Udiff=socdiff(Up(:,:,i),8,1);
        Vdiff=socdiff(Vp(:,:,i),8,2);
%         Udiff=gradient_compact_rich(Up(:,:,i),8,1);
%         Vdiff=gradient_compact_rich(Vp(:,:,i),8,2);
        Udiffstore(:,:,i)=Udiff;
        Vdiffstore(:,:,i)=Vdiff;
        
        Ixxt= real(sqrt(Ixx.^2 - (Autod.^2/16).*(Udiff).^2 ));
        Iyyt= real(sqrt(Iyy.^2 - (Autod.^2/16).*(Vdiff).^2 ));
        %scaling and bias correction
        UMCx(:,:,i)=sqrt(biasx.^2+(Ixxt.^2)./scaling);
        UMCy(:,:,i)=sqrt(biasy.^2+(Iyyt.^2)./scaling);
        
        
        %% Get IM uncertainty estimates
        UIMx(:,:,i)=imsol.imx./2;
        UIMy(:,:,i)=imsol.imy./2;
        
        %% Get CS uncertainty estimates
        Unrxtemp=cell2mat(davissol.Frames{1,1}.Components{22,1}.Planes)';
        Unrytemp=cell2mat(davissol.Frames{1,1}.Components{23,1}.Planes)';
        Unbxtemp=cell2mat(davissol.Frames{1,1}.Components{25,1}.Planes)';
        Unbytemp=cell2mat(davissol.Frames{1,1}.Components{26,1}.Planes)';
                
        Unrx=Unrxtemp(end-1:-1:2,2:end-1)./2;
        Unry=Unrytemp(end-1:-1:2,2:end-1)./2;
        Unbx=Unbxtemp(end-1:-1:2,2:end-1)./2;
        Unby=Unbytemp(end-1:-1:2,2:end-1)./2;
        
        UCSx(:,:,i)=sqrt(Unrx.^2+Unbx.^2);
        UCSy(:,:,i)=sqrt(Unry.^2+Unby.^2);
        
   
    end
    %% Save error and uncertainty values
    if savemat==1
        save(fullfile(outputdir{1},'PivChal05B.mat'),'err_up','err_vp','err_ud','err_vd','UMCx','UMCy','UIMx','UIMy','UCSx','UCSy','Up','Vp','Ud','Vd','Xp','Yp','Xt','Yt','Xd','Yd','Udiff','Vdiff');
    end
    clearvars -except casename inputdir outputdir savemat dn;
    fprintf('\nDone...\n');
    
end

%%
if strcmp(casename{3},'stagnation_flow')
    fprintf('\n stagnation_flow case uncertainty calculation \n');
    %% Load True Solution
    truesol=load(fullfile(inputdir.stagflow.truesoldir,[inputdir.stagflow.truesolbase,'.mat']));
    
    %% Set Parameters
    N=99;
    fstart=1;

    sxp=64;syp=80;
    sxd=64;syd=80;
    Up=zeros(sxp,syp,N);
    Vp=zeros(sxp,syp,N);
    Ud=zeros(sxd,syd,N);
    Vd=zeros(sxd,syd,N);
    
    Utrue=zeros(sxp,syp,N);
    Vtrue=zeros(sxp,syp,N);
%     Utrued=zeros(sxd,syd,N);
%     Vtrued=zeros(sxd,syd,N);
    
    err_up=zeros(sxp,syp,N);
    err_vp=zeros(sxp,syp,N);
    err_ud=zeros(sxd,syd,N);
    err_vd=zeros(sxd,syd,N);
    
    UMCx=zeros(sxp,syp,N);
    UMCy=zeros(sxp,syp,N);
    UIMx=zeros(sxp,syp,N);
    UIMy=zeros(sxp,syp,N);
    UCSx=zeros(sxd,syd,N);
    UCSy=zeros(sxd,syd,N);
    Udiffstore=zeros(sxp,syp,N);
    Vdiffstore=zeros(sxp,syp,N);
    
    %% Loop through N fields
    for i=fstart:N
        
        %% Load prana solution
        pranasol=load(fullfile(inputdir.stagflow.pranasol,[inputdir.stagflow.pranabaseMC, num2str(i,'%06.0f'),'.mat']));
        imsol=load(fullfile(inputdir.stagflow.IMdir,[inputdir.stagflow.pranabaseIM, num2str(i,'%06.0f'),'.mat']));
%         keyboard;
        Up(:,:,i)=pranasol.U;
        Vp(:,:,i)=pranasol.V;
        if i==fstart
            Xp=pranasol.X;
            Yp=pranasol.Y;
        end
        %% Load Davis Solution
        davissol=readimx(fullfile(inputdir.stagflow.davisol,[inputdir.stagflow.davisbase,num2str(i,'%05.0f'),'.vc7']));
%         keyboard;
        Utemp=cell2mat(davissol.Frames{1,1}.Components{1,1}.Planes)';
        Vtemp=cell2mat(davissol.Frames{1,1}.Components{2,1}.Planes)';
        
        Ud(:,:,i)=Utemp(end:-1:1,:);
        Vd(:,:,i)=-Vtemp(end:-1:1,:);
        %Davis X and Y grid 
        if i==fstart
            [D] = create2DVec(davissol.Frames{1,1});
            Xd=D.X(:,:)';
%             Xd=Xd(:,end:-1:1);
            Yd=D.Y(end:-1:1,:)';
            clear D;
%             Xt=truesol.Xt;
%             Yt=truesol.Yt;
        end
        
        %% Get True Solution
        Utrue(:,:,i)=truesol.Ut;
        Vtrue(:,:,i)=truesol.Vt;
        if i==fstart
            unfitx=truesol.unfitx;
            unfity=truesol.unfity;
        end
        %% Calculate Error
        err_up(:,:,i)=(Up(:,:,i)-Utrue(:,:,i));
        err_vp(:,:,i)=(Vp(:,:,i)-Vtrue(:,:,i));
        err_ud(:,:,i)=(Ud(:,:,i)-Utrue(:,:,i));
        err_vd(:,:,i)=(Vd(:,:,i)-Vtrue(:,:,i));
        
        %% Get MC uncertainty estimates
        Ixx=pranasol.ixx;
        Iyy=pranasol.iyy;
        scaling=pranasol.mi;
        biasx=reshape(pranasol.uncertainty(:,15),sxp,syp);
        biasy=reshape(pranasol.uncertainty(:,16),sxp,syp);
        Autod=reshape(pranasol.uncertainty(:,dn),sxp,syp);
       
        %Gradient correction
        Udiff=socdiff(Up(:,:,i),16,1);
        Vdiff=socdiff(Vp(:,:,i),16,2);
%         Udiff=gradient_compact_rich(Up(:,:,i),16,1);
%         Vdiff=gradient_compact_rich(Vp(:,:,i),16,2);
        Udiffstore(:,:,i)=Udiff;
        Vdiffstore(:,:,i)=Vdiff;
        
        Ixxt= real(sqrt(Ixx.^2 - (Autod.^2/16).*(Udiff).^2 ));
        Iyyt= real(sqrt(Iyy.^2 - (Autod.^2/16).*(Vdiff).^2 ));
        %scaling and bias correction
        UMCx(:,:,i)=sqrt(biasx.^2+(Ixxt.^2)./scaling +unfitx.^2);
        UMCy(:,:,i)=sqrt(biasy.^2+(Iyyt.^2)./scaling +unfity.^2);
        
        
        %% Get IM uncertainty estimates
        UIMx(:,:,i)=sqrt((imsol.imx).^2 +unfitx.^2);
        UIMy(:,:,i)=sqrt((imsol.imy).^2 +unfity.^2);
        
        %% Get CS uncertainty estimates
        Unrxtemp=cell2mat(davissol.Frames{1,1}.Components{22,1}.Planes)';
        Unrytemp=cell2mat(davissol.Frames{1,1}.Components{23,1}.Planes)';
        Unbxtemp=cell2mat(davissol.Frames{1,1}.Components{25,1}.Planes)';
        Unbytemp=cell2mat(davissol.Frames{1,1}.Components{26,1}.Planes)';
                
        Unrx=Unrxtemp(end:-1:1,:);
        Unry=Unrytemp(end:-1:1,:);
        Unbx=Unbxtemp(end:-1:1,:);
        Unby=Unbytemp(end:-1:1,:);
        
        UCSx(:,:,i)=sqrt(Unrx.^2+Unbx.^2 +unfitx.^2);
        UCSy(:,:,i)=sqrt(Unry.^2+Unby.^2 +unfity.^2);
%         keyboard;
   
    end
    Xt=Xp;
    Yt=Yp;
    %% Save error and uncertainty values
    if savemat==1
        save(fullfile(outputdir{1},'stagnation_flow.mat'),'err_up','err_vp','err_ud','err_vd','UMCx','UMCy','UIMx','UIMy','UCSx','UCSy','Up','Vp','Ud','Vd','Xp','Yp','Xt','Yt','Xd','Yd','Udiff','Vdiff');
    end
    clearvars -except casename inputdir outputdir savemat dn;
    fprintf('\nDone...\n');
    
end

%%
if strcmp(casename{4},'Vortex_Ring')
    fprintf('\n Vortex_Ring case uncertainty calculation \n');
    %% Load True Solution
    truesol=load(fullfile(inputdir.vortex.truesoldir,[inputdir.vortex.truesolbase,'.mat']));
    
    %% Set Parameters
    N=50;
    fstart=1;
    foffsetprana=24;
    foffsetdavis=1;
    
    spf=1;
    xscale=1/11.9690;%stereojob.datasave.scaling.wil.xscale;
    yscale=1/11.9653;%stereojob.datasave.scaling.wil.yscale;
    xorigin=576.8623;
    yorigin=395.8564;
    
    %true solution grid
    [Xt,Yt]=meshgrid(-34:0.45:20,-27:0.45:27);
    %mapping true solution grid to pixel grid using magnification and
    %origin shift
    Xt1=Xt./xscale+xorigin;
    Yt1=Yt./yscale+yorigin;
    % get limits of true solution grid in pixel coordinates
    ax=min(Xt1(:));
    bx=max(Xt1(:));
    ay=min(Yt1(:));
    by=max(Yt1(:));
    
    %Load 1 prana sol
    soltemp=load(fullfile(inputdir.vortex.pranasol,[inputdir.vortex.pranabaseMC, num2str(1+foffsetprana,'%05.0f'),'.mat']));
    %prana solution grid
    Xsol1=soltemp.X;
    Ysol1=soltemp.Y;
    %Load 1 Davis sol
    davissoltemp=readimx(fullfile(inputdir.vortex.davisol,[inputdir.vortex.davisbase,num2str(1+foffsetdavis,'%05.0f'),'.vc7']));
    %davis solution grid 
    D=create2DVec(davissoltemp.Frames{1,1});
    Xd1=D.X(:,:)';
    Yd1=D.Y(end:-1:1,:)';
    clear D soltemp davissoltemp;   
    % indidces corresponding to region of interest
    
    xvec=Xsol1(1,:);
    yvec=Ysol1(:,1);
    xvec1=Xd1(1,:);
    yvec1=Yd1(:,1);
    
    %Get indices which will be considered from the prana solution
    xind=find(xvec>ax & xvec<bx);
    yind=find(yvec>ay & yvec<by);
    xind1=find(xvec1>ax & xvec1<bx);
    yind1=find(yvec1>ay & yvec1<by);
    
    Xp=Xsol1(yind,xind);
    Yp=Ysol1(yind,xind);
    Xd=Xd1(yind1,xind1);
    Yd=Yd1(yind1,xind1);
    
    sxp=length(yvec);syp=length(xvec);
%     sxd=64;syd=80;
    sxp1=length(yind);syp1=length(xind);
    sxd1=length(yind1);syd1=length(xind1);
    
    Up=zeros(sxp1,syp1,N);
    Vp=zeros(sxp1,syp1,N);
    Ud=zeros(sxd1,syd1,N);
    Vd=zeros(sxd1,syd1,N);
    
    Utrue=zeros(sxp1,syp1,N);
    Vtrue=zeros(sxp1,syp1,N);
    Utrued=zeros(sxd1,syd1,N);
    Vtrued=zeros(sxd1,syd1,N);
    
    err_up=zeros(sxp1,syp1,N);
    err_vp=zeros(sxp1,syp1,N);
    err_ud=zeros(sxd1,syd1,N);
    err_vd=zeros(sxd1,syd1,N);
    
    UMCx=zeros(sxp1,syp1,N);
    UMCy=zeros(sxp1,syp1,N);
    UIMx=zeros(sxp1,syp1,N);
    UIMy=zeros(sxp1,syp1,N);
    UCSx=zeros(sxd1,syd1,N);
    UCSy=zeros(sxd1,syd1,N);
    Udiffstore=zeros(sxp1,syp1,N);
    Vdiffstore=zeros(sxp1,syp1,N);
%      keyboard;
    %% Loop through N fields
    for i=fstart:N
        
        %% Load prana solution
        pranasol=load(fullfile(inputdir.vortex.pranasol,[inputdir.vortex.pranabaseMC, num2str(i+foffsetprana,'%05.0f'),'.mat']));
        imsol=load(fullfile(inputdir.vortex.IMdir,[inputdir.vortex.pranabaseIM, num2str(i+foffsetprana,'%05.0f'),'.mat']));
%         keyboard;
        Up(:,:,i)=pranasol.U(yind,xind);
        Vp(:,:,i)=pranasol.V(yind,xind);
       
        %% Load Davis Solution
        davissol=readimx(fullfile(inputdir.vortex.davisol,[inputdir.vortex.davisbase,num2str(i+foffsetdavis,'%05.0f'),'.vc7']));
        Utemp=cell2mat(davissol.Frames{1,1}.Components{1,1}.Planes)';
        Vtemp=cell2mat(davissol.Frames{1,1}.Components{2,1}.Planes)';
        
        Utemp=Utemp(end:-1:1,:);
        Vtemp=-Vtemp(end:-1:1,:);
        Ud(:,:,i)=Utemp(yind1,xind1);
        Vd(:,:,i)=Vtemp(yind1,xind1);
       
%         keyboard;
        %% Get True Solution
        Uttemp=(spf/xscale).*truesol.Ut(:,:,i);
        Vttemp=(spf/yscale).*truesol.Vt(:,:,i);
%         keyboard;
        Utrue(:,:,i)=interp2(Xt1,Yt1,Uttemp,Xp,Yp,'linear',0);
        Vtrue(:,:,i)=interp2(Xt1,Yt1,Vttemp,Xp,Yp,'linear',0);
        Utrued(:,:,i)=interp2(Xt1,Yt1,Uttemp,Xd,Yd,'linear',0);
        Vtrued(:,:,i)=interp2(Xt1,Yt1,Vttemp,Xd,Yd,'linear',0);
        
        
        %% Calculate Error
        err_up(:,:,i)=(Up(:,:,i)-Utrue(:,:,i));
        err_vp(:,:,i)=(Vp(:,:,i)-Vtrue(:,:,i));
        err_ud(:,:,i)=(Ud(:,:,i)-Utrued(:,:,i));
        err_vd(:,:,i)=(Vd(:,:,i)-Vtrued(:,:,i));
        
        %% Get MC uncertainty estimates
        Ixx=pranasol.ixx(yind,xind);
        Iyy=pranasol.iyy(yind,xind);
        scaling=pranasol.mi(yind,xind);
        biasx1=reshape(pranasol.uncertainty(:,15),sxp,syp); 
        biasy1=reshape(pranasol.uncertainty(:,16),sxp,syp);
        Autod1=reshape(pranasol.uncertainty(:,dn),sxp,syp);
        biasx=biasx1(yind,xind);
        biasy=biasy1(yind,xind);
        Autod=Autod1(yind,xind);
        
       
        %Gradient correction
        Udiff=socdiff(Up(:,:,i),16,1);
        Vdiff=socdiff(Vp(:,:,i),16,2);
%         Udiff=gradient_compact_rich(Up(:,:,i),16,1);
%         Vdiff=gradient_compact_rich(Vp(:,:,i),16,2);
        Udiffstore(:,:,i)=Udiff;
        Vdiffstore(:,:,i)=Vdiff;
        
        Ixxt= real(sqrt(Ixx.^2 - (Autod.^2/16).*(Udiff).^2 ));
        Iyyt= real(sqrt(Iyy.^2 - (Autod.^2/16).*(Vdiff).^2 ));
        %scaling and bias correction
        UMCx(:,:,i)=sqrt(biasx.^2+(Ixxt.^2)./scaling);
        UMCy(:,:,i)=sqrt(biasy.^2+(Iyyt.^2)./scaling);
        
        
        %% Get IM uncertainty estimates
        UIMx(:,:,i)=imsol.imx(yind,xind);
        UIMy(:,:,i)=imsol.imy(yind,xind);
        
        %% Get CS uncertainty estimates
        Unrxtemp=cell2mat(davissol.Frames{1,1}.Components{22,1}.Planes)';
        Unrytemp=cell2mat(davissol.Frames{1,1}.Components{23,1}.Planes)';
        Unbxtemp=cell2mat(davissol.Frames{1,1}.Components{25,1}.Planes)';
        Unbytemp=cell2mat(davissol.Frames{1,1}.Components{26,1}.Planes)';
                
        Unrx=Unrxtemp(end:-1:1,:);
        Unry=Unrytemp(end:-1:1,:);
        Unbx=Unbxtemp(end:-1:1,:);
        Unby=Unbytemp(end:-1:1,:);
        
        UCSx(:,:,i)=sqrt(Unrx(yind1,xind1).^2+Unbx(yind1,xind1).^2);
        UCSy(:,:,i)=sqrt(Unry(yind1,xind1).^2+Unby(yind1,xind1).^2);
        
   
    end
    %% Save error and uncertainty values
    if savemat==1
        save(fullfile(outputdir{1},'Vortex_Ring.mat'),'err_up','err_vp','err_ud','err_vd','UMCx','UMCy','UIMx','UIMy','UCSx','UCSy','Up','Vp','Ud','Vd','Xp','Yp','Xt','Yt','Xd','Yd','Udiff','Vdiff');
    end
    clearvars -except casename inputdir outputdir savemat dn;
    fprintf('\nDone...\n');
    
    
end

%%
if strcmp(casename{5},'Jetdata')
    fprintf('\n Jetdata case uncertainty calculation \n');
    
     %% Load True Solution
    truesol=load(fullfile(inputdir.jet.truesoldir,[inputdir.jet.truesolbase,'.mat']));
    
    %% Set Parameters
    N=496;
    fstart=1;
    foffsetprana=2;
    foffsetdavis=0;
    
    % Previous
    %X 316 to 404 (88); Y 201 to 329 (128)
    % Now
    % X 320 to 392 (72); Y 205 to 325 (120)
    c1=5;c2=23; %303+18-1=320  303+90-1=392  X direction
    r1=5;r2=35; %188+18-1=205  188+138-1=325 Y direction
    
%     I=78;J=70;
    % Sx=66;Sy=46;
%     sxp=40;syp=30;
    
    %tsx=6:2:71;
    %tsy=5:2:50;
    tsx=8:2:69;
    tsy=7:2:43;
    
    yind=r2-r1+1;xind=c2-c1+1;   
    yind1=r2-r1+1;xind1=c2-c1+1;   
    
%     sxp=length(yvec);syp=length(xvec);
%     sxd=64;syd=80;
    sxp1=yind;syp1=xind;
    sxd1=yind1;syd1=xind1;
    
    Up=zeros(sxp1,syp1,N);
    Vp=zeros(sxp1,syp1,N);
    Ud=zeros(sxd1,syd1,N);
    Vd=zeros(sxd1,syd1,N);
    
    Utrue=zeros(sxp1,syp1,N);
    Vtrue=zeros(sxp1,syp1,N);
%     Utrued=zeros(sxd1,syd1,N);
%     Vtrued=zeros(sxd1,syd1,N);
    
    err_up=zeros(sxp1,syp1,N);
    err_vp=zeros(sxp1,syp1,N);
    err_ud=zeros(sxd1,syd1,N);
    err_vd=zeros(sxd1,syd1,N);
    
    UMCx=zeros(sxp1,syp1,N);
    UMCy=zeros(sxp1,syp1,N);
    UIMx=zeros(sxp1,syp1,N);
    UIMy=zeros(sxp1,syp1,N);
    UCSx=zeros(sxd1,syd1,N);
    UCSy=zeros(sxd1,syd1,N);
    Udiffstore=zeros(sxp1,syp1,N);
    Vdiffstore=zeros(sxp1,syp1,N);
    
    %% Loop through N fields
    for i=fstart:N
        
        %% Load prana solution
        pranasol=load(fullfile(inputdir.jet.pranasol,[inputdir.jet.pranabaseMC, num2str(i+foffsetprana,'%05.0f'),'.mat']));
        imsol=load(fullfile(inputdir.jet.IMdir,[inputdir.jet.pranabaseIM, num2str(i+foffsetprana,'%05.0f'),'.mat']));
        
        if i==1
            Xp=pranasol.X+302;
            Yp=pranasol.Y+187;
            sxp=size(Xp,1);syp=size(Yp,2);
            
            Xp=Xp(r1:r2,c1:c2);
            Yp=Yp(r1:r2,c1:c2);
        end
%         keyboard;
        Ustemp=pranasol.U(end:-1:1,:);
        Vstemp=-pranasol.V(end:-1:1,:);
        Up(:,:,i)=Ustemp(r1:r2,c1:c2);
        Vp(:,:,i)=Vstemp(r1:r2,c1:c2);
        
%        keyboard;
        %% Load Davis Solution
        davissol=readimx(fullfile(inputdir.jet.davisol,[inputdir.jet.davisbase,num2str(i+foffsetdavis,'%05.0f'),'.vc7']));
%         keyboard;
        if i==1
            [D] = create2DVec(davissol.Frames{1,1});
            Xd=D.X;
            Yd=D.Y;
            Xd=Xd';
            Yd=Yd';
            Xd=Xd(r1:r2,c1:c2)+302;
            Yd=Yd(r1:r2,c1:c2)+187;
            clear D;
        end
    
        Udtemp=cell2mat(davissol.Frames{1,1}.Components{1,1}.Planes)';
        Vdtemp=cell2mat(davissol.Frames{1,1}.Components{2,1}.Planes)';
        Ud(:,:,i)=Udtemp(r1:r2,c1:c2);
        Vd(:,:,i)=Vdtemp(r1:r2,c1:c2);
        
        
       
%         keyboard;
        %% Get True Solution
        Uttemp=truesol.Ut(:,:,i);
        Vttemp=truesol.Vt(:,:,i);
%         keyboard;
        Utrue(:,:,i)=Uttemp(tsx,tsy);
        Vtrue(:,:,i)=Vttemp(tsx,tsy);
        
        if i==fstart
            Xt=truesol.X(tsx,tsy);
            Yt=truesol.Y(tsx,tsy);
        end
        
        
        %% Calculate Error
        err_up(:,:,i)=(Up(:,:,i)-Utrue(:,:,i));
        err_vp(:,:,i)=(Vp(:,:,i)-Vtrue(:,:,i));
        err_ud(:,:,i)=(Ud(:,:,i)-Utrue(:,:,i));
        err_vd(:,:,i)=(Vd(:,:,i)-Vtrue(:,:,i));
        
        %% Get MC uncertainty estimates
        Ixxtemp=pranasol.ixx(end:-1:1,:);
        Iyytemp=pranasol.iyy(end:-1:1,:);
        scalingtemp=pranasol.mi(end:-1:1,:);
%         keyboard;
        biasx1=reshape(pranasol.uncertainty(:,15),sxp,syp); 
        biasy1=reshape(pranasol.uncertainty(:,16),sxp,syp);
        Autod1=reshape(pranasol.uncertainty(:,dn),sxp,syp);
        biasx1=biasx1(end:-1:1,:);
        biasy1=biasy1(end:-1:1,:);
        Autod1=Autod1(end:-1:1,:);
        
        Ixx=Ixxtemp(r1:r2,c1:c2);
        Iyy=Iyytemp(r1:r2,c1:c2);
        scaling=scalingtemp(r1:r2,c1:c2);
        biasx=biasx1(r1:r2,c1:c2);
        biasy=biasy1(r1:r2,c1:c2);
        Autod=Autod1(r1:r2,c1:c2);
        
       
        %Gradient correction
        Udiff=socdiff(Up(:,:,i),16,1);
        Vdiff=socdiff(Vp(:,:,i),16,2);
%         Udiff=gradient_compact_rich(Up(:,:,i),16,1);
%         Vdiff=gradient_compact_rich(Vp(:,:,i),16,2);
        Udiffstore(:,:,i)=Udiff;
        Vdiffstore(:,:,i)=Vdiff;
        
        Ixxt= real(sqrt(Ixx.^2 - (Autod.^2/16).*(Udiff).^2 ));
        Iyyt= real(sqrt(Iyy.^2 - (Autod.^2/16).*(Vdiff).^2 ));
        %scaling and bias correction
        UMCx(:,:,i)=sqrt(biasx.^2+(Ixxt.^2)./scaling);
        UMCy(:,:,i)=sqrt(biasy.^2+(Iyyt.^2)./scaling);
        
        
        %% Get IM uncertainty estimates
        imxtemp=imsol.imx(end:-1:1,:);
        imytemp=imsol.imy(end:-1:1,:);
        UIMx(:,:,i)=imxtemp(r1:r2,c1:c2);
        UIMy(:,:,i)=imytemp(r1:r2,c1:c2);
        
        %% Get CS uncertainty estimates
        Unrxtemp=cell2mat(davissol.Frames{1,1}.Components{22,1}.Planes)';
        Unrytemp=cell2mat(davissol.Frames{1,1}.Components{23,1}.Planes)';
        Unbxtemp=cell2mat(davissol.Frames{1,1}.Components{25,1}.Planes)';
        Unbytemp=cell2mat(davissol.Frames{1,1}.Components{26,1}.Planes)';
                
        Unrx=Unrxtemp(r1:r2,c1:c2);
        Unry=Unrytemp(r1:r2,c1:c2);
        Unbx=Unbxtemp(r1:r2,c1:c2);
        Unby=Unbytemp(r1:r2,c1:c2);
        
        UCSx(:,:,i)=sqrt(Unrx.^2+Unbx.^2);
        UCSy(:,:,i)=sqrt(Unry.^2+Unby.^2);
%         keyboard;
   
    end
    %% Save error and uncertainty values
    if savemat==1
        save(fullfile(outputdir{1},'Jetdata.mat'),'err_up','err_vp','err_ud','err_vd','UMCx','UMCy','UIMx','UIMy','UCSx','UCSy','Up','Vp','Ud','Vd','Xp','Yp','Xt','Yt','Xd','Yd','Udiff','Vdiff');
    end
    clearvars -except casename inputdir outputdir savemat dn;
    fprintf('\nDone...\n');
    
    
end


end