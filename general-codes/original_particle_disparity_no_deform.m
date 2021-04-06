function [deltax,deltay,Npartu,Npartv]=original_particle_disparity_no_deform(im1,im2,X,Y,U,V,ws)
% im1 and im2 are original images
% X,Y are vector grid points
%U,V are U and V velocities computed in that pass
%ws window size in prana

rsearch=1;  %search radius
nRmsLength=1; % kernel of the low pass filter (particle image diameter)
weight='sqrt'; %'peaks', 'I1I2', 'sqrt'
position_weight='gauss'; %'gauss', 'tophat'

ws=ws+1;

% Region Of Interest (ROI), where the uncertainty is computed
ROI = [0 10000 0 10000]; % [xmin xmax ymin ymax]

grdx = X(1,:);
grdy = Y(:,1);


[J,I]=size(im1);
im=zeros(J,I);
im(:,:,1)=im1; clear im1;
im(:,:,2)=im2; clear im2;

%% sliding avg subtraction
% fprintf('Sliding average subtraction: ')
im1Avg=slidingAvgDavis(im(:,:,1),ws);
im(:,:,1)=im(:,:,1)-im1Avg;
clear im1Avg
im2Avg=slidingAvgDavis(im(:,:,2),ws);
im(:,:,2)=im(:,:,2)-im2Avg;
clear im2Avg

% displacement in pixels
displ_px=cat(3,U,V);

% %% image deformation (of the image pair) through the sinc interpolation
% imdef=image_def(im,displ_px,X,Y);
% %% Gaussian Smoothing
% imdef1=squeeze(imdef(:,:,1));
% imdef2=squeeze(imdef(:,:,2));
% clear imdef
% imdef1=slidingAvgDavis(imdef1,nRmsLength);
% imdef2=slidingAvgDavis(imdef2,nRmsLength);
%% disparity matrices computation
% fprintf('Disparity computation: ')
%     keyboard;
% t5=toc;
% [xdispl,ydispl,~,peaks]=disparity_computation(imdef1,imdef2,rsearch,weight);
% [xdispl,ydispl,~,peaks]=disparity_computation(im1,im2,rsearch,weight);
[xdispl,ydispl,~,peaks] = disparity_computation(im(:,:,1), im(:,:,2), rsearch, weight);

% t6=toc;
% fprintf([num2str(t6-t5,'%4.1f') ' s.\n'])

%% Uncertainty computation
% fprintf('Uncertainty computation: ')
% [extot,exbias,exrms,Npartu,Cu]=error_window4(xdispl,peaks,grdx,grdy,ws,nRmsLength,position_weight,ROI);
% [eytot,eybias,eyrms,Npartv,Cv]=error_window4(ydispl,peaks,grdx,grdy,ws,nRmsLength,position_weight,ROI);

[extot,~,~,Npartu,~]=error_window4(xdispl,peaks,grdx,grdy,ws,nRmsLength,position_weight,ROI);
[eytot,~,~,Npartv,~]=error_window4(ydispl,peaks,grdx,grdy,ws,nRmsLength,position_weight,ROI);

deltax=extot;
deltay=eytot;


end

function [Xdisparity,Ydisparity,immult,peaks]=disparity_computation(im1,im2,rsearch,weight)
if nargin<4
    weight='peaks';
end

% new
% im1dum=slidingAvgDavis(im1,3);
% im2dum=slidingAvgDavis(im2,3);
% immult=im1dum.*im2dum;
% clear im1dum im2dum;
%

immult=im1.*im2;
[J,I]=size(im1);
im1pos=im1>0; im2pos=im2>0;
im1=im1+abs(min(im1(:)))+eps;
im2=im2+abs(min(im2(:)))+eps;

[X1,Y1]=meshgrid(1:I,1:J);
[X2,Y2]=meshgrid(1:I,1:J);
Xpos1=X1; Ypos1=Y1;
Xpos2=X2; Ypos2=Y2;
%% get the maxima
peaks=get_max(immult,weight);
peaksIndex=find(peaks>0);

%% subpixel position
[X1sub,Y1sub]=  subpixel_position(im1); X1sub(isnan(X1sub))=0; Y1sub(isnan(Y1sub))=0;
[X2sub,Y2sub]=  subpixel_position(im2); X2sub(isnan(X2sub))=0; Y2sub(isnan(Y2sub))=0;

%% main loop
for nIndex=1:numel(peaksIndex);
    peakIndexLoc=peaksIndex(nIndex);

    % discrete particles positions
    [ip1,jp1]=find_particles2(im1,peakIndexLoc,rsearch);
    [ip2,jp2]=find_particles2(im2,peakIndexLoc,rsearch);
    % subpixel position
    if (ip1*ip2*jp1*jp2>0)
        Xpos1(peakIndexLoc)=X1(jp1,ip1)+X1sub(jp1,ip1);
        Ypos1(peakIndexLoc)=Y1(jp1,ip1)-Y1sub(jp1,ip1);
        Xpos2(peakIndexLoc)=X2(jp2,ip2)+X2sub(jp2,ip2);
        Ypos2(peakIndexLoc)=Y2(jp2,ip2)-Y2sub(jp2,ip2);
    else
        peaks(peakIndexLoc)=0;
    end
end

% %% main loop with 2D Gaussian!!!!
% for nIndex=1:numel(peaksIndex);
%     peakIndexLoc=peaksIndex(nIndex);
% 
%     % discrete particles positions
%     [ip1 jp1]=find_particles2(im1,peakIndexLoc,rsearch);
%     [ip2 jp2]=find_particles2(im2,peakIndexLoc,rsearch);
%     % subpixel position with 2D Gaussian!!!!
%     if (ip1*ip2*jp1*jp2>0 & jp1>2 & ip1>2 & jp2>2 & ip2>2 & jp1<(J-1) & jp2<(J-1) & ip1<(I-1) & ip2<(I-1))
%         im1loc=im1(jp1-1:jp1+1,ip1-1:ip1+1);
%         im2loc=im2(jp2-1:jp2+1,ip2-1:ip2+1);
%         [dx1,dy1]=GaussFit2D(im1loc,5,'TH');
%         [dx2,dy2]=GaussFit2D(im2loc,5,'TH');
%         Xpos1(peakIndexLoc)=X1(jp1,ip1)+dx1;
%         Ypos1(peakIndexLoc)=Y1(jp1,ip1)-dy1;
%         Xpos2(peakIndexLoc)=X2(jp2,ip2)+dx2;
%         Ypos2(peakIndexLoc)=Y2(jp2,ip2)-dy2;
% %         figure(1), clf, imagesc(im1loc), axis equal, axis tight
% %         dx1,dy1, pause
%     else
%         peaks(peakIndexLoc)=0;
%     end
% end

%% Disparity matrices
Xdisparity=Xpos2-Xpos1; 
Ydisparity=Ypos2-Ypos1;

%% Thresholding
Thr= 0*std(immult(:))/3;
Xdisparity=Xdisparity.*(immult>Thr).*im1pos.*im2pos;
Ydisparity=Ydisparity.*(immult>Thr).*im1pos.*im2pos;
peaks=peaks.*(immult>Thr).*im1pos.*im2pos;
end

function [xsub,ysub]=  subpixel_position(im)
% function that computes the subpixel position of the particles

% to avoid non positive values
if min(im(:))<0
    im=im+abs(min(im(:)))+eps;
end

%logarithm of the image
im=log(im);

%shifted images
[J,I]=size(im);
imE=zeros(J,I)+eps; imW=zeros(J,I)+eps; imN=zeros(J,I)+eps; imS=zeros(J,I)+eps;
imE(:,1:end-1)=im(:,2:end);
imW(:,2:end)=im(:,1:end-1);
imN(2:end,:)=im(1:end-1,:);
imS(1:end-1,:)=im(2:end,:);

%subpixel position
xsub=(imW-imE)./2./(imE+imW-2*im);
ysub=(imS-imN)./2./(imS+imN-2*im);
end
function peaks = get_max(im,weight)
% function that get the peaks positions of the image im
% im: 2D matrix of size (J,I)
% peaks: 2D matrix of size(J,I) composed by 1 in the peaks positions and 0 otherwise

[J,I]=size(im);

imC=im(2:J-1,2:I-1); % image in the center (excludind a border of 1 pixel)
imW=im(2:J-1,1:I-2);
imE=im(2:J-1,3:I);
imN=im(1:J-2,2:I-1);
imS=im(3:J,2:I-1);

imNW=im(1:J-2,1:I-2); imNE=im(1:J-2,3:I);
imSW=im(3:J,1:I-2); imSE=im(3:J,3:I);

peaks=zeros(J,I);
peaks(2:J-1,2:I-1)=(imC>imW & imC>imE & imC>imN & imC>imS & imC>imNW & imC>imNE & imC>imSW & imC>imSE);

if strcmp(weight,'I1I2')
    peaks=peaks.*im;
elseif strcmp(weight,'sqrt')
    peaks=peaks.*sqrt(abs(im));
elseif strcmp(weight,'peaks')~=1
    error('weight must be peaks, sqrt or I1I2')
end %if strcmp

%figure,imagesc(peaks), axis equal, axis tight
% pause
end % end function

function [iPeak,jPeak]=find_particles2(im,index,rsearch)

% detect the index of the discrete particle position with a loop on the search radius
[J,I]=size(im);
[jp,ip]=ind2sub(size(im),index);
r=0; iPeak=0; jPeak=0;
imax=0;
while (r<=rsearch && imax==0);
    for j=max(2,jp-r):min(J-1,jp+r)
        for i=max(2,ip-r):min(I-1,ip+r)
            if (im(j,i)>im(j-1,i) && im(j,i)>im(j,i-1) && im(j,i)>im(j+1,i) && im(j,i)>im(j,i+1) && im(j,i)>imax)
                jPeak=j; iPeak=i; imax=im(j,i);
            end
        end
    end
    r=r+1;
end

end % function

function [etot,ebias,erms,Npart,C] = error_window4(displ,peaks,grdx,grdy,ws,~,position_weight,ROI)
% function that computes the "bias" and rms error averaging the errors in an interrogation window
% peaks: 2D matrix of size (J,I): it is zero where no peaks are present, in
% the other points it expresses the peak displacement
% grdx: vector of the grid points in the x direction
% grdy: vector of the grid points in the y direction
% ws: window size

%% initialization
lenx=length(grdx); leny=length(grdy); [J,I]=size(peaks);
% if grdy(end)<grdy(1) grdy=flipud(grdy); end %to have both grdx and grdy in increasing order
r=(ws-1)/2; %radius of the window
if strcmp(position_weight,'tophat')
    coeff=1;
    gaussweight=ones(2*r+1);
elseif strcmp(position_weight,'gauss')
    coeff=1.75;%coeff=1.4;
    gaussweight=gaussian2De_2(round(2*coeff*r+1),round(2*coeff*r+1));
else
    error('"position_weight" has to be "tophat" or "gauss"');
end

ebias=zeros(leny,lenx); erms=zeros(leny,lenx);
Npart=zeros(leny,lenx);
etot=Npart;
% C=cell(leny,lenx);
C=0;
% keyboard;
%% main loop
for jgrid=1:leny
    j=grdy(jgrid);
    jmin=max(1,j-coeff*r); jmax=min(J,j+coeff*r);
    for igrid=1:lenx
        i=grdx(igrid);
        
        if (i >= ROI(1) && i <= ROI(2) && j>=ROI(3) && j<=ROI(4))
            imin=max(1,i-coeff*r); imax=min(I,i+coeff*r);
            
            displwin=displ(jmin:jmax,imin:imax); %displacement in the window
            gaussloc=gaussweight(1+coeff*r-(j-jmin):1+coeff*r+jmax-j,1+coeff*r-(i-imin):1+coeff*r+(imax-i));
            peakswin=peaks(jmin:jmax,imin:imax); %peaks in the window
            peakswin=peakswin.*gaussloc;
            peaksold=peakswin;
            displold=displwin;
            displwin=displold(displold~=0);
            peakswin=peakswin(displold~=0);
            
            % outliers removal
            outliers=OutliersPercentile(displwin);
            displwin(outliers==1)=[];
            peakswin(outliers==1)=[];
            
            % computation
            ebiasloc=sum(displwin.*peakswin)./sum(peakswin);
            ermsloc=sqrt(sum(peakswin.*(displwin-ebiasloc).^2)/sum(peakswin));
            ebias(jgrid,igrid)=ebiasloc;
            Nweight=(peaksold>0).*gaussloc; Nweight=sum(Nweight(:));
            if Nweight<1; Nweight=1; end
            erms(jgrid,igrid)=ermsloc;
            etot(jgrid,igrid)=sqrt(ermsloc.^2/Nweight+ebiasloc.^2);
            Npart(jgrid,igrid)=Nweight;%N
        end % if ROI
    end
end
% keyboard;

ebias(isnan(ebias))=0;
erms(isnan(erms))=0;
etot(isnan(etot))=0;
end % end function

function outliers = OutliersPercentile(vec)
% outliers is 1 in case of outlier, 0 otherwise

outliers=abs(vec-mean(vec))>3*std(vec);
% outliers=1-(vec>prctile(vec,10) & vec<prctile(vec,90));
end

function [g]=gaussian2De_2(Ny,Nx)
if nargin<2, Nx=Ny; end
% 2D gaussian function
[X,Y]=meshgrid(-(Nx-1)/2:(Nx-1)/2,-(Ny-1)/2:(Ny-1)/2);
r=sqrt(X.^2+Y.^2);
r=r/(Nx-1)*2;
g=exp(-2*r.^2);
end %end function

function matrixOut = slidingAvgDavis(matrixIn, ker)

%Initial error statements and definitions
if nargin<2, error('Not enough input arguments'), end
if length(ker(1))~=1, error('Nr must be a scalar'), end
% keyboard;
% smoothing
if ker>1
    [J,I]=size(matrixIn);
    
    % rim construction
    if mod(I,2)==0; Ir=I+1; else Ir=I; end
    if mod(J,2)==0; Jr=J+1; else Jr=J; end
    matrixIn2=zeros(Jr,Ir);
    matrixIn2(1:J,1:I)=matrixIn;
    %[Xr Yr]=meshgrid(1:Ir,1:Jr);
    %z=find(matrixIn2==0); nz=find(matrixIn2~=0);
    %matrixIn2(z)=griddata(Xr(nz),Yr(nz),matrixIn2(nz),Xr(z),Yr(z),'nearest');
    
    %filt=zeros(size(matrixIn2));
    [X,Y]=meshgrid(-(Ir-1)/2:(Ir-1)/2,-(Jr-1)/2:(Jr-1)/2);
    sigma=-log(1-2/ker);
    r=abs(X)+abs(Y);
    filt=exp(-sigma*r);
    filt(J+1:Jr,I+1:Ir)=0;
    filt=filt/sum(filt(:));
    Ishift=(Ir+1)/2; Jshift=(Jr+1)/2;
    
    Mf=ifft2(fft2(filt).*fft2(matrixIn2));
    matrixOut=circshift(Mf,[Jshift,Ishift]);
    
    matrixOut=matrixOut(1:J,1:I);
else
    matrixOut=matrixIn;
end % if ker>1

end % end function

function imdef = image_def(im,displ,X,Y)
% function that performs the image deformation through the sinc
% interpolation
% im: image pair, having size (J,I,2)
% displ: displacement in the two directions, of size (leny,lenx,2)
% X and Y: vector positions, of size (J,I)

%% dense predictor (displacement in each pixel position)
displ(:,:,2)=-displ(:,:,2);
[J,I,~]=size(im);
densepred=zeros(J,I,2);
grdx=X(1,:); grdy=Y(:,1);
for h=1:2
    densepred(:,:,h)=build_predictor(grdx,grdy,displ(:,:,h),I,J);	% dense predictor 
end
[tmpx,tmpy]=meshgrid(1:I,1:J);

%% displacement to be imposed to the two images
xtransf=cell(1,2); ytransf=cell(1,2);
xtransf{1}=tmpx - 0.5*densepred(:,:,1);
xtransf{2}=tmpx + 0.5*densepred(:,:,1);
ytransf{1}=tmpy - 0.5*densepred(:,:,2);
ytransf{2}=tmpy + 0.5*densepred(:,:,2);
clear tmpx tmpy densepred
% keyboard;
%% deformation (sinc interpolation for the images, based on the predictor)
imdef=zeros(J,I,2);
for h=1:2
    %imdef(:,:,h)=sincinterp_even(im(:,:,h),xtransf{h},ytransf{h},6);
    %imdef(:,:,h)=sincinterp(im(:,:,h),xtransf{h},ytransf{h},3);
%     imdef(:,:,h)=Whittaker7(im(:,:,h),xtransf{h},ytransf{h},3);
%     imdef(:,:,h)=sincBlackmanInterp2(im(:,:,h),xtransf{h},ytransf{h},3);
imdef(:, :, h) = whittaker_blackman(im(:, :, h), xtransf{h}, ytransf{h}, 6, 1);
end
% [Xim Yim]=meshgrid(1:I,1:J);
% parfor h=1:2
%      imdef(:,:,h)=interp2(Xim,Yim,im(:,:,h),xtransf{h},ytransf{h},'*spline');
% end
% imdef(isnan(imdef))=0;

end
function densepred=build_predictor(grdxold,grdyold,displ,dimx,dimy) %iter,densepredold

% initialize densepred
flip='no';
if grdyold(end)<grdyold(1)
    flip='yes';
    grdyold=flipud(grdyold);
    displ=flipud(displ);
end
densepred=zeros(grdyold(end),grdxold(end));

% if iter==1
%     [X,Y]=meshgrid(grdxold,grdyold);
%     [XI,YI]=meshgrid(1:dimx,1:dimy);
%     
%     densepred=(interp2(X,Y,displ,XI,YI,'*linear'));
%     densepred(find(isnan(densepred)))=0;
% else
[X,Y]=meshgrid(grdxold,grdyold);

xin=grdxold(1);
xfin=grdxold(end);
yin=grdyold(1);
yfin=grdyold(end);

%If percentage of D is < 50% or iteration number is < mask iteration value then make denspred for whole image, else apply mask setting old values in D.
[Yline,Xline]=find(ones(dimy,dimx));
% keyboard;
%densepred=densepredold;
densepredline=(interp2(X,Y,displ,Xline,Yline,'*linear'));
bufdensepred=densepred(:);
bufdensepred(Yline+(dimy*(Xline-1)))=densepredline;
densepred=reshape(bufdensepred,dimy,dimx);

%(yin:yfin,xin:xfin)
%Left boundary
tmp2=densepred(yin,:);
densepred(1:yin-1,:)=tmp2(  ones( length(1:yin-1),1 ),:  );
%Right boundary
tmp2=densepred(yfin,:);
densepred(yfin+1:dimy,:)=tmp2(  ones( length(yfin+1:dimy),1 ) ,:  );
%Upper boundary
tmp2=densepred(:,xin)';
densepred(:,1:xin-1)=tmp2(  ones( length(1:xin-1),1 ) ,:  )';
%Lower boundary
tmp2=densepred(:,xfin)';
densepred(:,xfin+1:dimx)=tmp2(  ones( length(xfin+1:dimx),1 ) ,:  )';
%Take out NaN's
densepred(isnan(densepred))=0;

% end
if strcmp(flip,'yes')
    densepred=flipud(densepred);
end

%contourf(XI,YI,densepred), colorbar
end


