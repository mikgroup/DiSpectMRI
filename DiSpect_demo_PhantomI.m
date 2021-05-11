%% This is a demo processing of DiSpect for phantom I:
load PhantomI;
%% Displacement encoded echo seperation from 4-step phase cycling
PhaseCyclingMat4=[-1, -1i, +1, 1i; ...
                  -1, 1i, +1, -1i; ...
                  +1, +1, +1, +1];             
APDispEncoded_ksp=PhantomI.APDispEnc_ksp;
LRDispEncoded_ksp=PhantomI.LRDispEnc_ksp;
[Nx,Ny,Nd_AP,Nc,Np]=size(APDispEncoded_ksp);
Nd_LR=size(LRDispEncoded_ksp,3);
APDispDecoded_ksp=reshape(APDispEncoded_ksp,[],Np)*PhaseCyclingMat4.';
LRDispDecoded_ksp=reshape(LRDispEncoded_ksp,[],Np)*PhaseCyclingMat4.';
% Decoded echoes:#1->Displacement encoded echo; 
%                #2->T1 echo; 
%                #3->Conjugated Displacement encoded echo;
APDispDecoded_ksp=reshape(APDispDecoded_ksp(:,1),[Nx,Ny,Nd_AP,Nc]);
LRDispDecoded_ksp=reshape(LRDispDecoded_ksp(:,1),[Nx,Ny,Nd_LR,Nc]);
%% FT operation for both datasets 
imD_AP=cifftn(APDispDecoded_ksp,[1 2 3]);
imD_LR=cifftn(LRDispDecoded_ksp,[1 2 3]);
%% Shearing transformation to obtain the absolute physical coordinate from displacement spectra 
[im_AP,FOVd_AP]=DiSpect_shearing(imD_AP,[2 3],[PhantomI.FOVinf(2),PhantomI.FOVinf(4)]);
[im_LR,FOVd_LR]=DiSpect_shearing(imD_LR,[1 3],[PhantomI.FOVinf(1),PhantomI.FOVinf(3)]);
%% PhantomI has a very simple structure, we only used the two scans to create the source images
im_trace=bsxfun(@times,reshape(im_LR,[Nx,Ny,size(im_LR,3),1,Nc]),reshape(im_AP,[Nx,Ny,1,size(im_AP,3),Nc]));
%% Prepare the dataset for visulization
% The first two dimensions are target regions while the third-fourth
% dimensions are the sources
% However, the FOV information is different for the targets and sources
% Let's do some interplation (if pixel resolution is different), windowing, cropping for better visualization of the dataset
Nd_LR=size(im_trace,3);
Nd_AP=size(im_trace,4);
im_trace_show=crop(sos(im_trace),[Nx,Ny,Nx,Ny]);
im_trace_show=cifftn(bsxfun(@times,cfftn(im_trace_show,[3 4]),reshape((hamming(Nx)*hamming(Ny).'),[1 1 Nx,Ny])),[3 4]);
NNx=4*Nx;
NNy=4*Ny;
im_trace_show=imresize(im_trace_show,[NNx,NNy]);
imaging_target_all=sos(sos(im_trace_show,3),4);

Target_positions=[216 208 196 184 172 160 148 132 120 108 96 84 72 60 48 36 216 208 196 184 172 160 148 132 120 108  96  84  72 60  48 36;
                  82  48  39  39  39  24  18  18  18  33  39 39 39 64 82 82 82  116 125 125 125 140 146 146 146 131 125 125 125 100 82 82];             
Npos=size( Target_positions,2);  
BWbase=zeros(NNx,NNy);
[xx, yy]=ndgrid(1:NNx,1:NNy);
BWbase(((xx-NNx/2).^2+(yy-NNy/2).^2)<36)=1;      
im_target_overlay=zeros(NNx,NNy,Npos);
im_source_overlay=zeros(NNx,NNy,Npos);
for m=1:size( Target_positions,2)
    BWcurrent=circshift(BWbase,[Target_positions(1,m)-NNx/2,Target_positions(2,m)-NNy/2]);
    im_target_overlay(:,:,m)=imaging_target_all.*BWcurrent;
    im_source_overlay(:,:,m)=imresize(squeeze(sos(sos(bsxfun(@times,im_trace_show,BWcurrent),1),2)),[NNx,NNy]);
end
%% Generate a DiSpect tracing viewer of the target images with corresponding source images

winMnWidth=800/2;
winMnHeight=600/2;
ScreenSize=get(0,'ScreenSize');

handles.winMain=figure('Name','Diplacement MRI tracing',...
		'Position',[ScreenSize(3)/2-winMnWidth/2 ScreenSize(4)/2-winMnHeight/2 winMnWidth winMnHeight],...
		'NumberTitle','off',...
		'Tag','zzMain',...
        'Units','normalized',...
		'MenuBar','figure');

handles.imAxes=axes('Position',[0.01 0 0.63 0.80],...
		'Layer','top',...
        'HandleVisibility','on',...
		'Box','ON',...
        'color',[0 0 0],...
        'Units','normalized',...
		'Tag','imAxes',...
		'HitTest','on'); 
  
handles.imshow=imshow(imaging_target_all.',[],'InitialMagnification','fit','Parent',handles.imAxes) ;

pimshow = plotboxpos(handles.imAxes);
handles.xSpecAxes=axes('Position',[pimshow(1) pimshow(2)+pimshow(4) pimshow(3) pimshow(4)/2],...
		'Layer','top',...
        'HandleVisibility','on',...
		'Box','off',...
        'color',[1 1 1],...
        'Units','normalized',...
		'Tag','xSpecAxes',...
		'HitTest','on'); 
    
handles.ySpecAxes=axes('Position',[pimshow(1)+pimshow(3) pimshow(2) pimshow(3)/2 pimshow(4)],...
		'Layer','top',...
        'HandleVisibility','on',...
		'Box','off',...
        'color',[1 1 1],...
        'Units','normalized',...
		'Tag','ySpecAxes',...
		'HitTest','on'); 
handles.imAxes2=axes('Position',[pimshow(1)+pimshow(3) pimshow(2)+pimshow(4) pimshow(3)/2 pimshow(4)/2],...
		'Layer','top',...
        'HandleVisibility','on',...
		'Box','off',...
        'color',[1 1 1],...
        'Units','normalized',...
		'Tag','imAxes2',...
		'HitTest','on'); 

      
im_phantom= imread('phantom3Dprint.png');
handles.imshow=imagesc(im_phantom,'Parent',handles.imAxes2) ;axis off;

set(handles.imAxes2,'xtick',[])
set(handles.imAxes2,'ytick',[])
set(handles.xSpecAxes,'xtick',[])
set(handles.xSpecAxes,'ytick',[])
set(handles.ySpecAxes,'xtick',[])
set(handles.ySpecAxes,'ytick',[])

%%
cgreen=zeros(NNy,NNx,3);
cgreen(:,:,2)=1;
cred=zeros(NNy,NNx,3);
cred(:,:,1)=1;

filename = 'PhantomI_tracing_demo.gif';
set(handles.winMain,'color',[1 1 1])
for m=1:Npos   
    BWtarget=im_target_overlay(:,:,m).';
    BWtarget=BWtarget./max(BWtarget(:));
    BWsource=(im_source_overlay(:,:,m).').^2;
    BWsource=BWsource./max(BWsource(:));
 
    if isfield(handles,'htarget')
        delete(handles.htarget);
    end
    if isfield(handles,'hsource')
        delete(handles.hsource);
    end
    axes(handles.imAxes);
    hold on;handles.htarget = imshow(cgreen);hold off;
    set(handles.htarget , 'AlphaData',BWtarget);
    hold on;handles.hsource = imshow(cred);hold off;
    set(handles.hsource , 'AlphaData', BWsource);

    handles.xspectrum=plot(sos(BWsource,1),'Parent',handles.xSpecAxes,'LineWidth',2);
    xlim(handles.xSpecAxes,[1,NNx]);
    handles.yspectrum=plot(sos(BWsource,2),'Parent',handles.ySpecAxes,'LineWidth',2);
    xlim(handles.ySpecAxes,[1,NNy]);view(handles.ySpecAxes,[90 90]);
    set(handles.xSpecAxes,'xtick',[])
    set(handles.xSpecAxes,'ytick',[])
    set(handles.ySpecAxes,'xtick',[])
    set(handles.ySpecAxes,'ytick',[])
    
    frame=getframe(handles.winMain);
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
    if m == 1 
       imwrite(imind,cm,filename,'gif','DelayTime',1/8, 'Loopcount',inf); 
    else 
       imwrite(imind,cm,filename,'gif','DelayTime',1/8,'WriteMode','append'); 
    end    
    pause(0.05)
end
%% Create a color filters for color coded displacement maps
Nd=192;
gfilter=zpad(tukeywin(Nd/2,0.4),[Nd,1]);
bfilter=crop(circshift(zpad(gfilter,[2*Nd,1]),[-round(Nd*1.05/3) 0]),[Nd,1]);
rfilter=crop(circshift(zpad(gfilter,[2*Nd,1]),[+round(Nd*1.05/3) 0]),[Nd,1]);
color_filter=cat(2,rfilter,gfilter,bfilter);
figure();plot(1:Nd,rfilter,'r',1:Nd,gfilter,'g',1:Nd,bfilter,'b')

Ndx=40;
imD_LR_adj=ipermute(interp1(permute(crop(sos(imD_LR),[Nx,Ny,Ndx]),[3 1 2]),linspace(1,Ndx,Nd)),[3 1 2]);
im_LR_RGB=squeeze(sos(bsxfun(@times,repmat(flip(imD_LR_adj,3),[1 1 1 3]),reshape(color_filter,[1 1 size(imD_LR_adj,3) 3])),3));
im_LR_RGB=im_LR_RGB./max(im_LR_RGB(:));

Ndy=24;
imD_AP_adj=ipermute(interp1(permute(crop(sos(imD_AP),[Nx,Ny,Ndy]),[3 1 2]),linspace(1,Ndy,Nd)),[3 1 2]);
im_AP_RGB=squeeze(sos(bsxfun(@times,repmat(flip(imD_AP_adj,3),[1 1 1 3]),reshape(color_filter,[1 1 size(imD_AP_adj,3) 3])),3));
im_AP_RGB=im_AP_RGB./max(im_AP_RGB(:));

figure();
subplot(1,2,1);imshow(permute(im_LR_RGB,[2 1 3]));
subplot(1,2,2);imshow(permute(im_AP_RGB,[2 1 3]));
