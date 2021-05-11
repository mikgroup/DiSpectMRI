%% This is a demo processing of DiSpect for invivo data
% Please note that echo separation was previously performed on this dataset
% The rawdata set of the displacement encoding echo: [Nro,Npe,Nd_LR,Nd_SI,N_mixing, Nc];
% The FOV information: [FOVro, FOVpe, FOVd_LR,FOVd_SI]; 
load Braindata;
%% FT operation 
% Windowing the k space to reduce the gibbs ringing
% Remove the displacement spectrum dc component (corresponding mainly to static tissues)
kspDiSpect=Braindata.kspDiSpect;
[Nro,Npe,Nd_LR,Nd_SI,N_mixing, Nc]=size(kspDiSpect);
kwin1=1;
kwin2=(hamming(Nd_LR)*hamming(Nd_SI).').^1/4;
kwin=bsxfun(@times,kwin1,permute(kwin2,[3 4 1 2]));
im_DiSpect=cifftn(bsxfun(@times,kspDiSpect,kwin),[1 2 3 4]);
im_DiSpect(:,:,end/2+1,:,:)=0;
im_DiSpect(:,:,:,end/2+1,:,:)=0;
%% Shearing transformation to obtain the absolute physical coordinate from displacement spectra
% This shearing transformation is done along the first and third dimensions
% Note that we did a circshift to make z=0 at the top of the matrix
im_DiSpect_zshift=circshift(im_DiSpect,[0 0 0 Nd_SI/2 0]);
im_trace=sos(DiSpect_shearing(im_DiSpect_zshift,[1 3],[Braindata.FOVinf(1),Braindata.FOVinf(3)]));
im_trace_adj=permute(imresize(crop(permute(im_trace,[3 4 1 2 5]),[30 16 64 64 8]),[30*4 20*4]),[3 4 1 2 5]);
im_trace_adj=flip(permute(im_trace_adj,[2 1 4 3 5]),1);
%% Let's look at the source images at different mixing times by look at serveral ROIs (Fig. 7)
BW_target_4roi=Braindata.BW_target_4roi;
im_source_roi1=squeeze(sos(sos(bsxfun(@times,im_trace_adj(:,:,:,:,2:5),BW_target_4roi(:,:,1)),1),2));
im_source_roi2=squeeze(sos(sos(bsxfun(@times,im_trace_adj(:,:,:,:,2:5),BW_target_4roi(:,:,2)),1),2));
im_source_roi3=squeeze(sos(sos(bsxfun(@times,im_trace_adj(:,:,:,:,2:5),BW_target_4roi(:,:,3)),1),2));
im_source_roi4=squeeze(sos(sos(bsxfun(@times,im_trace_adj(:,:,:,:,2:5),BW_target_4roi(:,:,4)),1),2));
titles={'ROI#1@Tm=200ms', 'ROI#1@Tm=300ms','ROI#1@Tm=400ms','ROI#1@Tm=500ms','ROI#1Sum(Tm)',...
        'ROI#2@Tm=200ms', 'ROI#2@Tm=300ms','ROI#2@Tm=400ms','ROI#2@Tm=500ms','ROI#2Sum(Tm)',...
        'ROI#3@Tm=200ms', 'ROI#3@Tm=300ms','ROI#3@Tm=400ms','ROI#3@Tm=500ms','ROI#3Sum(Tm)',...
        'ROI#4@Tm=200ms', 'ROI#4@Tm=300ms','ROI#4@Tm=400ms','ROI#4@Tm=500ms','ROI#4Sum(Tm)'};
im_data=cat(4,cat(3,im_source_roi1,sos(im_source_roi1)),...
                cat(3,im_source_roi2,sos(im_source_roi2)),...
                    cat(3,im_source_roi3,sos(im_source_roi3)),...
                        cat(3,im_source_roi4,sos(im_source_roi4)));
im_data=im_data./max(im_data(:));
figure();imshowMRI(im_data.^2,0,[4 5],titles)
%% Referring to MRA, Let's look the source images from four vascular territories by integrating across the mixing times (Figure 8)
BW_4Territories=Braindata.BW_4Territories;
im_trace_adj_thd=im_trace_adj;
im_trace_adj_thd(im_trace_adj<0.03*max(im_trace_adj(:)))=0;
im_target_background=sos(sos(sos(im_trace_adj_thd)));
im_target_background=imresize(im_target_background./max(im_target_background(:)),5,'lanczos2');

im_source_roi1=abs(imresize(squeeze(sos(sos(sos(bsxfun(@times,abs(im_trace_adj_thd(:,:,:,:,2:6)),abs(BW_4Territories(:,:,1))),1),2),5)),5,'lanczos2'));
im_source_roi2=abs(imresize(squeeze(sos(sos(sos(bsxfun(@times,abs(im_trace_adj_thd(:,:,:,:,2:6)),abs(BW_4Territories(:,:,2))),1),2),5)),5,'lanczos2'));
im_source_roi3=abs(imresize(squeeze(sos(sos(sos(bsxfun(@times,abs(im_trace_adj_thd(:,:,:,:,2:6)),abs(BW_4Territories(:,:,3))),1),2),5)),5,'lanczos2'));
im_source_roi4=abs(imresize(squeeze(sos(sos(sos(bsxfun(@times,abs(im_trace_adj_thd(:,:,:,:,2:6)),abs(BW_4Territories(:,:,4))),1),2),5)),5,'lanczos2'));

% Manually coregirster the souces images with MRA image
im_MRA=double(Braindata.im_MRA);
im_MRA=im_MRA./max(im_MRA(:));
im_MRA=padarray(padarray(im_MRA,[60 45],'pre'),[60 40],'post');
im_MRA=imresize(imrotate(im_MRA,4),[506 694]);
im_left_right=padarray(flip(im_source_roi2+im_source_roi3,3),[135 25],0,'pre');
im_left_right=im_left_right./max(im_left_right(:));
im_up=padarray(flip(im_source_roi1,3),[130 24],0,'pre');
im_up=im_up./max(im_up(:));
im_down=padarray(flip(im_source_roi4,3),[130 24],0,'pre');
im_down=im_down./max(im_down(:));

figure();
[Nx,Ny]=size(im_target_background);
cyellow=zeros(Nx,Ny,3);  cyellow(:,:,[1 2])=1;
cblue=zeros(Nx,Ny,3);   cblue(:,:,3)=1;
cgreen=zeros(Nx,Ny,3);  cgreen(:,:,2)=1;
cred=zeros(Nx,Ny,3);    cred(:,:,1)=1;
subplot(1,3,1);imshow(im_target_background.^1/2,[0 0.4]);title({'Displacement spectrum energy', 'masked with four vascular territories '});
hold on;hsource = imshow(cblue);hold off;
set(hsource , 'AlphaData', 0.3*abs(imresize(abs(BW_4Territories(:,:,1)),5,'bilinear')));
hold on;hsource = imshow(cred);hold off;
set(hsource , 'AlphaData', 0.3*abs(imresize(abs(BW_4Territories(:,:,2)+BW_4Territories(:,:,3)),5,'bilinear')));
hold on;hsource = imshow(cyellow);hold off;
set(hsource , 'AlphaData', 0.3*abs(imresize(abs(BW_4Territories(:,:,4)),5,'lanczos2')));

subplot(1,3,2);
[Nx,Ny]=size(im_up);
cyellow=zeros(Nx,Ny,3);cyellow(:,:,[1 2])=1;
cblue=zeros(Nx,Ny,3);cblue(:,:,3)=1;
cgreen=zeros(Nx,Ny,3);cgreen(:,:,2)=1;
cred=zeros(Nx,Ny,3);cred(:,:,1)=1;

imshow(im_MRA,[0.15 1.5]);title({'Spectrum images from red ROIs','overlaid on Coronal MRA '});
hold on;hsource = imshow(cblue);hold off;
set(hsource , 'AlphaData', 0.6*im_up);
hold on;hsource = imshow(cyellow);hold off;
set(hsource , 'AlphaData', 0.6*im_down);

subplot(1,3,3);
[Nx,Ny]=size(im_left_right);
cyellow=zeros(Nx,Ny,3);cyellow(:,:,[1 2])=1;
cblue=zeros(Nx,Ny,3);cblue(:,:,3)=1;
cgreen=zeros(Nx,Ny,3);cgreen(:,:,2)=1;
cred=zeros(Nx,Ny,3);cred(:,:,1)=1;
imshow(im_MRA,[0.15 1.5]);title({'Spectrum images from blue and yellow ROIs','overlaid on Coronal MRA '});
hold on;hsource = imshow(cred);hold off;
set(hsource , 'AlphaData', 0.6*im_left_right);
%% Let's look at the source images and their vascular territories at different mixing durations (Figure 9)
im_brain=Braindata.im_brain;
BW_source_rois=Braindata.BW_source_rois;
BW_source_rois(:,:,5)=circshift(BW_source_rois(:,:,5),[4 0 0]);
BW_brain_rmCSF=Braindata.BW_brain_rmCSF;
% Note that we remove some signals that have aliased along z axis and do T1weighting correction along Tmixing
im_trace_adj_mask=bsxfun(@times,im_trace_adj,BW_brain_rmCSF);
im_trace_adj_mask(:,:,[1:10 64:end],:,:)=0;
im_trace_adj_thd(:,:,[1:10 64:end],:,:)=0;
T1=0.8;
T1weight=1./exp(-(0.1:0.1:0.8)/T1);
im_trace_adj_mask=bsxfun(@times,im_trace_adj_mask,reshape(T1weight,[1 1 1 1 8]));

im_Vessels=squeeze(sos(sos(sos(bsxfun(@times,im_trace_adj_thd(:,:,:,:,2:6),BW_brain_rmCSF)),1),2));

[Nx,Ny]=size(im_brain);
cyellow=zeros(Nx,Ny,3);cyellow(:,:,[1 2])=1;
cgreen=zeros(Nx,Ny,3);cgreen(:,:,2)=1;
cred=zeros(Nx,Ny,3);cred(:,:,1)=1;

figure();
% ROI#1 in the source image and  their vascular territories at Tmixing=600,400,200ms
subplot(3,4,1);imshow(abs(imresize(crop(im_Vessels,[64 80]),[Nx Ny],'lanczos2')),[]); title('Spectrum images');
hold on;hsource = imshow(cred); 
set(hsource , 'AlphaData', 0.2*imresize(crop(BW_source_rois(:,:,1),[64 80]),[Nx Ny]));
hold on;hsource = imshow(cgreen); 
set(hsource , 'AlphaData', 0.2*imresize(crop(BW_source_rois(:,:,2),[64 80]),[Nx Ny]));
[N1,N2,~]=size(BW_source_rois);
im_target_1=squeeze(sos(sos(bsxfun(@times,im_trace_adj_mask,reshape(BW_source_rois(:,:,1),[1 1 N1 N2])),3),4));
im_target_2=squeeze(sos(sos(bsxfun(@times,im_trace_adj_mask,reshape(BW_source_rois(:,:,2),[1 1 N1 N2])),3),4));
im_target_1=2*im_target_1/max(im_target_1(:));
im_target_2=2*im_target_2/max(im_target_2(:));

subplot(3,4,2);imshow(im_brain,[]);title('Tmixing=600ms');
hold on;hsource = imshow(cred); 
set(hsource , 'AlphaData', imresize(im_target_1(:,:,6),[Nx Ny]));
hold on;hsource = imshow(cgreen); 
set(hsource , 'AlphaData', imresize(im_target_2(:,:,6),[Nx Ny]));

subplot(3,4,3);imshow(im_brain,[]);title('Tmixing=400ms');
hold on;hsource = imshow(cred); 
set(hsource , 'AlphaData', imresize(im_target_1(:,:,4),[Nx Ny]));
hold on;hsource = imshow(cgreen); 
set(hsource , 'AlphaData', imresize(im_target_2(:,:,4),[Nx Ny]));

subplot(3,4,4);imshow(im_brain,[]);title('Tmixing=200ms');
hold on;hsource = imshow(cred); 
set(hsource , 'AlphaData', imresize(im_target_1(:,:,2),[Nx Ny]));
hold on;hsource = imshow(cgreen); 
set(hsource , 'AlphaData', imresize(im_target_2(:,:,2),[Nx Ny]));

% ROI#2 in the source image and  their vascular territories at Tmixing=600,400,200ms
subplot(3,4,5);imshow(abs(imresize(crop(im_Vessels,[64 80]),[Nx Ny],'lanczos2')),[]); 
hold on;hsource = imshow(cred); 
set(hsource , 'AlphaData', 0.2*imresize(crop(BW_source_rois(:,:,3),[64 80]),[Nx Ny]));
hold on;hsource = imshow(cgreen); 
set(hsource , 'AlphaData', 0.2*imresize(crop(BW_source_rois(:,:,4),[64 80]),[Nx Ny]));
[N1,N2,~]=size(BW_source_rois);
im_target_3=squeeze(sos(sos(bsxfun(@times,im_trace_adj_mask,reshape(BW_source_rois(:,:,3),[1 1 N1 N2])),3),4));
im_target_4=squeeze(sos(sos(bsxfun(@times,im_trace_adj_mask,reshape(BW_source_rois(:,:,4),[1 1 N1 N2])),3),4));
im_target_3=2*im_target_3/max(im_target_3(:));
im_target_4=2*im_target_4/max(im_target_4(:));

subplot(3,4,6);imshow(im_brain,[]);
hold on;hsource = imshow(cred); 
set(hsource , 'AlphaData', imresize(im_target_3(:,:,6),[Nx Ny]));
hold on;hsource = imshow(cgreen); 
set(hsource , 'AlphaData', imresize(im_target_4(:,:,6),[Nx Ny]));

subplot(3,4,7);imshow(im_brain,[]);
hold on;hsource = imshow(cred); 
set(hsource , 'AlphaData', imresize(im_target_3(:,:,4),[Nx Ny]));
hold on;hsource = imshow(cgreen); 
set(hsource , 'AlphaData', imresize(im_target_4(:,:,4),[Nx Ny]));

subplot(3,4,8);imshow(im_brain,[]);
hold on;hsource = imshow(cred); 
set(hsource , 'AlphaData', imresize(im_target_3(:,:,2),[Nx Ny]));
hold on;hsource = imshow(cgreen); 
set(hsource , 'AlphaData', imresize(im_target_4(:,:,2),[Nx Ny]));


% ROI#3 in the source image and their vascular territories at Tmixing=600,400,200ms
subplot(3,4,9);imshow(abs(imresize(crop(im_Vessels,[64 80]),[Nx Ny],'lanczos2')),[]); 
hold on;hsource = imshow(cyellow); 
set(hsource , 'AlphaData', 0.2*imresize(crop(BW_source_rois(:,:,5),[64 80]),[Nx Ny]));

im_target_5=squeeze(sos(sos(bsxfun(@times,im_trace_adj_mask,reshape(BW_source_rois(:,:,5),[1 1 N1 N2])),3),4));
im_target_5=im_target_5.^2;
im_target_5=2*im_target_5/max(im_target_5(:));



subplot(3,4,10);imshow(im_brain,[]);
hold on;hsource = imshow(cyellow); 
set(hsource , 'AlphaData', imresize((im_target_5(:,:,6)),[Nx Ny]));

subplot(3,4,11);imshow(im_brain,[]);
hold on;hsource = imshow(cyellow); 
set(hsource , 'AlphaData', imresize((im_target_5(:,:,4)),[Nx Ny]));

subplot(3,4,12);imshow(im_brain,[]);
hold on;hsource = imshow(cyellow); 
set(hsource , 'AlphaData', imresize((im_target_5(:,:,2)),[Nx Ny]));
