%% This is a demo processing of DiSpect for phantom II
% The rawdata set: [Nro,Npe,Nd_LR,Nd_SI,N_mixing, Nc,Nphasecycling];
% The FOV information: [FOVro, FOVpe, FOVd_LR,FOVd_SI]; 
load PhantomII;
%% Displacement encoded echo seperation from 3-step phase cycling
PhaseCyclingEnc3=[       1,              1,          1;
                    exp(-1i*pi*2/3), exp(1i*pi*2/3), 1;
                    exp(1i*pi*2/3), exp(-1i*pi*2/3), 1];             
PhaseCyclingDec3=inv(PhaseCyclingEnc3);
ksp_raw=PhantomII.ksp_raw;
dSize=size(ksp_raw);
kspDiSpect=reshape(ksp_raw,[],3)*PhaseCyclingDec3.';
kspDiSpect=reshape(kspDiSpect(:,1),dSize(1:end-1));
%% FT operation 
%  Window k space to reduce the gibbs ringing
[Nro,Npe,Nd_LR,Nd_SI,N_mixing, Nc]=size(kspDiSpect);
kwin1=1;
kwin2=hamming(Nd_LR)*hamming(Nd_SI).';
kwin=bsxfun(@times,kwin1,permute(kwin2,[3 4 1 2]));
im_DiSpect=cifftn(bsxfun(@times,kspDiSpect,kwin),[1 2 3 4]);

%% Shearing transformation to obtain the absolute physical coordinates from displacement spectra 
% The shearing transformation is done along the first and third dimensions
% Note that we did a circshift to make z=0 at the top of the matrix
im_DiSpect=circshift(im_DiSpect,[0 0 0 Nd_SI/2 0]);
im_trace=sos(DiSpect_shearing(im_DiSpect,[1 3],[PhantomII.FOVinf(1),PhantomII.FOVinf(3)]));
%% Let's look at the source images at different mixing times by summing up all the image voxels
BWtarget=PhantomII.BWtarget;
im_trace_show=im_trace;
im_source=squeeze(sos(sos(bsxfun(@times,im_trace_show,BWtarget),1),2));
Nx=256;
Ny=160;
im_source_3mixing=flip(flip(permute(imresize(im_source(:,3:end-4,[1 3 5]),[Nx Ny]),[2 1 3]),1),2);
im_source_3mixing=im_source_3mixing./repmat(max(max(im_source_3mixing,[],1),[],2),[Ny Nx 1]);
figure();imshowMRI(im_source_3mixing,[0,1;0.02 1;0.25 1],[],{'Tmixing=100 ms','Tmixing=200 ms','Tmixing=300 ms'})
%% Let's look at the target images from 2 sections at Tmixing=150
[~,~,Nd1,Nd2,~]=size(im_trace);
im_phantom=PhantomII.im_phantom;
BW_2section=PhantomII.BW_2section;
im_target1=squeeze(sos(sos(bsxfun(@times,im_trace(:,:,:,:,2),reshape(BW_2section(:,:,1),[1 1 Nd1,Nd2])),3),4));
im_target2=squeeze(sos(sos(bsxfun(@times,im_trace(:,:,:,:,2),reshape(BW_2section(:,:,2),[1 1 Nd1,Nd2])),3),4));

Nx=256;
Ny=256;
im_target1=flip(imresize(im_target1.',[Nx Ny]),2);
im_target2=flip(imresize(im_target2.',[Nx Ny]),2);
im_target1=im_target1./max(im_target1(:));
im_target2=im_target2./max(im_target2(:));

cgreen=zeros(Nx,Ny,3);
cgreen(:,:,2)=1;
cred=zeros(Nx,Ny,3);
cred(:,:,1)=1;
figure();
imagesc(im_phantom) ;title('Color-coded targeted image from 2 sections');axis off;
hold on;h1=imshow(cred); set(h1 , 'AlphaData',im_target1);
hold on;h2=imshow(cgreen); set(h2 , 'AlphaData',im_target2);

%% Let's look at the target images from 4 sections at Tmixing=150
BW_4section=PhantomII.BW_4section;
im_target1=squeeze(sos(sos(bsxfun(@times,im_trace(:,:,:,:,2),reshape(BW_4section(:,:,1),[1 1 Nd1 Nd2])),3),4));
im_target2=squeeze(sos(sos(bsxfun(@times,im_trace(:,:,:,:,2),reshape(BW_4section(:,:,2),[1 1 Nd1 Nd2])),3),4));
im_target3=squeeze(sos(sos(bsxfun(@times,im_trace(:,:,:,:,2),reshape(BW_4section(:,:,3),[1 1 Nd1 Nd2])),3),4));
im_target4=squeeze(sos(sos(bsxfun(@times,im_trace(:,:,:,:,2),reshape(BW_4section(:,:,4),[1 1 Nd1 Nd2])),3),4));

Nx=256;
Ny=256;
im_target1=flip(imresize(im_target1.',[Nx Ny]),2);
im_target2=flip(imresize(im_target2.',[Nx Ny]),2);
im_target3=flip(imresize(im_target3.',[Nx Ny]),2);
im_target4=flip(imresize(im_target4.',[Nx Ny]),2);
im_target1=im_target1./max(im_target1(:));
im_target2=im_target2./max(im_target2(:));
im_target3=im_target3./max(im_target3(:));
im_target4=im_target4./max(im_target4(:));


cyellow=zeros(Nx,Ny,3);
cyellow(:,:,[1 2])=1;
cblue=zeros(Nx,Ny,3);
cblue(:,:,3)=1;
cgreen=zeros(Nx,Ny,3);
cgreen(:,:,2)=1;
cred=zeros(Nx,Ny,3);
cred(:,:,1)=1;

figure();
imagesc(im_phantom) ;title('Color-coded targeted image from 4 sections');axis off;
hold on;h1=imshow(cyellow); set(h1 , 'AlphaData',0.8*im_target1);
hold on;h2=imshow(cblue); set(h2 , 'AlphaData',0.8*im_target2);
hold on;h3=imshow(cgreen); set(h3 , 'AlphaData',0.8*im_target3);
hold on;h4=imshow(cred); set(h4 , 'AlphaData',0.8*im_target4);