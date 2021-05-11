function [im,rdrFOV]=DiSpect_shearing(im,r_dr_dims,r_dr_FOVs)
% In order to achiveing p(...,r,...,dr)-->p(...,r,...r+dr)
%          fft             exp(1i*kd*r)                           ifft
% p(r,d)-------> m(r,kd) ---------------> m(r,kd)*exp(1i*kd*r) --------->p(r,r+d)

%         zpad             fft             exp(1i*kd*r)                           ifft
% p(r,d)------->prd(r,d)-------> m(r,kd) ---------------> m(r,kd)*exp(1i*kd*r) --------->p(r,r+d)

%
% input:  im  ->  multiple dimensional image dataset
%         r_dr_dims -> indicates the dimensions of r (regular imaging dimension)
%                      and dr (displacement encoded dimension) in the input image;
%         r_dr_FOVs -> indicates the FOV information for r (regular imaging dimension)
%                      and dr (displacement encoded dimension) in the input image;

% Output: im_shearing  ->  image after this shearing transformation
%        rdrFOVs   ->  FOV information in absolute physical coordinate;

% (c) zhiyong zhang, 2021, zhiyongxmu@gmail.com

%%
rdim=r_dr_dims(1);
drdim=r_dr_dims(2);

rFOV=r_dr_FOVs(1);
drFOV=r_dr_FOVs(2);
rdrFOV=abs(rFOV)+abs(drFOV);

Npad=2*round(size(im,rdim)*abs(rdrFOV/rFOV)/2);
rdrFOV=Npad*drFOV/size(im,drdim);
kFOV=Npad*1/rdrFOV;

imsize=size(im);
imsize(drdim)=Npad;
im=zpad(im,imsize);

kim=cfft(im,drdim);
Nk=size(kim,drdim);
Nr=size(kim,rdim);

k=reshape((-Nk/2:Nk/2-1)*kFOV/Nk,[ones(1,drdim-1) Nk 1]);
r=reshape((-Nr/2:Nr/2-1)*rFOV/Nr,[ones(1,rdim-1) Nr 1]);
ph=exp(-1i*2*pi*bsxfun(@times,k,r));

kimph=bsxfun(@times,kim,ph);
im=cifft(kimph,drdim);


