function [ a ] = see4_DiSpect(IM, noabs, rescale, figpointer, count)

use_abs = false;
if nargin < 2 || noabs == 0
    use_abs = true;
end

if nargin < 3
    rescale = false;
end

if nargin < 4
    figpointer = 22;
end

if nargin < 5
    count = inf
end

button=1;
im = sqrt(sum(abs(IM(:,:,:).^2),3));


if rescale
    imsort = sort(abs(IM(:)), 'ascend');
    ymin = imsort(1);
    ymax = imsort(round(end - end/1000));
end

figmain = figpointer;
figpointer = figpointer + 1;

h0=figure(figmain); imshow(im,[],'InitialMagnification','fit'), impixelinfo

%     h1 = imellipse;
%     BW1 = createMask(h1);
%     h2 = imellipse;
%     BW2 = createMask(h2);
%     
while count > 0
    count = count - 1;
    figure(figmain)

%     h = impoly();
%     h=imrect();
    h=drawpolygon;
    BW = createMask(h);
    
    a=squeeze(sum(sum(bsxfun(@times,IM,logical(BW)),1),2));
    
%     imS1=cat(3,a(:,:,1),a(:,:,3).^2,sum(a(:,:,5).^(4),3));
%     win1=hamming(size(imS1,1)).^(0);
%     win2=hamming(size(imS1,2)).^(1);
%     imS2=abs(cfftn(bsxfun(@times,cifftn(imS1,[1 2]),win1*win2.'),[1 2]));
%     figure();imshowMRI(imresize(crop(imS2(:,1:24,:),[48 24 3]),[192 96]).^(1),0);
    
%     figure(figpointer), 
    imSize=size(a);
    imSize=imSize(1:2);
%     if imSize(1)<64 || imSize(2)<64
%         imSize=imSize(1:2)*10;
%     end
%     imshow3(imresize(a,imSize));
    figure();imshowMRI(cat(3,a(:,:,:),sos(a(:,:,:))),0)
%      imshowMRI(imresize(crop(imS2(:,1:24,:),[48 24 3]),[192 96]).^(1),0);
%     ylim([0 0.02*max(sum(a(:,:),2))]);
    
%     title(sprintf('(%d, %d)',round(y),round(x)));
    if button == 3
        figpointer = figpointer + 1;
    end

end