function [HP, LP ] = getpatches( im, patch_size, zooming,ratio )
%   Detailed explanation goes here
% Function to get patches from the dictionary images

%Inputs
%   -im:           High Resolution image from dictionary of images
%   -patch_size:   patch size for the low resolution image
%   -zooming:      Desired super resolution ratio
%   -ratio:        Downsampling ratio between the 2 cameras computed after
%   alignment procedure
% Ouputs
%   -HP:           High Resolution Patch
%   -LP:           Low  Resolution Patch


lz = zooming;

if size(im, 3) == 3,
    hIm = rgb2gray(im);
else
    hIm = im;
end;

if rem(size(hIm,1),zooming)
    nrow = floor(size(hIm,1)/zooming)*zooming;
    hIm = hIm(1:nrow,:);
end;
if rem(size(hIm,2),zooming)
    ncol = floor(size(hIm,2)/zooming)*zooming;
    hIm = hIm(:,1:ncol);
end;
%H = fspecial('gaussian',7,2);
%lIm = imfilter(hIm,H,'same');
lIm = imresize(hIm,1/(zooming*ratio));
% blurring can be done here

[nrow, ncol] = size(lIm);

x = randperm(nrow-patch_size-lz-1);
y = randperm(ncol-patch_size-lz-1);
[X,Y] = meshgrid(x,y);

xrow = X(:);
ycol = Y(:);


%display(size(xrow));
%display(size(ycol));

patch_num = 20000;

if(length(xrow) < patch_num )
    patch_num = length(xrow);
end;
%xrow = xrow(1:patch_num);
%ycol = ycol(1:patch_num);

% zoom the original image
lIm = imresize(lIm, lz,'bicubic');
hIm = double(hIm);
lIm = double(lIm);

HP = zeros(zooming^2*patch_size^2,patch_num);
LP = zeros(lz^2*4*patch_size^2,patch_num);

%HP = zeros(patch_size^2,patch_num);
%LP = zeros(4*patch_size^2,patch_num); 
% compute the first and second order gradients
hf1 = [-1,0,1];
vf1 = [-1,0,1]';
 
lImG11 = conv2(lIm,hf1,'same');
lImG12 = conv2(lIm,vf1,'same');
 
hf2 = [1,0,-2,0,1];
vf2 = [1,0,-2,0,1]';
 
lImG21 = conv2(lIm,hf2,'same');
lImG22 = conv2(lIm,vf2,'same');

count = 1;
for pnum = 1:patch_num,
    
    hrow = (xrow(pnum)-1)*zooming + 1;
    hcol = (ycol(pnum)-1)*zooming + 1;
    Hpatch = hIm(hrow:hrow+zooming*patch_size-1,hcol:hcol+zooming*patch_size-1);
    
    lrow = (xrow(pnum)-1)*lz + 1;
    lcol = (ycol(pnum)-1)*lz + 1;
    
    Lpatch1 = lImG11(lrow:lrow+lz*patch_size-1,lcol:lcol+lz*patch_size-1);
    Lpatch2 = lImG12(lrow:lrow+lz*patch_size-1,lcol:lcol+lz*patch_size-1);
    Lpatch3 = lImG21(lrow:lrow+lz*patch_size-1,lcol:lcol+lz*patch_size-1);
    Lpatch4 = lImG22(lrow:lrow+lz*patch_size-1,lcol:lcol+lz*patch_size-1);
     
    Lpatch = [Lpatch1(:),Lpatch2(:),Lpatch3(:),Lpatch4(:)];
    Lpatch = Lpatch(:);
     
    HP(:,count) = Hpatch(:)-mean(Hpatch(:));
    LP(:,count) = Lpatch;
    
    count = count + 1;
 end
fprintf('patching_done!\n');
