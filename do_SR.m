function [hIm, ww] = do_SR(lIm, zooming, patch_size, overlap, Dh, Dl,regres,n)
% Use sparse representation as the prior for image super-resolution
% Usage
%       [hIm] = L1SR(lIm, zooming, patch_size, overlap, Dh, Dl, lambda)
% 
% Inputs
%   -lIm:           low resolution input image, single channel, e.g.
%   illuminance
%   -zooming:       zooming factor, e.g. 3
%   -patch_size:    patch size for the low resolution image
%   -overlap:       overlap among patches, e.g. 1
%   -Dh:            dictionary for the high resolution patches
%   -Dl:            dictionary for the low resolution patches
%   -regres:       'L1' use the sparse representation directly to high
%                   resolution dictionary;
%                   'L2' use the supports found by sparse representation
%                   and apply least square regression coefficients to high
%                   resolution dictionary.
% Ouputs
%   -hIm:           the recovered image, single channel
%


[lhg, lwd] = size(lIm);
hhg = lhg*zooming;
hwd = lwd*zooming;

%tmpb = imread(bicubic_name);
%tmpb2=rgb2hsv(tmpb);
%mIm = tmpb2(:,:,3);
%mIm = double(mIm);

mIm = imresize(lIm, zooming,'bicubic');
[mhg, mwd] = size(mIm);
hpatch_size = patch_size*zooming;
mpatch_size = patch_size*zooming;

% extract gradient feature from lIm
hf1 = [-1,0,1];
vf1 = [-1,0,1]';
hf2 = [1,0,-2,0,1];
vf2 = [1,0,-2,0,1]';

lImG11 = conv2(lIm,hf1,'same');
lImG12 = conv2(lIm,vf1,'same');
lImG21 = conv2(lIm,hf2,'same');
lImG22 = conv2(lIm,vf2,'same');

lImG11 = imresize(lImG11, zooming,'bicubic');
lImG12 = imresize(lImG12, zooming,'bicubic');
lImG21 = imresize(lImG21, zooming,'bicubic');
lImG22 = imresize(lImG22, zooming,'bicubic');

lImfea(:,:,1) = lImG11;
lImfea(:,:,2) = lImG12;
lImfea(:,:,3) = lImG21;
lImfea(:,:,4) = lImG22;

lgridx = 2:patch_size-overlap:lwd-patch_size;
lgridx = [lgridx, lwd-patch_size];
lgridy = 2:patch_size-overlap:lhg-patch_size;
lgridy = [lgridy, lhg-patch_size];

mgridx = (lgridx - 1)*zooming + 1;
mgridy = (lgridy - 1)*zooming + 1;

% using linear programming to find sparse solution
bhIm = imresize(lIm, zooming, 'bicubic');
%bhIm = mIm;
hIm = zeros([hhg, hwd]);
nrml_mat = zeros([hhg, hwd]);

hgridx = (lgridx-1)*zooming + 1;
hgridy = (lgridy-1)*zooming + 1;

disp('Processing the patches sequentially...');
count = 0;

% loop to recover each patch
for xx = 1:length(mgridx),
    for yy = 1:length(mgridy),
        
        mcolx = mgridx(xx);
        mrowy = mgridy(yy);
        
        count = count + 1;
        if ~mod(count, 1000),
            fprintf('.\n');
        end;
        if ~mod(count, 100),
            fprintf('.');
        end;
        
        mpatch = mIm(mrowy:mrowy+mpatch_size-1, mcolx:mcolx+mpatch_size-1);

        mmean = mean(mpatch(:));
        
        mpatchfea = lImfea(mrowy:mrowy+mpatch_size-1, mcolx:mcolx+mpatch_size-1, :);
        mpatchfea = mpatchfea(:);
        

        mnorm = sqrt(sum(mpatchfea.^2));
        
        if mnorm > 1,
            y = mpatchfea./mnorm;
        else
            y = mpatchfea;
        end;

        
        %w = SolveLasso(Dl, y, size(Dl, 2), 'nnlasso', [], lambda);
  	  w = doProject(y, Dl',n);
%         w = feature_sign(Dl, y, lambda);
        
        if isempty(w),
            display('w is zero');
            display(xx);
            display(yy);
            w = zeros(size(Dl, 2), 1);
        end;

        switch regres,
            case 'L1'
                if mnorm > 1,
                    hpatch = Dh*w'*mnorm;
                else
                    hpatch = Dh*w';
                end;
            case 'L2'
                idx = find(w);
                lsups = Dl(:, idx);
                hsups = Dh(:, idx);
                w = inv(lsups'*lsups)*lsups'*mpatchfea;
                hpatch = hsups*w;
            otherwise
                error('Unknown fitting!');
        end;
      
% 		  display(size(hpatch));
        hpatch = reshape(hpatch, [hpatch_size, hpatch_size]);
        hpatch = hpatch + mmean;
        
        hcolx = hgridx(xx);
        hrowy = hgridy(yy);
        
        hIm(hrowy:hrowy+hpatch_size-1, hcolx:hcolx+hpatch_size-1)...
            = hIm(hrowy:hrowy+hpatch_size-1, hcolx:hcolx+hpatch_size-1) + hpatch;
        nrml_mat(hrowy:hrowy+hpatch_size-1, hcolx:hcolx+hpatch_size-1)...
            = nrml_mat(hrowy:hrowy+hpatch_size-1, hcolx:hcolx+hpatch_size-1) + 1;
fprintf('out................ \n');
xx
yy

    end;
end;

fprintf('done!\n');

% fill the empty
hIm(1:3, :) = bhIm(1:3, :);
hIm(:, 1:3) = bhIm(:, 1:3);

hIm(end-2:end, :) = bhIm(end-2:end, :);
hIm(:, end-2:end) = bhIm(:, end-2:end);

nrml_mat(nrml_mat < 1) = 1;
hIm = hIm./nrml_mat;
%hIm = uint8(hIm);

