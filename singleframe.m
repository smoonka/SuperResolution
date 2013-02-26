% This file does single frame SR for LLE 
% It receives as input a single image, it considers this to be the HR.
% It downsamples this, and generates LR-HR pairs of patches.
% Then does k-medioid on the patches and generates a dictionary for LLE
% Then LLE is done on the original input image to upsample it further

% do_config;

% maxNumCompThreads(1);
clear all;
Dhm = [];
Dlm = [];
G = [];
clear all;
Dhm = [];
Dlm = [];
G = [];
HP=[];
LP=[];
n = 2;
% =====================================================================
% specify the parameter settings

patch_size = 3; % patch size for the low resolution input image
overlap = 2; % overlap between adjacent patches
lambda = 0.1; % sparsity parameter
zooming = 3; % zooming factor, if you change this, the dictionary needs to be retrained.

codebook_size = 1024; % size of the dictionary

% =====================================================================
%filename = '/home/smoonka/Downloads/lle_sr/test_run3/training_image.png';
for i =2:n
str = int2str(i);
filename = '/home/smoonka/Downloads/lle_sr/';
%filename = '/home/smoonka/Downloads/lle_sr/test11withsameCameradifferentframes/';
if n< 10
calibdata = strcat(filename,'000',str,'.cam');
filename = strcat(filename,'000',str,'.png');
else
calibdata = strcat(filename,'00',str,'.cam');
filename = strcat(filename,'00',str,'.png');
end
image = imread(filename);
ratio = alignment(calibdata, image);
image = imresize(image,1/ratio);
ratio = 1;
[HPi,LPi] = getpatches(image,patch_size, zooming,ratio);
if (size(HPi,2)< codebook_size) 
HP = [HP HPi];
LP = [LP LPi];
if (size(HP,2)> codebook_size) 
[Dhi, Dli] = coupled_dic_train(HP, LP, codebook_size, lambda);
HP=[];
LP=[];
end
else
[Dhi, Dli] = coupled_dic_train(HPi, LPi, codebook_size, lambda);
end
Dhm = [Dhm,Dhi];
Dlm = [Dlm,Dli];
end

save Dhm.mat Dhm;
save Dlm.mat Dlm;

%fullimage_fn='/home/smoonka/Downloads/lle_sr/lrimagetest.png';
%fullimage_fn='/home/smoonka/Downloads/lle_sr/testsetwithgoodimages/test.png';
fullimage_fn='/home/smoonka/Downloads/lle_sr/HHI_180_ChelseaWolves_Clip1_0000.png';
%fullimage_fn='/home/smoonka/Downloads/lle_sr/0020.png';
fullimage = imread(fullimage_fn);
%lr_color = imresize(fullimage, 1/zooming, 'bicubic');
%imwrite(lr_color,'lr.png');

%image = imread('lr.png');
image = fullimage;
%image = image(500:1300,3000:5000,:);
image = image(500:1100,4000:4500,:);
image(:,:,1)=fliplr(image(:,:,1));
image(:,:,2)=fliplr(image(:,:,2));
image(:,:,3)=fliplr(image(:,:,3));
%image = imresize(image,1/zooming);
bicubic = imresize(image, zooming,'bicubic');
%bicubic = imread(bicubic_name);
if( size(image,3) == 3)
    interp2=rgb2ycbcr(bicubic);
    hImcb = interp2(:,:,2);
    hImcr = interp2(:,:,3);
    lr2=rgb2ycbcr(image);
    lr = lr2(:,:,1);

    
%     interp2=rgb2hsv(bicubic);
%     hImcb = interp2(:,:,1);
%     hImcr = interp2(:,:,2);
%     lr2=rgb2hsv(image);
%     lr = lr2(:,:,3);
else
    lr = image;
end;
lr = double(lr);
[hImy] = sparse_coding(lr, zooming, patch_size, overlap, Dhm, Dlm,'L1',n);

% imwrite(hImy,'greyresult.png');
if( size(image,3) == 3)
    
    hr = bicubic;
    hr(:,:,1) = uint8(hImy);
    hr(:,:,2) = hImcb;
    hr(:,:,3) = hImcr;
    hr = ycbcr2rgb(hr);
% 
%     
%     
%     hr = bicubic;
%     hr(:,:,3) = (hImy);
%     hr(:,:,1) = hImcb;
%     hr(:,:,2) = hImcr;
%     %hr = ycbcr2rgb(hr)
%     hr = hsv2rgb(hr);
else
    hr = hImy;
end;

imwrite(hr,'lle_resultp.png');
imwrite(bicubic,'bicubic_resultp.png');


%addpath('Solver');
%addpath('Sparse coding');

% training coupled dictionaries for super-resolution

% if ~skip_smp_training,
%     disp('Sampling image patches...');
%     [Xh, Xl] = rnd_smp_dictionary(tr_dir, patch_size, zooming, num_patch);
%     save('Data/Dictionary/smp_patches.mat', 'Xh', 'Xl');
%     skip_dictionary_training = false;
% end;
% 
% if ~skip_dictionary_training,
%     load('/users/visics/vnambood/programs/sr/sr_cs/lle_sr/Data/Dictionary/smp_patches.mat');
%     [Dh, Dl] = coupled_dic_train(Xh, Xl, codebook_size, lambda);
%     save('Data/Dictionary/Dictionary.mat', 'Dh', 'Dl');
% else
%     load('/users/visics/vnambood/programs/sr/sr_cs/lle_sr/Data/Dictionary/Dictionary.mat');
% end;

% =====================================================================
% Process the test image 


% startval=str2num(start_index);
% endval=str2num(end_index);
% logfn=sprintf('/esat/hadar/vnambood/log/sr/lle_sr/tugraz_up_logfile_%d_%d',startval,endval);
% %logfn=sprintf('/users/visics/vnambood/programs/sr/sr_cs/lle_sr/logfile_%d_%d',startval,endval);
% display(startval);
% display(endval);
% display(logfn);
% fid=fopen(logfn,'w');
% database = retr_database_dir(database_path);
% %save('database.mat','database');
% nFea = length(database.path);
% for ii=startval:endval
%     display(ii);
%     fprintf(fid,'id = %d\n',ii);
%     img_fn=database.path{ii};
%     fprintf(fid,'fn = %s\n',img_fn);
%     display(img_fn);
%     lr_img=imread(img_fn);
%     bicubic = imresize(lr_img, zooming,'bicubic');
%     if( size(lr_img,3) == 3)
%         interp2=rgb2ycbcr(bicubic);
%         hImcb = interp2(:,:,2);
%         hImcr = interp2(:,:,3);
%         lr2=rgb2ycbcr(lr_img);
%         lr = lr2(:,:,1);
%     else
%         lr = lr_img;
%     end;
%     lr = double(lr);
%     %[hImy] = L1SR(lr, zooming, patch_size, overlap, Dh, Dl, lambda, regres);
%     [hImy] = do_SR(lr, zooming, patch_size, overlap, Dh, Dl,'L1');
%     if( size(lr_img,3) == 3)
%        hr = bicubic;
%         hr(:,:,1) = uint8(hImy);
%         hr(:,:,2) = hImcb;
%         hr(:,:,3) = hImcr;
%         hr = ycbcr2rgb(hr);
%     else
%         hr = hImy;
%     end;
%     
%     [dir, fn, ext] = fileparts(img_fn);
%     wdir=strrep(dir,database_path,result_path); 
% %     bdir=strrep(dir,database_path,bicubic_path);
% %     bic_fn=fullfile(bdir,[fn ext]);
% %     imwrite(bicubic,bic_fn);
%     hr_fn = fullfile(wdir,[fn ext]);
%     imwrite(hr,hr_fn);
% end;
% fclose(fid);
