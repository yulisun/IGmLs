clear all
close all
warning('off')
addpath(genpath(pwd))
%% load dataset

str1 = 'Dataset2_im1.bmp';
str2 = 'Dataset2_im2.bmp';
strgt = 'Dataset2_gt.bmp';
image_1 = imread(str1); %Image 1 serves as the reference image, on which the changes are overlaid
image_2 = imread(str2);
Ref_gt = imread(strgt); % Note that the ground-truth map is annotated with Image 1 as the reference.
figure;
subplot(131);imshow(image_1);title('image1')
subplot(132);imshow(image_2);title('image2')
subplot(133);imshow(Ref_gt,[]);title('Refgt')
% Use a checkerboard map to visualize the registration error.
grid_num=8;
grid_size=floor(min(size(image_1,1),size(image_1,2))/grid_num);
[~,~,img3] = mosaic_map(image_1,image_2,grid_size);
figure; imshow(img3,[]);title('Fused image of the board');
fprintf(['\n Data loading is completed...... ' '\n'])
%% Parameter setting

opt.type_1 = 'optical'; % optical or sar, based on the type of image 1 
opt.type_2 = 'optical'; % optical or sar, based on the type of image 2
opt.Nf = 2500; 
opt.Nc = 500;
opt.Niter = 2;
opt.wh = 18; % the height of the search window
if strcmp(str1,'Dataset1_im1.bmp') == 1
    opt.wh = 18; 
elseif strcmp(str1,'Dataset2_im1.bmp') == 1
    opt.wh = 9;
end
opt.ww = opt.wh; % the width of the search window
opt.ws = 3; % the search stride. A larger ws can improve efficiency but may degrade accuracy
opt.alfa = 0.01; 
opt.gama = 2; 
opt.kmax = round(sqrt(opt.Nc));
%% IGmLs

t_o = clock;
fprintf(['\n IGmLs is running...... ' '\n'])
[CM_iter,DI_iter] = IGmLs_main(image_1,image_2,opt);
fprintf(['\n' '====================================================================== ' '\n'])
fprintf('\n');fprintf('The total computational time of IGmLs is %.3f s\n', etime(clock, t_o));
fprintf(['\n' '====================================================================== ' '\n'])
%% Displaying results

fprintf(['\n Displaying the results...... ' '\n'])
Ref_gt = double(Ref_gt(:,:,1));
Ref_gt = Ref_gt/max(Ref_gt(:));
for i= 1:opt.Niter
    [tp,fp,tn,fn,fplv,fnlv,~,~,OA(i),kappa(i),imw]=performance(CM_iter(:,:,i),Ref_gt);
    F1(i) = 2*tp/(2*tp+fp+fn);
    fprintf('Iteration %d: OA = %.3f; Kc = %.3f; F1 = %.3f\n', i, OA(i), kappa(i), F1(i));
end
figure;
for iter = 1:opt.Niter    
    subplot(opt.Niter,3,1+3*(iter-1));imshow(remove_outlier(DI_iter(:,:,iter)),[]),title(['DI of the' num2str(iter) '-th iteration'])
    subplot(opt.Niter,3,2+3*(iter-1));imshow(CMplotRGB(CM_iter(:,:,iter),Ref_gt)),title(['CM of the' num2str(iter) '-th iteration'])
    bw = boundarymask(CM_iter(:,:,iter),4);
    img1_CM = imoverlay(image_1,bw,[1 0 0]);
    subplot(opt.Niter,3,3+3*(iter-1));imshow(img1_CM),title(['CM of the' num2str(iter) '-th iteration overlaid on image 1'])
end
