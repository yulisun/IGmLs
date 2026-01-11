function [sup_pixel,N] = GMMSP_segmentation(image_t1,seg_scal)
[h, w, b1]=size(image_t1);
v_x = floor(sqrt(h*w/seg_scal));
if v_x < 8
    v_x = 8;
end
v_y = v_x;
new_image = zeros(h, w,3);
if b1 ==3
    new_image = image_t1;
elseif b1==1
    new_image(:,:,1) = image_t1;
    new_image(:,:,2) = image_t1;
    new_image(:,:,3) = image_t1;
elseif b1==2
    new_image(:,:,1) = image_t1(:,:,1);
    new_image(:,:,2) = image_t1(:,:,2);
elseif b1>3
    [pc,score,latent,tsquare] = pca(double(image_t1),'NumComponents',3);
    for i = 1:3
        tmep = score(:,i);
        new_image(:,:,i) = reshape(tmep,[h w]);
    end
end
new_image = im2uint8(new_image);
sup_pixel = mx_GMMSP(new_image, v_x, v_y);
N = max(sup_pixel(:));