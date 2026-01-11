function [image_feature] = PrcfeatureExtraction(Ns,idx_sup,image_vector) 
[hw,b] = size(image_vector);
image_feature = zeros(Ns,4*b);
for i = 1:Ns
    index_vector = idx_sup{i};
    sub_superpixel = image_vector(index_vector,:);
    mean_feature = (mean(sub_superpixel));
    median_feature = (median(sub_superpixel));
    m75_feature = (prctile(sub_superpixel,75));
    m25_feature = (prctile(sub_superpixel,25));
   
    image_feature(i,1:b) = mean_feature;
    image_feature(i,b+1:2*b) = m25_feature;
    image_feature(i,2*b+1:3*b) = median_feature;
    image_feature(i,3*b+1:4*b) = m75_feature;
end
image_feature(find(isnan(image_feature)==1)) = 0;    
