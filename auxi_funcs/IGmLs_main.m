function [CM_map_iter,DI_map_iter] = IGmLs_main(image_1,image_2,opt)
%% pre-processing
t_pre = clock;
image_1 = double(image_1);
image_2 = double(image_2);
ht1 = fspecial('average',5);
image_1 = imfilter(image_1, ht1,'symmetric');
image_2 = imfilter(image_2, ht1,'symmetric');
[image_1,~] = image_normlized(image_1,opt.type_1);
[image_2,~] = image_normlized(image_2,opt.type_2);
image_1 = image_1(opt.wh+1:end-opt.wh,opt.ww+1:end-opt.ww,:);
[Cosup_fine,Nf] = GMMSP_segmentation(image_1,opt.Nf);
[Cosup_coarse,Nc] = GMMSP_segmentation(image_1,opt.Nc);
idxsup_fine = label2idx(Cosup_fine);
idxsup_coarse = label2idx(Cosup_coarse);
[ht1, wt1, bt1] = size(image_1);
imgt1b = reshape(image_1,[],bt1);
t1feature_fine = zeros(Nf,4*bt1);
[t1feature_fine] = PrcfeatureExtraction(Nf,idxsup_fine,imgt1b) ;
[t1feature_coarse] = PrcfeatureExtraction(Nc,idxsup_coarse,imgt1b) ;
fprintf('\n');fprintf('The computational time of Preprocessing (t_pre) is %i \n', etime(clock, t_pre));
%% Img2 feature extraction and KNN template
t_template = clock;
[ht2, wt2, bt2] = size(image_2);
imgt2b = reshape(image_2,[],bt2);
img2_o = image_2(opt.wh+1:opt.wh+ht1,opt.ww+1:opt.ww+wt1,:);
img2_o = reshape(img2_o,[],bt2);
[t2feature_coarse] = PrcfeatureExtraction(Nc,idxsup_coarse,img2_o) ;
idsort_t2 = cell(floor(2*opt.wh/opt.ws)+1,floor(2*opt.ww/opt.ws)+1);
dist_t2 = cell(floor(2*opt.wh/opt.ws)+1,floor(2*opt.ww/opt.ws)+1);
t2_feature = cell(floor(2*opt.wh/opt.ws)+1,floor(2*opt.ww/opt.ws)+1);
f = zeros(floor(2*opt.wh/opt.ws)+1,floor(2*opt.ww/opt.ws)+1,Nf);
for dh = -opt.wh:opt.ws:opt.wh
    for dw = -opt.ww:opt.ws:opt.ww
        dh_s = 1+(opt.wh + dh)/opt.ws;
        dw_s = 1+(opt.ww + dw)/opt.ws;
        img2 = image_2(1+opt.wh+dh:ht1+opt.wh+dh,1+opt.ww+dw:wt1+opt.ww+dw,:);
        img2b = reshape(img2,[],bt2);
        [t2_feature_temp] = PrcfeatureExtraction(Nf,idxsup_fine,img2b) ;
        t2_feature{dh_s,dw_s} = t2_feature_temp;
    end
end
fprintf('\n');fprintf('The computational time of calculating feature templates t2 is %i \n', etime(clock, t_template));
%% Iterative framework
labels = zeros(opt.Niter,Nf);
iter = 1;
idex_unchange = 1:Nc;
CM = zeros(1,ht1*wt1);
library_changeP = zeros(Nc,1);
while iter<=(opt.Niter)
    %---------------------  template update ----------------------%
    t_di = clock;
    t1_feature_library = t1feature_coarse(idex_unchange,:);
    t2_feature_library = t2feature_coarse(idex_unchange,:);
    %---------------------  f calculation----------------------%
    
    [idsort_t1, dist_t1] = knnsearch(t1_feature_library,t1feature_fine,'k',Nf);
    for dh = -opt.wh:opt.ws:opt.wh
        for dw = -opt.ww:opt.ws:opt.ww
            dh_s = 1+(opt.wh + dh)/opt.ws;
            dw_s = 1+(opt.ww + dw)/opt.ws;
            t2_feature_temp = t2_feature{dh_s,dw_s};
            [idsort_t2_temp, dist_t2_temp] = knnsearch(t2_feature_library,t2_feature_temp,'k',opt.kmax);
            idsort_t2{dh_s,dw_s} = idsort_t2_temp;
            dist_t2{dh_s,dw_s} = dist_t2_temp;
            for i = 1: Nf
                id1 = idsort_t1(i,2:opt.kmax);
                distx = dist_t1(i,2:opt.kmax);
                id2 = idsort_t2_temp(i,2:opt.kmax);
                disty = dist_t2_temp(i,2:opt.kmax);
                dist_xy = pdist2(t1_feature_library(id2,:),t1feature_fine(i,:));
                dist_yx = pdist2(t2_feature_library(id1,:),t2_feature_temp(i,:));
                fx = (mean(dist_xy)-mean(distx))/bt1;
                fy = (mean(dist_yx)-mean(disty))/bt2;
                f(dh_s,dw_s,i) =  fx + fy ;
            end
        end
    end
    fprintf('\n');fprintf('The computational time of DI is %i \n', etime(clock, t_di));
    
    %---------------------  Dre calculation----------------------%
    t_cm = clock;
    Dre = zeros(Nf,2);
    fmin = zeros(Nf,1);
    for i = 1:Nf
        fv = f(:,:,i);
        [minValue, linearIndex] = min(fv(:));
        fmin(i) = minValue;
        [Dre(i,1), Dre(i,2)] = ind2sub(size(fv), linearIndex);
    end
    
    Dre(:,1) = (Dre(:,1)-1)*opt.ws -opt.wh;
    Dre(:,2) = (Dre(:,2)-1)*opt.ws -opt.ww;
    
    %---------------------  labels calculation----------------------%
    [labels(iter,:)] = solveEnergyMRF(fmin,Dre,Cosup_fine,opt);
    
    %---------------------  template update----------------------%
    img2_new = zeros(ht1*wt1,bt2);
    for i = 1:Nf
        index_vector = idxsup_fine{i};
        [h1, w1] = ind2sub([ht1 wt1],index_vector);
        h1new = h1 + opt.wh + Dre(i,1);
        w1new = w1 + opt.ww + Dre(i,2);
        idx_new = sub2ind([ht2 wt2],h1new,w1new);
        img2_new(index_vector,1:bt2) = imgt2b(idx_new,1:bt2);
        CM(index_vector) = labels(iter,i);
    end
    [t2feature_coarse] = PrcfeatureExtraction(Nc,idxsup_coarse,img2_new) ;
    for i = 1:Nc
        index_vector = idxsup_coarse{i};
        library_changeP(i) = mean(CM(index_vector));
    end
    idex_unchange = library_changeP<0.5;
    Dreiter(iter,:,:) = Dre;
    fminiter(iter,:) = (fmin -min(fmin)) /(max(fmin)-min(fmin));
    iter = iter+1;
    fprintf('\n');fprintf('The computational time of CM is %i \n', etime(clock, t_cm));
    
end

CM_map=zeros(ht1*wt1,opt.Niter);
DI_map=zeros(ht1*wt1,opt.Niter);
for i = 1:length(idxsup_fine)
    index_vector = idxsup_fine{i};
    for iter = 1:opt.Niter
        CM_map(index_vector,iter) = labels(iter,i);
        DI_map(index_vector,iter) = fminiter(iter,i);
    end
end
CM_map_iter = zeros(ht1+opt.wh*2,wt1+opt.ww*2,opt.Niter);
DI_map_iter = zeros(ht1+opt.wh*2,wt1+opt.ww*2,opt.Niter);
for iter = 1:opt.Niter
    CM_map_iter(opt.wh+1:opt.wh+ht1,opt.ww+1:opt.ww+wt1,iter)  =reshape(CM_map(:,iter),[ht1 wt1]);
    DI_map_iter(opt.wh+1:opt.wh+ht1,opt.ww+1:opt.ww+wt1,iter)  =reshape(DI_map(:,iter),[ht1 wt1]);
end
