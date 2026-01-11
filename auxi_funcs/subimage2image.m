function DI = subimage2image(segDI,image,height,width,Niter)
[H, W, ~] = size(image);
hp = ceil(H/height);
wp = ceil(W/width);
DI = zeros(H,W,Niter);
for row = 1:hp-1
    for col = 1:wp-1
        DI((row-1)*height+1:row*height,(col-1)*width+1:col*width,:) = segDI{row,col};
    end
end
for col = 1:wp-1
    DI((hp-1)*height+1:H,(col-1)*width+1:col*width,:)=segDI{hp,col};
end
for row = 1:hp-1
    DI((row-1)*height+1:row*height,(wp-1)*width+1:W,:) = segDI{row,wp};
end
DI((hp-1)*height+1:H,(wp-1)*width+1:W,:) = segDI{hp,wp};