function seg = image2subimage(image,height,width)
[H W B] = size(image);
hp = ceil(H/height);
wp = ceil(W/width);
seg = cell(hp,wp);

for row = 1:hp-1
    for col = 1:wp-1
        seg(row,col)= {image((row-1)*height+1:row*height,(col-1)*width+1:col*width,:)};
    end
end
for col = 1:wp-1
    seg(hp,col)= {image((hp-1)*height+1:H,(col-1)*width+1:col*width,:)};
end
for row = 1:hp-1
    seg(row,wp)= {image((row-1)*height+1:row*height,(wp-1)*width+1:W,:)};
end
seg(hp,wp)= {image((hp-1)*height+1:H,(wp-1)*width+1:W,:)};
