RGB = reshape(ones(64,1)*reshape(jet(64),1,192),[64,64,3]);
R = RGB(:,:,1);
G = RGB(:,:,2);
B = RGB(:,:,3);
imshow(R)
figure, imshow(G)
figure, imshow(B)
figure, imshow(RGB)