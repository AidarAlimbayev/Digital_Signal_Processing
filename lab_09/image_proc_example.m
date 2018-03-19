L = imread('aidar.jpg');
L = rgb2gray(L);
%% L = imrotate(L,35,'bicubic');
%%imshow(L);
%%figure, imhist(L);
%%L = histeq(L);
%%figure, imshow(L);
%H = fspecial('unsharp');
%L = imfilter(L,H,'replicate');
%h = fspecial('prewitt');
%L1 = imfilter(L, h, 'replicate');
%figure, imshow(L1)
%%imcrop;