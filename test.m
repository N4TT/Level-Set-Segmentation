clear all;
C = double(rgb2gray(imread('ct.png')));
mask = zeros(size(C));
mask(115:125,220:230) = 1;
imshow(mask);
simpleseg(C, mask, 1500, 30, 90, 0.5);