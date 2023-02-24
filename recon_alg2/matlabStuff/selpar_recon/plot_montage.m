function plot_montage(path, Nimg, saveFlag)
% plot_montage(path, Nimg, saveFlag)
%
% Script to create a montage of SoS (sum of squares) selpar reconstructed
% images (all slices)
% 
% path = the relative path to SoS selpar recon file "sosimg.mat" 
%
% Nimg = the number of images (slices) per row
%
% saveFlag = save figure flag:
% 'no' (none, default)
% 'yes' (save)
%

close all, clc

if nargin == 2
   saveFlag = 'no'; 
end

load(strcat(path, '/sosimg.mat'))
fig_handle = figure('Position', get(0, 'Screensize'));

% set path to bias field maps
load('.../data/BiasMatrix.mat')
bias = bias./max(bias(:));
img = img./bias;
limit = 3.0;
img = img./mean(img(:));

montage(img, 'Size', [nan, Nimg], 'DisplayRange', [0 limit]), colorbar
set(gca, 'fontSize', 20)

if strcmp(saveFlag, 'yes') == 1
    saveas(fig_handle, strcat(path, '/montage.jpg'),'jpg');
end