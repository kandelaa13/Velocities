%% A script to segment cells from BF images

% Name of the file, change!
rootname = 'D:\Bendings\11112022\20221111_noY27\20221111_BF_10x\20221111_BF_10x_f0001.tif';

% Parameters of the segmentation algorithm, to be optimized
area_min=10000;
min_object_size=area_min;
min_hole_size=800;
treshold_finetune=0.01;

% This commented lines and the ones after the plot are meant to generate a
% video from a sequences of plots, if the plot section is inside a loop,

%writerObj = VideoWriter('myVideo_corrected.avi');
%writerObj.FrameRate = 10;
%open(writerObj);

Fluo_tiff_info = imfinfo(rootname);
K = length(Fluo_tiff_info);


% Matlab function imread accepts a second input. imread(image) will only
% extract the first slice if the image is a tif stack. imread(image,k) will
% read and extract the k-th slice of the stack
imBF1 = double(imread(rootname));
res_EGT = EGT_Segmentation(imBF1, min_object_size,min_hole_size,treshold_finetune);

% This function grabs the Black&White image and assigns an integer number
% to each simply connected region in it.
CC = bwconncomp(res_EGT);

figure(1);
% imshow plots the image, visboundaries draws the red outlines
imshow(imBF1,[]); axis equal; colormap gray;hold on;
visboundaries(res_EGT)
drawnow
hold off
pause(0.05)
title(['EGT'])
% To generate a video. These two lines should be inside the loop.
%frame = getframe(gcf);
%writeVideo(writerObj,frame);

% This line should be after the loop has ended
%close(writerObj);

%{
if CC.NumObjects ~=1
    s = regionprops(res_EGT,'Centroid');
    centroids = cat(1,s.Centroid);
    pos_nuc = [xtrack(i),ytrack(i)];
    % Nearest neighbor search to find which region is the real cell
    Idx = knnsearch(centroids, pos_nuc);
    % Loop through all simply connected regions of the mask, to
    % eliminate all but the good one
    for ii = 1:CC.NumObjects
        if ii == Idx
            continue
        end
        res_EGT(CC.PixelIdxList{ii}) = 0;
    end
end
%}