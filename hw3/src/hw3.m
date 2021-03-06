%% Computer Vision, spring 2013
%
%  Assignment 3
%  Aviv Rosenberg

%% NOTE: This code assumes that 'SIFT for MATLAB' exists
% SIFT for MATLAB: http://www.vlfeat.org/~vedaldi/code/sift.html
% o Should be in the MATLAB path
% o Should be compiled

close all;
clear;

% NOTE: vars ending with '_' are temporary.
framesPrefix__ = 'frames';
descriptorsPrefix__ = 'descr';
matchesPrefix__ = 'matches';

%% Load all images as grayscale
%
basepath_ = fileparts(mfilename('fullpath'));
imfolder_ = [basepath_ '/../img/'];
imageFiles_ = strsplit(ls(imfolder_));
imageFiles_ = imageFiles_(1:end-1);
imVarNames__ = cell(1, length(imageFiles_));

for i_ = 1:length(imageFiles_)
    imFileName_ = imageFiles_{i_};
    % build a variable name in matlab wrokspace.
    varname_ = strsplit(imFileName_,'.');
    varname_ = varname_{1};
    varname_ = ['im' genvarname(varname_)];
    imVarNames__(i_) = {varname_};
    
    % build full path to the image file
    filepath_ = [imfolder_ imFileName_];
    
    % read image as grayscale, and scale it to [0,1]
    im_ = imreadbw(filepath_);
    im_ = im_ - min(im_(:));
    im_ = im_ / max(im_(:));
    
    % set image into a new variable named according to 'varname'.
    assignin('base', varname_, im_);
    
    % also read as color
    im_ = imread(filepath_);
    assignin('base', [varname_ 'c'], im_);
    
end
clearvars *_ -except *__;

%% Apply SIFT
%
for i_ = 1:length(imVarNames__);
    imVarName_ = imVarNames__{i_};
    im_ = eval(imVarName_);
    
    fprintf(2,'\nComputing frames and descriptors for "%s"... ', imVarName_); tic;
    %[frames_, descr_] = sift(im_, 'Threshold', 0.075, 'FirstOctave', 0);
    [frames_, descr_] = sift(im_, 'FirstOctave', 0);
    fprintf(2,'done (%.3fs)\n',toc);
    
    assignin('base', [framesPrefix__ num2str(i_)], frames_);
    assignin('base', [descriptorsPrefix__  num2str(i_)], descr_);
end
clearvars *_ -except *__;

%% q2: Finding keypoint frames and descriptors
%
for i_ = 1:length(imVarNames__);
    imVarName_ = imVarNames__{i_};
    im_ = eval(imVarName_);
    
    figure(i_); clf;
    imagesc(im_); colormap gray;
    hold on;
    
    frames_ = eval([framesPrefix__ num2str(i_)]);
    descr_ = eval([descriptorsPrefix__ num2str(i_)]);
    
    % randomly take some of the frames
    ind_ = randperm(size(frames_,2));
    ind_ = ind_(1:10);
    
    h_ = plotsiftframe(frames_(:,ind_)); set(h_,'LineWidth',1,'Color','g');
    h_ = plotsiftdescriptor(descr_(:,ind_), frames_(:,ind_));
end
clearvars *_ -except *__;
%
%% q3 & q4: Finding and plotting matches
%
im1_ = eval(imVarNames__{1});
descr1_ = eval([descriptorsPrefix__ '1']);
frames1_ = eval([framesPrefix__ '1']);
for i_ = 2:length(imVarNames__);
    % Obtain the two descriptor matrices to match
    descr2_ = eval([descriptorsPrefix__ num2str(i_)]);
    
    % Match them
    matches_ = findMatches(descr1_, descr2_, 'thresh', 1.0);
    
    % Assign to a variable in the workspace
    assignin('base', [matchesPrefix__ '1' num2str(i_)], matches_);
    
    % plot matches
    im2_ = eval(imVarNames__{i_});
    frames2_ = eval([framesPrefix__ num2str(i_)]);
    figure(str2double(['11' num2str(i_)])); axis image; clf;
    plotMatches(im1_,im2_,frames1_(1:2,:),frames2_(1:2,:),matches_,'Placement','horz');
end
clearvars *_ -except *__;
%
%% q5: Find affine matches
%
im1_ = eval(imVarNames__{1});
frames1_ = eval([framesPrefix__ '1']);
for i_ = 2:length(imVarNames__);
    % Obtain SIFT frame of image i_
    frames2_ = eval([framesPrefix__ num2str(i_)]);
    
    % Obtain the matches we found previously
    matches_ = eval([matchesPrefix__ '1' num2str(i_)]);
    
    % extract the image coordinates of the matching points
    pts1_ = frames1_(1:2,matches_(1,:));
    pts2_ = frames2_(1:2,matches_(2,:));
    
    % find affine matches
    [H_, inliers_, ~] = affineMatches(pts1_, pts2_);
    
    % extract inliers
    matches_ = matches_(:, inliers_);
    
    % Assign new matches to a variable in the workspace
    assignin('base', [matchesPrefix__ 'Affine1' num2str(i_)], matches_);
    assignin('base', ['HA1' num2str(i_)], H_);
    
    % plot new matches
    im2_ = eval(imVarNames__{i_});
    figure(str2double(['21' num2str(i_)])); axis image; clf;
    plotMatches(im1_,im2_,frames1_(1:2,:),frames2_(1:2,:),matches_,'Placement','horz');
end
clearvars *_ -except *__;
%
%% q6: Align Affine
%
% Obtain first color image
im1_ = eval([imVarNames__{1} 'c']);
for i_ = 2:length(imVarNames__);
    % Obtain second color image
    im2_ = eval([imVarNames__{i_} 'c']);
    
    % Obtain affine transform between them
    H_ = eval(['HA1' num2str(i_)]);
    
    % Compute aligned image
    imAligned_ = alignImages(im1_, im2_, H_);
    
    % Show it
    figure(str2double(['31' num2str(i_)]));
    imshow(imAligned_);
end
clearvars *_ -except *__;
%
%% q8: Find homography matches
%
im1_ = eval(imVarNames__{1});
im1c_ = eval([imVarNames__{1} 'c']);
frames1_ = eval([framesPrefix__ '1']);
for i_ = 2:length(imVarNames__);
    % Obtain SIFT frame of image i_
    frames2_ = eval([framesPrefix__ num2str(i_)]);
    
    % Obtain the matches we found previously
    matches_ = eval([matchesPrefix__ '1' num2str(i_)]);
    
    % extract the image coordinates of the matching points
    pts1_ = frames1_(1:2,matches_(1,:));
    pts2_ = frames2_(1:2,matches_(2,:));
    
    % find homography matches
    [H_, inliers_, ~] = homographyMatches(pts1_, pts2_);
    
    % extract inliers
    matches_ = matches_(:, inliers_);
    
    % Assign new matches to a variable in the workspace
    assignin('base', [matchesPrefix__ 'Homography1' num2str(i_)], matches_);
    assignin('base', ['HH1' num2str(i_)], H_);
    
    % plot new matches
    im2_ = eval(imVarNames__{i_});
    figure(str2double(['41' num2str(i_)])); axis image; clf;
    plotMatches(im1_,im2_,frames1_(1:2,:),frames2_(1:2,:),matches_,'Placement','horz');
    
    % Align images
    % Obtain second color image
    im2c_ = eval([imVarNames__{i_} 'c']);
    
    % Compute aligned image
    imAligned_ = alignImages(im1c_, im2c_, H_);
    
    % Show it
    figure(str2double(['51' num2str(i_)]));
    imshow(imAligned_);
end
clearvars *_ -except *__;
%% Clean up variables
%
clearvars *_;
