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
    im_ = (im_ - min(im_(:))) / max(im_(:));
    
    % set image into a new variable named according to 'varname'.
    assignin('base', varname_, im_);
end
clearvars *_ -except *__;

%% Apply SIFT
%
for i_ = 1:length(imVarNames__);
    imVarName_ = imVarNames__{i_};
    im_ = eval(imVarName_);
    
    fprintf('\nComputing frames and descriptors for "%s"... ', imVarName_); tic;
    [frames_, descr_] = sift(im_);
    fprintf('done (%.3fs)\n',toc);
    
    assignin('base', [framesPrefix__ num2str(i_)], frames_);
    assignin('base', [descriptorsPrefix__  num2str(i_)], descr_);
end
clearvars *_ -except *__;

%% q2
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
    ind_ = ind_(1:25);
    
    h_ = plotsiftframe(frames_(:,ind_)); set(h_,'LineWidth',1,'Color','g');
    h_ = plotsiftdescriptor(descr_(:,ind_), frames_(:,ind_));
end
clearvars *_ -except *__;

%% Clean up variables
%
clearvars *_;
