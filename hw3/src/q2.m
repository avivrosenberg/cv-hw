%% Computer Vision, spring 2013
%  
%  Assignment 3
%  Aviv Rosenberg

%% NOTE: This code assumes that 'SIFT for MATLAB' exists
% SIFT for MATLAB: http://www.vlfeat.org/~vedaldi/code/sift.html
% o Should be in the MATLAB path
% o Should be compiled

%% Load all images as grayscale
% vars ending with '_' are temporary.
close all;
clear;
basepath_ = fileparts(mfilename('fullpath'));
imfolder_ = [basepath_ '/../img/'];
imageFiles_ = strsplit(ls(imfolder_));
imageFiles_ = imageFiles_(1:end-1);
imVarNames = {};

for imFileName_ = imageFiles_
    % build a variable name in matlab wrokspace.
    varname_ = strsplit(imFileName_{:},'.');
    varname_ = varname_{1};
    varname_ = ['im' genvarname(varname_)];
    imVarNames = [imVarNames {varname_}];
    
    % build full path to the image file
    filepath_ = [imfolder_ imFileName_{:}];
    
    % read image as grayscale, and scale it to [0,1]
    im_ = imreadbw(filepath_);
    im_ = (im_ - min(im_(:))) / max(im_(:));
    
    % set image into a new variable named according to 'varname'.
    assignin('base', varname_, im_); 
end
clearvars *_;

%% Apply SIFT
%
[frames1,descr1] = sift(imStopSign1);
[frames2,descr2] = sift(imStopSign2);
[frames3,descr3] = sift(imStopSign3);
[frames4,descr4] = sift(imStopSign4);
