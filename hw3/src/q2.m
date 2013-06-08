%% Computer Vision, spring 2013
%  
%  Assignment 3
%  Aviv Rosenberg

%% NOTE: This code assumes that 'SIFT for MATLAB' exists
% SIFT for MATLAB: http://www.vlfeat.org/~vedaldi/code/sift.html
% o Should be in the MATLAB path
% o Should be compiled

%% Load all images as grayscale
% vars beginning with t_ are temporary.
close all;
clear;
t_basepath = fileparts(mfilename('fullpath'));
t_imfolder = strcat(t_basepath,'/../img/');
t_imageFiles = strsplit(ls(t_imfolder));
t_imageFiles = t_imageFiles(1:end-1);

for t_imFileName = t_imageFiles
    % build a variable name in matlab wrokspace.
    t_varname = strsplit(t_imFileName{:},'.');
    t_varname = t_varname{1};
    t_varname = strcat('im', genvarname(t_varname));
    
    % build full path to the image file
    t_filepath = strcat(t_imfolder, t_imFileName{:});
    
    % read image as grayscale, and scale it to [0,1]
    t_im = imreadbw(t_filepath);
    t_im = (t_im - min(t_im(:))) / max(t_im(:));
    
    % set image into a new variable named according to 'varname'.
    assignin('base', t_varname, t_im); 
end
clearvars t_*;

%% Apply SIFT
%
[frames1,descr1] = sift(imStopSign1);
[frames2,descr2] = sift(imStopSign2);
[frames3,descr3] = sift(imStopSign3);
[frames4,descr4] = sift(imStopSign4);
