%% Computer Vision, spring 2013
%
%  Assignment 4
%  Aviv Rosenberg

%% Load Image Sets
%
close all;
clear;

basepath = fileparts(mfilename('fullpath'));
imfolder = [basepath '/../img/'];
imHeight = 500; % desired height in pixels to load the images in

supermanSet = ImageSet([imfolder 'Superman'], imHeight);
starbucksSet = ImageSet([imfolder 'Starbucks'], imHeight);
roadsignSet = ImageSet([imfolder 'Roadsign'], imHeight);