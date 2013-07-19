%% Computer Vision, spring 2013
%
%  Assignment 4
%  Aviv Rosenberg

%% Load Image Sets
%
close all;
clearvars -except *Set;

basepath_ = fileparts(mfilename('fullpath'));
imfolder_ = [basepath_ '/../img/'];
imHeight_ = 300; % desired height in pixels to load the images in

if (~exist('supermanSet', 'var'))
    supermanSet = ImageSet([imfolder_ 'Superman'], imHeight_);
end
if (~exist('roadsignSet', 'var'))
    roadsignSet = ImageSet([imfolder_ 'Roadsign'], imHeight_);
end
if (~exist('starbucksSet', 'var'))
    starbucksSet = ImageSet([imfolder_ 'Starbucks'], imHeight_);
end

clearvars *_;

