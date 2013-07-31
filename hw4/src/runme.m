%% Computer Vision, spring 2013
%
%  Assignment 4
%  Aviv Rosenberg

%% Set up path
%
basepath_ = fileparts(mfilename('fullpath'));
imfolder_ = [basepath_ '/../img/'];
libsfolder_ = [basepath_ '/../../common/'];

addpath(basepath_);
addpath(genpath(libsfolder_));

%% Compile external libs
%
fprintf('\n\n');

cd([libsfolder_ 'GCmex2.0'])
fprintf(2,'Compiling GCmex2.0...\n'); tic
compile_gc;
fprintf(2,'\n ... done. (%.3f s)\n\n', toc);

cd([libsfolder_ 'sift'])
fprintf(2,'Compiling SIFT...\n'); tic
sift_compile;
fprintf(2,'\n ... done. (%.3f s)\n\n', toc);

cd([libsfolder_ 'ssdesc-cpp-1.1.1'])
fprintf(2,'Compiling ssdesc-cpp-1.1.1...\n'); tic
mex mexCalcSsdescs.cc ssdesc.cc;
fprintf(2,'\n ... done. (%.3f s)\n\n', toc);

cd([basepath_ '/../../']);
%}

%% Load Image Sets
%
close all;
clearvars -except *Set *_;

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

%% Example run

starbucksSet.getLabels();

roadsignSet.getLabels();
