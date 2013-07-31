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
%
starbucksSet.getLabels();

roadsignSet.getLabels();
%}

%% Seif 5
%

res = starbucksSet.getLabels(9);
fprintf('\n');

label_sb9 = res{1,9};
tform_sb9 = res{2,9};
sb9 = starbucksSet.cImages{9};

tform = maketform('projective',tform_sb9.T);

rs_tmpl = roadsignSet.cTemplate;
rs_tmpl_xformed = ...
    imtransform(rs_tmpl, tform, 'XData', [1 size(sb9,2)], 'YData' , [1 size(sb9,1)]);
bg_ind = (rs_tmpl_xformed(:,:,1) >= 240) & (rs_tmpl_xformed(:,:,2) >= 240) & (rs_tmpl_xformed(:,:,3) >= 240);
bg_ind = repmat(bg_ind, [1 1 3]);
rs_tmpl_xformed(bg_ind) = 0;

gray_ind = (rs_tmpl_xformed(:,:,1) == rs_tmpl_xformed(:,:,2)) & (rs_tmpl_xformed(:,:,2) == rs_tmpl_xformed(:,:,3));
gray_ind = repmat(gray_ind, [1 1 3]);
rs_tmpl_xformed(gray_ind) = 0;

new_im = sb9;
label_sb9 = logical(repmat(label_sb9, [1 1 3]));
new_im(label_sb9) = 0;

figure;
imshowpair(new_im, rs_tmpl_xformed, 'blend');
