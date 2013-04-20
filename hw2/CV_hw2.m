%% Computer Vision, spring 2013
%  
%  Assignment 2
%  Aviv Rosenberg

%%
% Load all images
close all;
clear;

im_colorful1 = imread('./img/colorful1.jpg');
im_colorful2 = imread('./img/colorful2.jpg');
im_colorful3 = imread('./img/colorful3.png');
im_ladybug = imread('./img/ladybug.jpg');
im_mooncraters = imread('./img/MoonCraters.jpg');
im_planets = imread('./img/Planets.jpeg');

%% Question 1
% 
% 
detectCircles(im_mooncraters);
detectCircles(im_ladybug);
detectCircles(im_colorful3);
detectCircles(im_planets);
