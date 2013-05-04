function [ ] = hsvHist( im, varargin)
%HSVHIST Show histogram of Hue before and after quantization.

%% Paramters Parsing
parser = inputParser;
parser.addRequired('im', @(x) ~isrow(x) && ~iscolumn(x) && ~isscalar(x) && size(x,3) == 3);
parser.addParamValue('k',16, @isscalar);
parser.parse(im, varargin{:});

%% Calculate original and quantized HSV values
[imSz1, imSz2, imSz3] = size(im);

im_hsv = rgb2hsv(im);
h_orig = reshape(im_hsv(:,:,1), imSz1 * imSz2, 1);

[im_q, levels] = kHSVquantizer(im, 'k', parser.Results.k);

%% Display histograms for original image
figure;
[binCounts, binCenters] = hist(h_orig, length(levels));
bar(binCenters, binCounts,'hist');
axis([0 1 0 max(binCounts)]);
title(['Original image: Hue histogram, ' num2str(length(levels)) ' uniform bins'], 'FontSize', 13);

figure; hist(h_orig, sort(levels));
axis([0 1 0 max(binCounts)]);
title(['Original image: Hue histogram, ' num2str(length(levels)) ' bins @ quantization levels'], 'FontSize', 13);

%% Display histograms for quantized image
im_hsv_q = rgb2hsv(im_q);
h_q = reshape(im_hsv_q(:,:,1), imSz1 * imSz2, 1);

figure;
[binCounts, binCenters] = hist(h_q, length(levels));
bar(binCenters, binCounts,'hist');
axis([0 1 0 max(binCounts)]);
title(['Quantized image: Hue histogram, ' num2str(length(levels)) ' uniform bins'], 'FontSize', 13);

figure; hist(h_q, sort(levels));
axis([0 1 0 max(binCounts)]);
title(['Quantized image: Hue histogram, ' num2str(length(levels)) ' bins @ quantization levels'], 'FontSize', 13);

%% Report the error
fprintf('SSD Error between original nd quantized image: %d\n\n',rgbSSD(im,im_q));
end

