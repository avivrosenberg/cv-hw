function [ out_im ] = kTextureCluster( im, varargin )
%KTEXTURECLUSTER Cluster image regions according to texture

%% Paramters Parsing
parser = inputParser;
parser.addRequired('im', @(x) ~isrow(x) && ~iscolumn(x) && ~isscalar(x));
parser.addParamValue('k',8, @isscalar);
parser.addParamValue('angles',[0 60 120 180 240 300], @(x) isrow(x) || iscolumn(x));
parser.addParamValue('sizes',[3 7 11], @(x) isrow(x) || iscolumn(x));
parser.addParamValue('sigma',1, @isscalar);
parser.addParamValue('rgb', false, @islogical);
parser.parse(im, varargin{:});

%% Apply filters to image
fb = gaussianFilterBank(parser.Results.angles, parser.Results.sizes, 'sigma', parser.Results.sigma);

if (size(im,3) > 1 && ~parser.Results.rgb)
    im = rgb2gray(im);
end

%% Nested function to create a clustered image
    function im_clustered = clusterImage(im)
        observations = zeros(numel(im), numel(fb));
        
        % Calculate responce to all filters
        for i=1:numel(fb)
            im_filtered = imfilter(im, fb{i});
            observations(:,i) = im_filtered(:);
        end
        
        % Cluster according to filter response
        idx = kmeans(observations, parser.Results.k,...
            'start','cluster',...
            'emptyaction','drop');
        
        cMax = max(idx); cMin = min(idx);
        % Use the indices themselves as pixel values for the output image.
        im_clustered = reshape(idx, size(im,1), size(im,2));
        % Scale the output.
        im_clustered = uint8(255 .* (im_clustered - cMin) ./ (cMax - cMin));
    end

%% Create output image: superimpose all cluster images
out_im = uint8(zeros(size(im)));
for j=1:size(im,3)
    out_im(:,:,j) = clusterImage(im(:,:,j));
end

end