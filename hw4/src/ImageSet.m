classdef ImageSet
    
    properties (SetAccess=private)
        % color (original)
        cTemplate;
        cImages;
        % grayscale
        gTemplate;
        gImages;
        % L*a*b
        labTemplate;
        labImages;
    end
    
    properties (Dependent = true)
        size;
    end
    
    methods
        function set = ImageSet( imfolder , scaleToHeight )
            %% IMAGESET Loads image set with template image.
            %   Loads the contents of a directory containing an image set.
            %   An image set directory should contain .jpg files, with one named
            %   template.jpg, and the rest named arbitrarily.
            %   The returned set is structured in the following way:
            %       o set.cTemplate - color template image for the set
            %       o set.cImages - cell array with images in color
            %       o set.gTemplate - grayscale template image for the set
            %       o set.gImages - cell array with images in grayscale
            %       o set.labTemplate - L*a*b template image for the set
            %       o set.labImages - cell array with images in L*a*b
            fprintf(1,'\nImageSet: Loading set from %s...', imfolder); tic;
            imageFiles = strsplit(ls(imfolder)); % all file names in folder
            imageFiles = imageFiles(1:end-1); % remove empty string
            
            % initialize sets
            numImages = length(imageFiles) - 1; % subtract one for the template image
            set.cImages = cell(1, numImages);
            set.gImages = cell(1, numImages);
            
            idx = 1;
            
            for i_ = 1:length(imageFiles)
                imFileName = imageFiles{i_};
                
                % build full path to the image file
                filepath = [imfolder '/' imFileName];
                fprintf(1,'\n\tReading image "%s"...', imFileName);
                
                % read image as color and scale it
                imColor = imread(filepath);
                imColor = set.resize(imColor, scaleToHeight);
                
                % convert image to grayscale, and scale it to [0,1]
                if (ndims(imColor) == 3)
                    imGray = double(rgb2gray(imColor));
                else
                    imGray = double(imColor);
                end
                imGray = imGray - min(imGray(:));
                imGray = imGray / max(imGray(:));
                
                % convert image to L*a*b
                imLab = applycform(imColor, makecform('srgb2lab'));
                
                [~, name, ~] = fileparts(filepath);
                if (strcmpi(name, 'template'))
                    set.cTemplate = imColor;
                    set.gTemplate = imGray;
                    set.labTemplate = imLab;
                else
                    set.cImages{1,idx} = imColor;
                    set.gImages{1,idx} = imGray;
                    set.labImages{1,idx} = imLab;
                    idx = idx + 1;
                end
            end
            fprintf(1,'\nImageSet: ...Done (%.3fs).\n', toc);
        end
        
        function value = get.size(obj)
            %% SIZE returns the number of images in the set.
            value = length(obj.gImages);
        end
        
        function labelSet = getLabels(set, varargin)
            %% GETLABELS Find and label the template image in the other images.
            
            parser = inputParser;
            parser.addOptional('idx', 1:set.size, @isvector);
            parser.addParamValue('RecogThresh', 10, @isscalar);
            parser.addParamValue('MatchThresh', 0.9, @isscalar);
            parser.addParamValue('Descriptor', 'sift', @isstr);
            parser.addParamValue('TransformType', 'homography', @isstr);
            parser.parse(varargin{:});
            
            % disable some warnings
            warn_state = warning;
            warning('off','stats:kmeans:EmptyClusterRep');
            warning('off','stats:kmeans:EmptyCluster');
            warning('off','MATLAB:singularMatrix');
            warning('off','MATLAB:nearlySingularMatrix');
            
            labelSet = cell(2, set.size);
            
            % Compute homogeneous xy-coordinates of the corner pixels in the template
            [nrow, ncol,~] = size(set.cTemplate);
            templateHCoords = [0 0; ncol 0; ncol nrow; 0 nrow]; % flip columns/rows because it's an image...
            
            % Loop over all images, find object and apply segmentation
            templateKeypoints = [];     % template keypoints and descriptors will be cached
            templateDescriptors = [];   % to prevent re-computation each iteration.
            for i=parser.Results.idx
                
                fprintf(1,'\nImageSet:getLabels Processing image #%d...', i); tic;
                
                % Find some points in the test image that could belong to the object
                % and a transformation from template image coordinate to test image coordinates.
                [objectPoints, tform, templateKeypoints, templateDescriptors] = ...
                    findTemplate(set.cTemplate, set.cImages{i}, ...
                                 templateKeypoints, templateDescriptors, ...
                                 'transformtype', parser.Results.TransformType, ...
                                 'matchthresh', parser.Results.MatchThresh, ...
                                 'descriptor', parser.Results.Descriptor);
                
                if (size(objectPoints,2) < parser.Results.RecogThresh)
                    fprintf(2,' object not found, skipping.');
                    labelSet{1,i} = []; labelSet{2,i} = [];
                    continue;
                else
                    fprintf(1,' object found (%.3fs). Applying segmentation...', toc); tic;
                end
                
                % Compute the corner's xy-coordinates in the test image
                % of the region that might belong to the object.
                objectBounds = tform.transformPointsForward(templateHCoords)';
                
                % Create an approximate labeling - all pixels in the object
                % region will be labeled '1', the rest '0'.
                approximateLabels = ...
                    poly2mask(objectBounds(1,:), objectBounds(2,:), ...
                              size(set.gImages{i},1), size(set.gImages{i},2));
                approximateLabels = uint8(approximateLabels);
                
                % Now label all the keypoints that are on the object with '2'.
                objectInd = sub2ind(size(approximateLabels), round(objectPoints(2,:)), round(objectPoints(1,:)));
                approximateLabels(objectInd) = 2;

                
                % Use GraphCut to compute the exact labels from the approximate labels
                objectLabel = labelTemplate(set.cTemplate, set.cImages{i}, tform, approximateLabels);
                labelSet{1,i} = objectLabel;
                labelSet{2,i} = tform;
                fprintf(1,' done (%.3fs).', toc);
                             
                % DEBUG
                if (nargout == 0)
                    figure;
                    subplot_tight(2,2,1);
                    subimage(set.cImages{i}); axis off;
                    patch(objectBounds(1,:), objectBounds(2,:),ones(1,3),'EdgeColor','r','FaceColor','none', 'LineWidth', 1.5);
                    
                    subplot_tight(2,2,2);
                    subimage(mat2gray(imdilate(approximateLabels, strel('disk',2)))); axis off;
                    
                    subplot_tight(2,2,[3 4]);
                    subimage(imfuse(set.cImages{i}, labelSet{1,i})); axis off;
                    
                    drawnow;
                end
            end
            
            % restore warnings
            warning(warn_state);
        end
    end
    
    methods (Static, Access=private)
        function imResized = resize(im, scaleToHeight)
            % resize the image so it's height is scaleToHeight px high.
            [height, ~] = size(im);
            
            scale = scaleToHeight / height;
            imResized = imresize(im, scale);
        end
    end
    
end