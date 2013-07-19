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
    
        function labelSet = getLabels(set)
            %% GETLABELS Find and label the template image in the other images.
            
            labelSet = cell(1, set.size);
            templateKeypoints = [];
            templateDescriptors = [];
            for i=1:set.size
                
                fprintf(1,'\nImageSet:getLabels Processing image #%d...', i); tic;
                
                [objectPoints, tform, templateKeypoints, templateDescriptors] = ...
                    findTemplate(set.gTemplate, set.gImages{i}, ...
                    templateKeypoints, templateDescriptors, 'maxclusters',9,'matchthresh',0.9,'transformtype','affine');
                
                if (size(objectPoints,2) < 10)
                    fprintf(2,' object not found, skipping.');
                    labelSet{1,i} = [];
                    continue;
                else
                    fprintf(1,' object found (%.3fs). Applying segmentation...', toc); tic;
                end             
                
                backgroundRegion = [];
                
                objectLabel = labelTemplate(objectPoints, backgroundRegion);
                labelSet{1,i} = objectLabel;
                fprintf(1,' done (%.3fs).', toc);
            end

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