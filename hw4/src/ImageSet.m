classdef ImageSet
    %IMAGESET Loads image set with template image.
    %   Loads the contents of a directory containing an image set.
    %   An image set directory should contain .jpg files, with one named
    %   template.jpg, and the rest named arbitrarily.
    %   The returned set is structured in the following way:
    %       o set.cTemplate - color template image for the set
    %       o set.cImages - cell array with images in color
    %       o set.gTemplate - grayscale template image for the set
    %       o set.gImages - array with images in grayscale
    
    properties (SetAccess=private)
        cTemplate;
        cImages;
        gTemplate;
        gImages;
    end
    
    properties (Dependent = true)
        size;
    end
    
    methods
        function set = ImageSet( imfolder , scaleToHeight )
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
                
                % read image as grayscale, and scale it to [0,1]
                imGray = imreadbw(filepath);
                imGray = imGray - min(imGray(:));
                imGray = imGray / max(imGray(:));
                imGray = set.resize(imGray, scaleToHeight);
                
                % read image as color
                imColor = imread(filepath);
                imColor = set.resize(imColor, scaleToHeight);
                
                [~, name, ~] = fileparts(filepath);
                if (strcmpi(name, 'template'))
                    set.cTemplate = imColor;
                    set.gTemplate = imGray;
                else
                    set.cImages{1,idx} = imColor;
                    set.gImages{1,idx} = imGray;
                    idx = idx + 1;
                end
                
            end
        end
        
        function value = get.size(obj)
            value = length(obj.gImages);
        end
    end
    
    methods (Static, Access=private)
        function imResized = resize(im, scaleToHeight)
            % resize the image so it's height is 500px.
            [height, ~] = size(im);
            
            scale = scaleToHeight / height;
            imResized = imresize(im, scale);
        end
    end
    
end