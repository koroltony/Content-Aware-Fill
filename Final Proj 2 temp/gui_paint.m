function gui_paint
global patch_size;

if isfile('hole.png')
    delete('hole.png');
end

fig = figure('Name', 'Create Alpha Hole', 'Position', [100, 100, 800, 600]);

% Allows me to choose the image from my files
[file, path] = uigetfile({'*.jpg;*.jpeg;*.png;*.bmp', 'Image Files'});
img = imread(fullfile(path, file));
dim = size(img);
h = dim(1);
w = dim(2);
if h>400 || w>400
    r = h/w;
    img = imresize(img,[400*r 400]);
end

% Make the axes and show the image in the gui viewer
ax = axes('Parent', fig, 'Position', [0.05, 0.1, 0.7, 0.8]);
imshow(img);
title('Original Image');

% make a slider for paint width selection. Whenever it is changed, the
% callback runs the functions below to update values and logical arrays
uicontrol('Style', 'slider', 'Min', 1, 'Max', 200, 'Value', 10, 'Position', [700, 200, 100, 20], 'Callback', @brushSizeCallback);
PaintSize = uicontrol('Style', 'text', 'String', 'Brush Size: 10', 'Position', [700, 230, 100, 20]);

% Same idea for undo and done button
uicontrol('Style', 'pushbutton', 'String', 'Undo', 'Position', [700, 300, 100, 30], 'Callback', @undoCallback);
uicontrol('Style', 'pushbutton', 'String', 'Done', 'Position', [700, 350, 100, 30], 'Callback', @doneCallback);

uicontrol('Style', 'slider', 'Min', 3, 'Max', 15, 'Value', 3, 'Position', [700, 160, 100, 20], 'SliderStep',[1/6,1/6],'Callback', @patchSizeCallback);
patch_s = uicontrol('Style', 'text', 'String', 'Patch Size: 3', 'Position', [700, 180, 100, 20]);

% Initialize brush size at 10, set all actions to false
brushSize = 10;
patch_size = 3;
isPainting = false;

% this is a logical array to show whether the pixels had been painted
% on

paintedImg = false(size(img, 1), size(img, 2));



% update brush size whenever the size slider is changed
    function brushSizeCallback(size, ~)
        % set the actual brush size
        brushSize = round(get(size, 'Value'));
        % set the PaintSize text
        set(PaintSize, 'String', ['Brush Size: ', num2str(brushSize)]);
    end

% update patch size whenever the size slider is changed
    function patchSizeCallback(size, ~)
        % set the actual brush size
        patch_size = round(get(size, 'Value'));
        % set the PaintSize text
        set(patch_s, 'String', ['Patch Size: ', num2str(patch_size)]);
    end

% Updates the logical values which show painted pixels when undo button
% is clicked
    function undoCallback(~, ~)
        % sets all logical values to 0, so it's more like clear, not undo
        paintedImg = false(size(img, 1), size(img, 2));
        updateImage();
    end

% Saves the picture with alpha channel of 0 wherever the user painted
    function doneCallback(~, ~)

        imgSize = size(img);
        height = imgSize(1);
        width = imgSize(2);

        updatedImg = img;
        updatedImgr = updatedImg(:,:,1);
        updatedImgg = updatedImg(:,:,2);
        updatedImgb = updatedImg(:,:,3);
        
        updatedImgr(paintedImg) = 0;
        updatedImgg(paintedImg) = 0;
        updatedImgb(paintedImg) = 0;

        updatedImg(:,:,1) = updatedImgr;
        updatedImg(:,:,2) = updatedImgg;
        updatedImg(:,:,3) = updatedImgb;

        % make an alpha channel matrix with size of image
        alphaChannel = ones(height, width);
        alphaChannel(paintedImg) = 0;

        % Save the image with alpha channel
        imwrite(updatedImg, 'hole.png', 'Alpha', alphaChannel);
        disp('Image saved successfully');
    end



% Function to update the displayed image any time the mouse is clicked
% or any buttons are clicked. This will use the logical matrix to show
% red wherever the user paints
    function updateImage()
        updatedImg = img;

        % store a value of 255 for all pixels that were painted on
        redChannel = uint8(paintedImg) * 255;
        % everywhere that wasn't painted on is multiplied by 1 and the
        % painted part is multiplied by 0
        % Then, the red patch is added to the red channel of the image
        updatedImg(:,:,1) = updatedImg(:,:,1) .* uint8(~paintedImg) + redChannel;
        updatedImg(:,:,2) = updatedImg(:,:,2) .* uint8(~paintedImg);
        updatedImg(:,:,3) = updatedImg(:,:,3) .* uint8(~paintedImg);
        imshow(updatedImg, 'Parent', ax);
    end

% when mouse is clicked, toggle the painting enable
    function MouseDown(~, ~)
        isPainting = true;
    end

% when mouse is released don't paint anymore
    function MouseUp(~, ~)
        isPainting = false;
    end

% track the location of the user's mouse and determine the pixels
% affected by the brush
    function MouseMove(~, ~)
        % if the mouse is pressed
        if isPainting
            % finds the current position of the mouse using ax.CurrentPoint
            currentPoint = round(ax.CurrentPoint(1, 1:2));

            % Get brush size bounds by adding the brush size around the
            % current pixel location in the image. Make sure that the
            % drawing is in bounds (Kinda similar to the patch match
            % checks!
            minX = max(1, currentPoint(1) - brushSize);
            minY = max(1, currentPoint(2) - brushSize);
            maxX = min(size(img, 2), currentPoint(1) + brushSize);
            maxY = min(size(img, 1), currentPoint(2) + brushSize);

            % set the logical matrix to be true wherever the brush is
            % present
            paintedImg(minY:maxY, minX:maxX) = true;

            % Update image
            updateImage();
        end
    end

% matlab WindowButtonDown checks mouse movement
% which then runs MouseDown
set(fig, 'WindowButtonDownFcn', @MouseDown);

% same idea for mouse up and mose move
set(fig, 'WindowButtonUpFcn', @MouseUp);
set(fig, 'WindowButtonMotionFcn', @MouseMove);
end

