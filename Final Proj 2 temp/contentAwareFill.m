% Run this file to use the GUI

global patch_size;
patch_size = 3;
gui_paint;



% Start monitoring the creation of "hole.png" in the background
waitForFileCreation('hole.png');

% Now that the file exists, proceed to read and process it
[hole_im, ~, alpha] = imread('hole.png');
% isolate the alpha and rgb channels of the image
hole_im = im2double(hole_im);
alpha = im2double(alpha);
out = proj(hole_im,alpha);
imwrite(out, "output.png");

[R,G,B] = imsplit(hole_im);

R = regionfill(R, ~alpha);
G = regionfill(G, ~alpha);
B = regionfill(B, ~alpha);

% Add the channels back into the RGB image

interp_img = cat(3, R, G, B);

imwrite(interp_img,'interpResult.png');

function waitForFileCreation(filename)
while ~isfile(filename)
    % Wait until the file is created
    pause(0.1); % Adjust the duration of the pause as needed
end
end



