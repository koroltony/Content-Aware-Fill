% run this on a predetermined hole: setHolecows
% if you want to use the GUI, run contentAwareFill

global patch_size;
patch_size = 5;

% Now that the file exists, proceed to read and process it

[hole_im, ~, alpha] = imread('setHolecows.png');

% isolate the alpha and rgb channels of the image
hole_im = im2double(hole_im);
alpha = im2double(alpha);
out = proj(hole_im,alpha);
imwrite(out, "output1.png");

[R,G,B] = imsplit(hole_im);

R = regionfill(R, ~alpha);
G = regionfill(G, ~alpha);
B = regionfill(B, ~alpha);

% Add the channels back into the RGB image

interp_img = cat(3, R, G, B);

imwrite(interp_img,'interpResult.png');