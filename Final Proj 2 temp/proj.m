function [out] = proj(hole_im, alpha)
% Find the maximum dimension of the initial hole
maxDim = holeMaxDim(alpha);

global patch_size;
pad_size = (patch_size - 1) / 2;

% Number of scaled images you want
numScales = 10;

% First index is the RGB, second is the alpha
scales = cell(numScales, 2);

% How much we need to scale down to get to the coarsest scale
total_scale = maxDim / patch_size;

% How much we need to scale down in each step
step_scale = total_scale / (numScales - 1);

% Now iterate until the hole size is smaller than or equal to the patch size
scales{1,1} = hole_im;
scales{1,2} = alpha;

for i = 2:numScales
    % Pad the hole image and alpha mask before resizing
    padded_hole_im = padarray(hole_im, [pad_size, pad_size], 'symmetric');
    padded_alpha = padarray(alpha, [pad_size, pad_size], 'symmetric');

    % Resize the padded images
    scaled_hole_im = imresize(padded_hole_im, 1 / (step_scale * (i - 1)), 'nearest');
    scaled_alpha = imresize(padded_alpha, 1 / (step_scale * (i - 1)), 'nearest');

    % Store the scaled images
    scales{i,1} = scaled_hole_im;
    scales{i,2} = scaled_alpha;
end

% Now do the search and vote
% Start at the coarsest scale and interpolate the pixels in the
% hole region before running patch match
% Start by making a mask which includes the hole pixels (pixels with alpha
% values less than 1)
target_mask = scales{numScales,2} < 1;

% Create a ones matrix of size patch_size
ones_patch = ones(patch_size);

% Perform 2D convolution with the target_mask
overlap_region = conv2(target_mask, ones_patch, 'same');

% Check if the overlap_region is greater than zero to denote the region where patches overlap with the hole
extended_target_mask = overlap_region > 0;
extended_source_mask = ~(extended_target_mask);


% Interpolate using regionfill for the coarse scale hole
% To interpolate, we need to first make grayscale versions of each channel
[R,G,B] = imsplit(scales{numScales,1});
R = regionfill(R, target_mask);
G = regionfill(G, target_mask);
B = regionfill(B, target_mask);

% Add the channels back into the RGB image
interp_img = cat(3, R, G, B);
%figure;
%imshow(imresize(interp_img, [400, 400], 'nearest'));

% Now that we have the interpolated coarse scale image, we can start the
% patch-match
% The target image can just be the interpolated image
for k = 1:5
    coarse_NNF = patchMatchNNF1(interp_img, scales{numScales,1}, extended_target_mask, extended_source_mask);
    out = voteNNF(coarse_NNF, scales{numScales,1},interp_img);
    interp_img = out;
    %figure;
    %imshow(imresize(interp_img, [400, 400], 'nearest'));
end
%figure;
%imshow(imresize(interp_img, [400, 400], 'nearest'));

% Now upsample and blend with the hole region of the next image
for i = numScales-1:-1:1
    disp('Iteration: ');
    disp(numScales-i+1);
    % Upsample previous smaller filled image to the size of the current iteration
    upsample = imresize(out, [size(scales{i,1}, 1), size(scales{i,1}, 2)], 'nearest');

    % Extract the target mask of the current unfilled image
    new_target_mask = scales{i,2} < 1;
    % Create a ones matrix of size patch_size
    ones_patch = ones(patch_size);

    % Perform 2D convolution with the target_mask
    overlap_region = conv2(new_target_mask, ones_patch, 'same');

    % Check if the overlap_region is greater than zero to denote the region where patches overlap with the hole
    extended_target_mask = overlap_region > 0;
    extended_source_mask = ~(extended_target_mask);


    % Replace the regions with alpha less than one in the current unfilled image
    % with the corresponding regions from the upsampled image
    new_target = scales{i,1};
    % figure;
    % imshow(imresize(upsample, [400, 400], 'nearest'));
    upsample = new_target_mask .* upsample;
    % figure;
    % imshow(imresize(new_target_mask, [400, 400], 'nearest'));
    % figure;
    % imshow(imresize(upsample, [400, 400], 'nearest'));
    new_target = (~new_target_mask) .* new_target;
    % figure;
    % imshow(imresize(new_target, [400, 400], 'nearest'));
    new_target = new_target + upsample;
    % figure;
    % imshow(imresize(new_target, [400, 400], 'nearest'));

    % Run NNF and vote on the new image set
    for k = 1:5
        new_NNF = patchMatchNNF1(new_target, scales{i,1}, extended_target_mask, extended_source_mask);
        out = voteNNF(new_NNF, scales{i,1},new_target);
        new_target = out;
    end
end

% Output the final result
out = new_target;
end