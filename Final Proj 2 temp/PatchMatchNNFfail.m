% computes the NNF between patches in the target image and those in the source image
function NNF = PatchMatchNNFfail(target_image, source_image,target_mask,valid_source)
global patch_size;

fprintf("Computing NNF using PatchMatch...\n");

target_size = size(target_image);
source_size = size(source_image);

% initialize the NNF
NNF = zeros(target_size(1), target_size(2), 2);

% set number of iterations:

% seems to converge after 5 iterations
num_iter = 10;

% create masked target and source images:

padSize = (patch_size-1)/2;
pad_target_image = padarray(target_image,[padSize padSize],nan,"both");
pad_source_image = padarray(source_image,[padSize padSize],"symmetric");
pad_target_mask = padarray(target_mask,[padSize padSize],0);
imwrite(pad_target_mask,'targetmask.png');

tic

% --------------------- Step 1: Random Initialization -----------------

% for every pixel in the target image, we randomly assign a pixel in
% the source image, which basically means we create a random NNF to
% start with. Remember, our NNF shows offsets not coordinates.

% It's easier to first assign the coordinate and then get the offset

% NNF should map from 1 to size(source):

NNF(:,:,1) = randi([1 source_size(1)],[target_size(1) target_size(2)]);
NNF(:,:,2) = randi([1 source_size(2)],[target_size(1) target_size(2)]);

sourcerow = NNF(:,:,1);
sourcecol = NNF(:,:,2);

% extract the coordinate matrix of each target pixel by creating a
% meshgrid (it gives the row and column for each pixel to be used in
% calculating offset)

[targetcol,targetrow] = meshgrid(1:target_size(2), 1:target_size(1));

% now find the offsets: (source minus target)

NNF(:,:,1) = sourcerow - targetrow;
NNF(:,:,2) = sourcecol - targetcol;

% modified to make offsets outside of target 0

% create temporary NNF to store NNF during each iteration

NNF(:,:,1) = NNF(:,:,1).*target_mask;
NNF(:,:,2) = NNF(:,:,2).*target_mask;
newNNF = NNF;

imwrite(newNNF(:,:,1),'nnf.png');
imwrite(newNNF(:,:,2),'nnf2.png');

% --------------------- Step 2: Propagation ---------------------------

% determine indices for iteration:

for iter = 1:num_iter

    % check if we are propogating from top left or bottom right
    if mod(iter,2)==1
        % if we are on an odd iteration, we start in the top left
        i_seq = 1:target_size(1);
        j_seq = 1:target_size(2);
    else
        % on an even iteration start at the bottom right
        i_seq = target_size(1):(-1):1;
        j_seq = target_size(2):(-1):1;
    end

    % now actually iterate through the target and check matches

    for i = i_seq
        for j = j_seq
                if ~target_mask(i, j)
                    continue;
                end
            % extract patch we are looking at in target
            target_patch = pad_target_image(i:i+patch_size-1,j:j+patch_size-1,:);

            % find corresponding source patch using existing NNF
            source_row = i + NNF(i,j,1);
            source_col = j + NNF(i,j,2);
            source_patch = pad_source_image(source_row:source_row+patch_size-1,source_col:source_col+patch_size-1,:);

            % for odd iterations
            if mod(iter,2) == 1

                % find patch on top of current target patch
                % if we are in the first row of the target, just put the
                % original source patch

                target_top_row = max([1 i-1]);
                target_top_col = j;

                % find the shift for source patch using original NNF
                % shift_top = NNF(target_top_row,target_top_col,:);

                % make sure that the offsets used follow the NNF from the
                % first step. In row one it should be the same but in later
                % rows it should be the one above, so at i = 1 we just use
                % a row value of 1

                shift_top = NNF(max(1, i-1), j, :);

                % find source locations of top pixel
                source_row_top = target_top_row + shift_top(1);
                source_col_top = target_top_col + shift_top(2);
                source_patch_top = pad_source_image(source_row_top:source_row_top+patch_size-1,source_col_top:source_col_top+patch_size-1,:);

                % find the left pixel location in target and then find the
                % left pixel in the source

                target_left_row = i;
                target_left_col = max([1 j-1]);

                % find the shift for source patch using original NNF
                % shift_left = NNF(target_left_row,target_left_col,:);

                % make sure that the offsets used follow the NNF from the
                % first step. In row one it should be the same but in later
                % rows it should be the one to the left.
                shift_left = NNF(i, max(1, j-1), :);

                % find source locations of left pixel
                source_row_left = target_left_row + shift_left(1);
                source_col_left = target_left_col + shift_left(2);

                source_patch_left = pad_source_image(source_row_left:source_row_left+patch_size-1,source_col_left:source_col_left+patch_size-1,:);
                % compute the distances for all soure patches and target
                % patch to identify which one matches most closely
                currDist = patchDistance(source_patch,target_patch);
                topDist = patchDistance(source_patch_top,target_patch);
                leftDist = patchDistance(source_patch_left,target_patch);

                % compute the minimum distance between source patches
                [~,lowestInd] = min([currDist topDist leftDist]);

                % if the minimum distance is top or left, we update the NNF
                % to be close to that pixel
                switch lowestInd
                    case 1

                        newNNF(i,j,:) = NNF(i,j,:);

                    case 2
                        % store NNF of top pixel in best_NNF

                        best_NNF = shift_top;

                        % check if we can use the NNF below it to assign to
                        % the current pixel

                        if i + best_NNF(1) + 1 <= source_size(1)
                            best_NNF(1) = best_NNF(1) + 1;
                        elseif i + best_NNF(1) <= source_size(1)
                            best_NNF(1) = best_NNF(1);
                        else
                            best_NNF = NNF(i,j,:);
                        end

                        newNNF(i,j,:) = best_NNF;

                    case 3
                        % store NNF of left pixel in best_NNF

                        best_NNF = shift_left;

                        % check if we can use the NNF to the right of it to
                        % assign the current pixel

                        if j + best_NNF(2) + 1 <= source_size(2)
                            best_NNF(2) = best_NNF(2) + 1;
                        elseif j + best_NNF(2) <= source_size(2)
                            best_NNF(2) = best_NNF(2);
                        else
                            best_NNF = NNF(i,j,:);
                        end

                        newNNF(i,j,:) = best_NNF;
                end

                % if we are on an even iteration
            else

                % find patch on bottom of current target patch
                % if we are in the last row of the target, just put the
                % original source patch

                target_bot_row = min([target_size(1) i+1]);
                target_bot_col = j;

                % find the shift for source patch using original NNF
                % make sure that the offsets used follow the NNF from the
                % first step. In last row it should be the same as before

                shift_bot = NNF(min(1+i,target_size(1)),j,:);

                source_row_bot = target_bot_row + shift_bot(1);
                source_col_bot = target_bot_col + shift_bot(2);
                source_patch_bot = pad_source_image(source_row_bot:source_row_bot+patch_size-1,source_col_bot:source_col_bot+patch_size-1,:);

                % find the right pixel location in target and then find the
                % right pixel in the source

                target_right_row = i;
                target_right_col = min([target_size(2) j+1]);

                % make sure that the offsets used follow the NNF from the
                % first step. In last row it should be the same but in later
                % rows it should be the one to the right.
                shift_right = NNF(i,min(1+j,target_size(2)),:);

                % find source locations of right pixel
                source_row_right = target_right_row + shift_right(1);
                source_col_right = target_right_col + shift_right(2);
                source_patch_right = pad_source_image(source_row_right:source_row_right+patch_size-1,source_col_right:source_col_right+patch_size-1,:);
                % compute the distances for all soure patches and target
                % patch to identify which one matches most closely
                currDist = patchDistance(source_patch,target_patch);
                botDist = patchDistance(source_patch_bot,target_patch);
                rightDist = patchDistance(source_patch_right,target_patch);

                % compute the minimum distance between source patches
                [~,lowestInd] = min([currDist botDist rightDist]);

                % if the minimum distance is top or left, we update the NNF
                % to be close to that pixel
                switch lowestInd
                    case 1

                        newNNF(i,j,:) = NNF(i,j,:);

                    case 2
                        % store NNF of bot pixel in best_NNF

                        best_NNF = shift_bot;

                        % check if we can use the pixel above it to assign
                        % to the current pixel

                        if i + best_NNF(1) - 1 >= 1
                            best_NNF(1) = best_NNF(1) - 1;
                        elseif i + best_NNF(1) >= 1
                            best_NNF(1) = best_NNF(1);
                        else
                            best_NNF = NNF(i,j,:);
                        end

                        newNNF(i,j,:) = best_NNF;

                    case 3
                        % store NNF of right pixel in best_NNF

                        best_NNF = shift_right;

                        % check if we can use the pixel to the left of it
                        % to assign the current pixel

                        if j + best_NNF(2) - 1 >= 1
                            best_NNF(2) = best_NNF(2) - 1;
                        elseif j + best_NNF(2) >= 1
                            best_NNF(2) = best_NNF(2);
                        else
                            best_NNF = NNF(i,j,:);
                        end

                        newNNF(i,j,:) = best_NNF;
                end
            end

            % --------------------- Step 3: Random Search ---------------------


            % start with a window that covers the entire image
            dim = max(source_size);

            % coordinates of the center source pixel are:
            curs(1) = newNNF(i, j, 1) + i;
            curs(2) = newNNF(i, j, 2) + j;

            % initialize distance
            distance = patchDistance(pad_source_image(curs(1):curs(1)+patch_size-1, curs(2):curs(2)+patch_size-1, :), target_patch);


            % while the window size is greater than or equal to 3x3
            z = 0;
            while dim >= hole_multiplier*patch_size
                while true
                    % choose a random index within the window around the
                    % current source pixel. This should be centered around the
                    % curs(1) and curs(2) patch coordinates for optimal
                    % searching
                    source_row = max(1, min(source_size(1), randi([curs(1)-floor(dim/2), curs(1)+floor(dim/2)])));
                    source_col = max(1, min(source_size(2), randi([curs(2)-floor(dim/2), curs(2)+floor(dim/2)])));

                    source_patch= pad_source_image(source_row:source_row+patch_size-1, source_col:source_col+patch_size-1, :);

                    if ~any(pad_target_mask(source_row:source_row+patch_size-1, source_col:source_col+patch_size-1))
                        % If no hole pixels found in the source patch, break the loop
                        break;
                    end
                    z = z+1;
                    if z > 20
                        break;
                    end
                end

                dist = patchDistance(source_patch, target_patch);

                % check if distance is better or worse
                if dist <= distance
                    % update minimum distance
                    distance = dist;
                    % update source coordinates
                    curs = [source_row; source_col];
                end

                % update dimensions to be smaller (divide by 2)
                dim = dim / 2;
            end

            % update NNF based on the computed source coordinates
            if z > 20
                NNF(i, j, :) = [curs(1)-i, curs(2)-j];
            else
                NNF(i, j, :) = newNNF(i,j,:);
            end
        end
    end
end

toc

    function distance = patchDistance(source,target)

        % Calculate the squared difference of each channel of patches
        diff_squared = (source - target).^2;

        % Sum the squared differences over all channels without the NaN
        % values, which represents the valid region of the target
        distance = sum(diff_squared(:), 'omitnan');
    end

fprintf("Done!\n");
end

