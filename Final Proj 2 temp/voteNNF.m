function output = voteNNF(NNF, source_image, target_image)
    global patch_size;
    padsize = (patch_size-1)/2;
    pad_source_image = padarray(source_image,[padsize padsize],"symmetric");
    pad_target_image = padarray(target_image,[padsize padsize],nan,"both");

    fprintf("Voting to reconstruct the final result...\n");

    target_size = size(NNF);

    output = zeros(target_size(1)+2*padsize, target_size(2)+2*padsize, 3);
    weights = zeros(target_size(1)+2*padsize, target_size(2)+2*padsize);

    % loop through all indices in NNF and find the source patches. Then
    % combine them in the target image
    for i = 1:target_size(1)
        for j = 1:target_size(2)
            offset = NNF(i,j,:);
            yoffset = offset(1);
            xoffset = offset(2);

            % find source patch using shifted index of target patch
            source_patch = pad_source_image(i+yoffset:i+yoffset+patch_size-1,j+xoffset:j+xoffset+patch_size-1,:);
            target_patch = pad_target_image(i:i+patch_size-1,j:j+patch_size-1,:);

            distance = patchDistance(source_patch,target_patch);
            % find the closest patch on the border of the target region
            % (AKA find the closest patch where the NNF is not 0)
            closest_patch = pad_target_image(i:i+patch_size-1,j:j+patch_size-1,:);
            distance = distance + patchDistance(closest_patch,source_patch);
            % weight smaller distances more than large distances (use eps
            % to avoid divide by 0 errors in case of tiny difference
            weight = 1 / (distance + eps); 

            % Weight the source patch by its distance and add it to the output
            output(i:i+patch_size-1,j:j+patch_size-1,:) = output(i:i+patch_size-1,j:j+patch_size-1,:) + weight * source_patch;
            % add to the cumulative weight of that pixel
            weights(i:i+patch_size-1,j:j+patch_size-1) = weights(i:i+patch_size-1,j:j+patch_size-1) + weight;
        end
    end

    % Normalize the output by dividing by the sum of weights
    output = output ./ weights;

    output = output(padsize+1:end-padsize,padsize+1:end-padsize,:);

    fprintf("Done!\n");

        function distance = patchDistance(source,target)

        % Calculate the squared difference of each channel of patches
        diff_squared = (source - target).^2;

        % Sum the squared differences over all channels without the NaN
        % values, which represents the valid region of the target
        distance = sum(diff_squared(:), 'omitnan');
    end
end

% % use the NNF to vote the source patches
% function output = voteNNF(NNF, source_image)
%     global patch_size;
%     padsize = (patch_size-1)/2;
%     pad_source_image = padarray(source_image,[padsize padsize],"symmetric");
% 
%     fprintf("Voting to reconstruct the final result...\n");
% 
%     target_size = size(NNF);
% 
%     output = zeros(target_size(1)+2*padsize, target_size(2)+2*padsize, 3);
%     avg = zeros(target_size(1),target_size(2));
% 
%     % write your code here to reconstruct the output using source image
%     % patches
% 
%     % loop through all indices in NNF and find the source patches. Then
%     % combine them in the target image
% 
%     for i = 1:target_size(1)
%         for j = 1:target_size(2)
%             offset = NNF(i,j,:);
%             yoffset = offset(1);
%             xoffset = offset(2);
% 
%             % find source patch using shifted index of target patch
%             source_patch = pad_source_image(i+yoffset:i+yoffset+patch_size-1,j+xoffset:j+xoffset+patch_size-1,:);
% 
%             % splat the patch in the output image
% 
%             output(i:i+patch_size-1,j:j+patch_size-1,:) = output(i:i+patch_size-1,j:j+patch_size-1,:) + source_patch;
% 
%             % figure out total number of overlapping patches at this pixel
%             m = patch_size;
%             n = patch_size;
% 
%             if i-(padsize+1) < 1
%                 m = patch_size - abs(i-(padsize+1));
%             elseif target_size(1)-i < padsize
%                 m = patch_size-padsize+abs(target_size(1)-i);
%             end
% 
%             if j-(padsize+1) < 1
%                 n = patch_size - abs(j-(padsize+1));
%             elseif target_size(2)-j < padsize
%                 n = patch_size-padsize+abs(target_size(2)-j);
%             end
% 
%             avg(i,j) = n*m;
%         end
%     end
%     output = output(padsize+1:end-padsize,padsize+1:end-padsize,:);
%     output = (output)./avg;
% 
%     fprintf("Done!\n");
% end

% function output = voteNNF(NNF, source_image, target_image)
%     global patch_size;
%     padsize = (patch_size-1)/2;
%     pad_source_image = padarray(source_image,[padsize padsize],"symmetric");
%     pad_target_image = padarray(target_image,[padsize padsize],nan,"both");
% 
%     fprintf("Voting to reconstruct the final result...\n");
% 
%     target_size = size(NNF);
% 
%     % Precompute target patches
%     target_patches = zeros(target_size(1), target_size(2), patch_size, patch_size, 3);
%     for i = 1:target_size(1)
%         for j = 1:target_size(2)
%             offset = NNF(i, j, :);
%             yoffset = offset(1);
%             xoffset = offset(2);
%             target_patches(i, j, :, :, :) = pad_target_image(i:i+patch_size-1, j:j+patch_size-1, :);
%         end
%     end
% 
%     output = zeros(target_size(1)+2*padsize, target_size(2)+2*padsize, 3);
%     weights = zeros(target_size(1)+2*padsize, target_size(2)+2*padsize);
% 
%     % loop through all indices in NNF and find the source patches. Then
%     % combine them in the target image
%     for i = 1:target_size(1)
%         for j = 1:target_size(2)
%             offset = NNF(i, j, :);
%             yoffset = offset(1);
%             xoffset = offset(2);
% 
%             % Find source patch using shifted index of target patch
%             source_patch = pad_source_image(i+yoffset:i+yoffset+patch_size-1, j+xoffset:j+xoffset+patch_size-1, :);
% 
%             % Calculate distance between source and target patches
%             distance_to_target = patchDistance(source_patch, target_patches(i, j, :, :, :));
% 
%             % Coherence factor
%             coherence = coherenceFactor(NNF, i, j, target_image);
% 
%             % Combine distance to target and coherence into weight
%             weight = 1 / (distance_to_target + eps) * coherence;
% 
%             % Weight the source patch by its distance and add it to the output
%             output(i:i+patch_size-1, j:j+patch_size-1, :) = output(i:i+patch_size-1, j:j+patch_size-1, :) + weight * source_patch;
%             % Add to the cumulative weight of that pixel
%             weights(i:i+patch_size-1, j:j+patch_size-1) = weights(i:i+patch_size-1, j:j+patch_size-1) + weight;
%         end
%     end
% 
%     % Normalize the output by dividing by the sum of weights
%     output = output ./ weights;
% 
%     output = output(padsize+1:end-padsize, padsize+1:end-padsize, :);
% 
%     fprintf("Done!\n");
% 
%     function distance = patchDistance(source, target)
%         % Calculate the squared difference of each channel of patches
%         diff_squared = (source - target).^2;
% 
%         % Sum the squared differences over all channels without the NaN
%         % values, which represents the valid region of the target
%         distance = sum(diff_squared(:), 'omitnan');
%     end
% 
%     function coherence = coherenceFactor(NNF, i, j, target_image)
%         pad_target_image = padarray(target_image,[padsize padsize],nan,"both");
% 
%         % Initial search radius
%         initial_radius = 3;
% 
%         % Initialize coordinates of nearest neighbor
%         nearest_i = i;
%         nearest_j = j;
% 
%         % Perform radius search until a zero offset patch is found
%         radius = initial_radius;
%         found = false;
%         while ~found
%             % Iterate over the neighborhood within the current radius
%             for di = -radius:radius
%                 for dj = -radius:radius
%                     % Calculate coordinates in NNF
%                     ni = i + di;
%                     nj = j + dj;
% 
%                     % Ensure within bounds
%                     if ni >= 1 && ni <= size(NNF, 1) && nj >= 1 && nj <= size(NNF, 2)
%                         % Check if the offset is 0
%                         if all(NNF(ni, nj, :) == [0 0])
%                             % Update nearest neighbor
%                             nearest_i = ni;
%                             nearest_j = nj;
%                             found = true;
%                             break;
%                         end
%                     end
%                 end
%                 if found
%                     break;
%                 end
%             end
%             % Expand search radius if zero offset patch not found
%             if ~found
%                 radius = radius + 1;
%             end
%         end
% 
%         % Calculate distance to the source patch
%         distance = patchDistance(pad_target_image(i:i+patch_size-1, j:j+patch_size-1, :), ...
%                                  pad_target_image(nearest_i:nearest_i+patch_size-1, nearest_j:nearest_j+patch_size-1, :));
% 
%         % Coherence factor inversely proportional to closest distance
%         coherence = 1 / (distance + eps);
%     end
% end


