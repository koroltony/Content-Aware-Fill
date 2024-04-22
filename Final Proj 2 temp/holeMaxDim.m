function [maxDim] = holeMaxDim(alpha)
% read the hole image in
% now find the size of the hole in the image

% extract alpha channel

if all(alpha)
    maxDim = 0;
    return;
end

alpha_dim = size(alpha);

height = alpha_dim(1);
width = alpha_dim(2);

% iterate through the image to find the place where alpha starts being 0

% store the initial columns and rows bordering the hole
row_top = 1;
row_bot = height;

col_left = 1;
col_right = width;

% iterate through all rows and see bounds
i = 1;
% if the current row contains all ones, there is no hole, and the
% boundary should be moved down
while all(alpha(i,:)) && i <= height
    row_top = i+1;
    i = i+1;
end

% do the same thing from bottom to top
i = height;
% if the current row contains all ones, there is no hole, and the
% boundary should be moved down
while all(alpha(i,:)) && i >= 1
    row_bot = i-1;
    i = i-1;
end

% now subtract to get height of hole

hole_height = row_bot-row_top+1;


% now do this for the width

% iterate through all cols and see bounds
i = 1;
% if the current col contains all ones, there is no hole, and the
% boundary should be moved right
while all(alpha(:,i)) && i <= width
    col_left = i+1;
    i = i+1;
end

% do the same thing from right to left
i = width;
% if the current col contains all ones, there is no hole, and the
% boundary should be moved left
while all(alpha(:,i)) && i >= 1
    col_right = i-1;
    i = i-1;
end

% now subtract to get width of hole

hole_width = col_right-col_left+1;



maxDim = max(hole_width,hole_height);
end

