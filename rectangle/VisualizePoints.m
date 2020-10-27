
function handle = VisualizePoints(x, y, handle)
% input:
% x - x coordinates of points
% y - y coordinates of points
% handle - struct with pointers to the figure, axis and graphic objects
% output:
% updated struct with pointers

% remove old points from picture if points exists

if(handle.circles ~= 0)
    delete(handle.circles)
end

hold on;
subplot(handle.subplot(1));
handle.circles = scatter(x, y, 20, 'k', 'filled');
end