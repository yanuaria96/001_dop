function [x, y] = GeneratePointsInRectangle(N, w, h)
% Generate random distribution of points insede the rectangle
% intput: N - points cnt; w - rectangle width; h - rectangle height
% oubput: x, y - x and y arrays for points
x = []; y = [];
for k = 1:N 
    x = [x w*rand()];
    y = [y h*rand()];
end  
end
