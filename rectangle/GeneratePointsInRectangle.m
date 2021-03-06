function [x, y] = GeneratePointsInRectangle(N, w, h)
% Generate random distribution of points insede the rectangle
% intput: N - points cnt; w - rectangle width; h - rectangle height
% oubput: x, y - x and y arrays for points
b = sqrt(2)/(1+sqrt(2));
x = []; y = [];

for k = 1:N
    a = rand();
    if a<b
        x = [x w/2*rand()];
        y = [y h*rand()];
    else
        x = [x w/2*rand()+w/2];
        y = [y h*rand()];
    end
end  
end
