function [x, y] = GeneratePointsInRectangle(N, w, h)
x = []; y = [];
for k = 1:N 
    x = [x w*rand()];
    y = [y h*rand()];
end  
end
