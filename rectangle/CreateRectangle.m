function handles = CreateRectangle(w, h)
% In this function we create a rectangle and a struct of pointers
% input: w - width, h - height
% output: handlers
handles.fig = figure(20201018); % a pointer that points to the figure
handles.subplot(1) = subplot(1,1,1); % a pointer that points to the axis
handles.rect = rectangle('Position', [0, 0, w, h]); % a point that points to the rectangle lines
handles.circles = 0; % There are no any circles, set 0
end