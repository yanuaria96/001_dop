
function handles = GenerateSphere(R,x,y,z)
handles.fig = figure(20200225);
handles.subplot(1) = subplot(1,2,1);
[x,y,z]=sphere;
handles.sphere = surf(x,y,z);
colormap([0.9  0.9 0.9]);
hold on;
handles.scatter(1) = scatter3(0,0,0,30,[0 0 1],'filled');
handles.scatter(2) = scatter3(0,0,0,30,[1 0 0],'filled');
end