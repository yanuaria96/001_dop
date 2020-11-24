
function handles = GenerateSphere(R,x,y,z)
handles.fig = figure(20200225);
handles.subplot(1) = subplot(1,2,1);
[x,y,z]=sphere;
handles.sphere = surf(x,y,z);
colormap([0.9  0.9 0.9]);
hold on;
handles.scatter(1) = scatter3(0,0,0,30,[0 0 1],'filled');
handles.scatter(2) = scatter3(0,0,0,30,[1 0 0],'filled');
xlabel('x');
ylabel('y');
zlabel('z');

% Show 2d circle
handles.subplot(2) = subplot(1,2,2); % subplot for 2D circle
hold on;
r=1;
x0=0;
y0=0;
theta = linspace(0,2*pi,100);
plot(x0 + r*cos(theta),y0 + r*sin(theta),'k-');
plot(x0 + r*cos(theta)*sqrt(3)/2, y0 + r*sin(theta)*sqrt(3)/2,'k-');

handles.scatter(3) = scatter(0,0,30, [0 0 1], 'filled');
handles.scatter(4) = scatter(0,0,30, [1 0 0], 'filled');
xlabel('x'); 
ylabel('z');
end