
function handle = VisualizePoints(points,handle)
delete(handle.scatter);
hold on;
N = length(points);
X = [];
Y = [];
Z = [];
Xd = [];
Yd = [];
Zd = [];
for k = 1:N
    [x,y,z] = Sph2Cart(1,points(k).theta,points(k).phi);
    
    if ~points(k).isDuplicate
        X = [X x];
        Y = [Y y];
        Z = [Z z];
    else
        Xd = [Xd x];
        Yd = [Yd y];
        Zd = [Zd z];
    end
end
subplot(handle.subplot(1));
handle.scatter(1) = scatter3(X,Y,Z,30,[0 0 1],'filled');
handle.scatter(2) = scatter3(Xd,Yd,Zd,30,[1 0 0],'filled');

subplot(handle.subplot(2)); % subplot for 2D circle
handle.scatter(3) = scatter(X,Z,30, [0 0 1], 'filled');
handle.scatter(4) = scatter(Xd,Zd,30, [1 0 0], 'filled');
end