function rect_pdf()
% PDF - probability distribudion function
t=tic;
global D_X D_Y D_T ETTA W H W_2 H_2
N = 2;
c1_const = 0.000001;
D_X = c1_const; %perturbation incriment
D_Y = c1_const;
D_T = .01*N^(-2);
ETTA = 1;
W = 2; % picture width
H = 1; % picture height
W_2 = W / 2;
H_2 = H / 2;

[x, y] = GeneratePointsInRectangle(N, W, H);
vx = zeros(1,N); vy = zeros(1,N);

handler = CreateRectangle(W, H);
handler = VisualizePoints(x, y, handler);

sumPsi_old = 1e100;
for iter = 1:10000
    % Find potentials for each point
    maxVelocity = 0;
    sumPsi = 0;
    for k = 1:N               
        x_cur = x(k);
        y_cur = y(k);
        x_other = [x(1:k-1) x(k+1:end)];
        y_other = [y(1:k-1) y(k+1:end)];
%         x_other = x(1:k-1);
%         y_other = y(1:k-1);

        
        [vx_cur, vy_cur, sumPsi_cur] = FindVelocity(x_cur, y_cur, x_other, y_other);
        
        vx(k) = vx_cur;
        vy(k) = vy_cur;
        
        sumPsi = sumPsi + sumPsi_cur;        
        v2 = vx_cur^2 + vy_cur^2;
        maxVelocity = max([v2 maxVelocity]);
    end
   
    if sumPsi>sumPsi_old
        D_T = D_T/1.1;
%         x = x_old;
%         y = y_old;
    else
        D_T = D_T*1.1;
    
%     if maxVelocity < 3e-4
%         D_T = .04*N^(-2)*(maxVelocity/3e-4)^(-1/2);
%     end
    
    if maxVelocity < 3e-10
        break;
    end
    
    % Move each object according to its velocity
    x_old = x;
    y_old = y;
    sumPsi_old = sumPsi;
    
    x = x + vx*D_T;
    y = y + vy*D_T;

    end
    
    
    if rem(iter,1) == 0
        handler = VisualizePoints(x, y, handler);
        leftCnt = length(find(x < W_2));
        disp(['LEFT SIDE: ' num2str(leftCnt) ', RIGHT SIDE: ' num2str(N - leftCnt)]);
        pause(0.01);
    end
    
    if rem(iter,1) == 0
        disp([...
            'iter: ' num2str(iter) ...
            ', maxVelocity: ' num2str(maxVelocity) ...
            ', sumPsi: ' num2str(sumPsi)]);
    end
    
    if rem(iter,1) == 0
        SaveData(x, y);
    end
end
tt = toc(t);
disp(tt);
SaveData(x, y);

end

function SaveData(x, y)
dlmwrite('200_anis_x_y_z.txt', [x' y'],'delimiter','\t');
save 'rect_pdf_output';
end

function obj = MoveObject(obj)
    v = obj.velocity;
    obj.phi = obj.phi + v.phi;
    obj.theta = obj.theta + v.theta;
end

function [vx, vy, sumPsi] = FindVelocity(x_cur, y_cur, x_other, y_other)
global D_X D_Y ETTA

N = length(x_other);
sumPsi = 0;
vx = 0;
vy = 0;
for k = 1:N
    ox = x_other(k); % other x
    oy = y_other(k); % other y
    
    psiUnmoved = PSI(x_cur, y_cur, ox, oy);
    
    dPsidX = -(1/ETTA)*(1/D_X)*(PSI(x_cur+D_X, y_cur, ox, oy) - psiUnmoved);
    dPsidY = -(1/ETTA)*(1/D_Y)*(PSI(x_cur, y_cur+D_Y, ox, oy) - psiUnmoved);

    vx = vx + dPsidX;
    vy = vy + dPsidY;

    sumPsi = sumPsi + psiUnmoved;
end
PsiPenUnmoved = PsiPenalty(x_cur, y_cur);
dPsidX = -(1/ETTA)*(1/D_X)*(PsiPenalty(x_cur+D_X, y_cur) - PsiPenUnmoved);
dPsidY = -(1/ETTA)*(1/D_Y)*(PsiPenalty(x_cur, y_cur+D_Y) - PsiPenUnmoved);
vx = vx + dPsidX;
vy = vy + dPsidY;
sumPsi = sumPsi + PsiPenUnmoved;
end

function psi = PSI(x1, y1, x2, y2)
global W_2 H_2
% input: (x1, y1) - point #1, (x2, y2) - point #2
% output: psi

if(x1 < W_2)
    f1 = 2;
else
    f1 = 1;
end
if(x2 < W_2)
    f2 = 2;
else
    f2 = 1;
end

rij = sqrt((x2-x1)^2 + (y2-y1)^2);

psi = 1/(rij^2)/sqrt(f1*f2);
%psi = 1/(rij^2)/(f1*f2)^(1/4);
%psi = -log(rij)/sqrt(f1*f2);
end



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

function handles = CreateRectangle(w, h)
% In this function we create a struct of pointers
handles.fig = figure(20201018); % a pointer that points to the figure
handles.subplot(1) = subplot(1,1,1); % a pointer that points to the axis
handles.rect = rectangle('Position', [0, 0, w, h]); % a point that points to the rectangle lines
handles.circles = 0; % There are no any circles, set 0
end