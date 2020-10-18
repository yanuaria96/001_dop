function rect_pdf()
% PDF - probability distribudion function
t=tic;
%global D_THETA D_PHI ETTA D_T
global D_X D_Y D_T ETTA W H W_2 H_2
N = 80;
c1_const = 0.0001;
% D_THETA = c1_const*pi/180;
% D_PHI = c1_const*pi/180;
D_X = c1_const;
D_Y = c1_const;
D_T = .04*N^(-2);
ETTA = 1;
W = 2; % picture width
H = 1; % picture height
W_2 = W / 2;
H_2 = H / 2;

[x, y] = GeneratePointsInRectangle(N, W, H);
vx = zeros(1,N); vy = zeros(1,N);

handler = CreateRectangle(W, H);
handler = VisualizePoints(x, y, handler);

for iter = 1:10000
    % Find potentials for each point
    maxVelocity = 0;
    sumPsi = 0;
    for k = 1:N               
        x_cur = x(k);
        y_cur = y(k);
        x_other = [x(1:k-1) x(k+1:end)];
        y_other = [y(1:k-1) y(k+1:end)];
        
        [vx_cur, vy_cur, sumPsi_cur] = FindVelocity(x_cur, y_cur, x_other, y_other);
        
        vx(k) = vx_cur;
        vy(k) = vy_cur;
        
        sumPsi = sumPsi + sumPsi_cur;        
        v2 = vx_cur^2 + vy_cur^2;
        maxVelocity = max([v2 maxVelocity]);
    end
   
    if maxVelocity < 3e-4
        D_T = .04*N^(-2)*(maxVelocity/3e-4)^(-1/2);
    end
    
    if maxVelocity < 3e-10
        break;
    end
    
    % Move each object according to its velocity
    x = x + vx;
    y = y + vy;
    
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
%global D_THETA D_PHI ETTA D_T 
global D_X D_Y D_T ETTA

N = length(x_other);
sumPsi = 0;
vx = 0;
vy = 0;
for k = 1:N
    ox = x_other(k); % other x
    oy = y_other(k); % other y
    
    psiUnmoved = PSI(x_cur, y_cur, ox, oy);

    dPsidX = -(D_T/ETTA).*(1/D_X).*(PSI(x_cur+D_X, y_cur, ox, oy) - psiUnmoved);
    dPsidY = -(D_T/ETTA).*(1/D_Y).*(PSI(x_cur, y_cur+D_Y, ox, oy) - psiUnmoved);

    vx = vx + dPsidX;
    vy = vy + dPsidY;

    sumPsi = sumPsi + psiUnmoved;
end
end

function psi = PSI(x1, y1, x2, y2)
global W_2 H_2
% input: (x1, y1) - point #1, (x2, y2) - point #2
% output: psi

% myu=[0;0;1];
% k=2;
% ri=[sin(thetai).*cos(phii);sin(thetai).*sin(phii);cos(thetai)];
% rj=[sin(thetaj).*cos(phij);sin(thetaj).*sin(phij);cos(thetaj)];
% rij=norm(ri-rj);
%psi=1/(rij^2);
%psi=1/(rij^2)/sqrt(Mises_Fisher(k,myu,ri))/sqrt(Mises_Fisher(k,myu,rj));
%psi=1/(rij^2)/sqrt(Valee_Poussin(k,myu,ri))/sqrt(Valee_Poussin(k,myu,rj));
% psi=1/(rij^2)/sqrt(PDF_Matrix(ri))/sqrt(PDF_Matrix(rj));

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

rij = sqrt((x2-x1).^2 + (y2-y1).^2);

psi = 1/(rij^2)/sqrt(f1*f2);
%psi = 1/(rij^2)/(f1*f2)^(1/4);
%psi = -log(rij)/sqrt(f1*f2);

E = 0.1;
penetration_x = MaccellyBrace(abs(x1 - W_2) - W_2);
penetration_y = MaccellyBrace(abs(y1 - H_2) - H_2);
penetration2 = penetration_x^2 + penetration_y^2;
penalty = -penetration2 * E * 1/2;

psi = psi + penalty;
end

function value = MaccellyBrace(x)
if(x >= 0)
    value = x;
else
    value = 0;
end
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