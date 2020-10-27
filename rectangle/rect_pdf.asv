function rect_pdf()
% PDF - probability distribudion function
t=tic;
global D_X D_Y D_T ETTA W H W_2 H_2
N = 100;
c1_const = 0.000001;
D_X = c1_const; %perturbation incriment
D_Y = c1_const;
D_T = .001*N^(-2);
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
        x = x_old;
        y = y_old;
    else
        D_T = D_T*1.1;
    
%     if maxVelocity < 3e-4
%         D_T = .04*N^(-2)*(maxVelocity/3e-4)^(-1/2);
%     end
    

    % Move each object according to its velocity
    x_old = x;
    y_old = y;
    sumPsi_old = sumPsi;
    
    x = x + vx*D_T;
    y = y + vy*D_T;
    end
    
    if maxVelocity < 3e-10
        break;
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