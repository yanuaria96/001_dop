function rect_pdf()
% PDF - probability distribudion function
t=tic;
global D_X D_Y ETTA W H W_2 H_2
N = 400;
c1_const = 0.000001;
D_X = c1_const; %perturbation incriment
D_Y = c1_const;
ETTA = 1;
W = 2; % picture width
H = 1; % picture height
W_2 = W / 2;
H_2 = H / 2;

[x, y] = GeneratePointsInRectangle(N, W, H);
vx = zeros(1,N); vy = zeros(1,N);

handler = CreateRectangle(W, H);
handler = VisualizePoints(x, y, handler);

position = struct('x', x, 'y', y, 'vx', vx, 'vy', vy, 'sumPsi', 1e100, 'maxVelocity', 0.0, 'D_T', 0.001*N^(-2));
% Vector of positions.  
positions = [position];

for iter = 1:10000

    % Find vx, vy, sumPsi, maxVelosity for each position in positions vector
    for j = 1:length(positions)
        pos = positions(j); 
    	[pos.vx, pos.vy, pos.sumPsi, pos.maxVelocity] = FindVelocities(pos.x, pos.y);
        positions(j) = pos;
    end

    % Sort all positions by sumPsi
    [tmp,indices] = sort([positions.sumPsi]);
    positions = positions(indices);

    % Min
    pos_best = positions(1);
    x = pos_best.x;
    y = pos_best.y;
    vx = pos_best.vx;
    vy = pos_best.vy;
    maxVelocity = pos_best.maxVelocity;
    sumPsi = pos_best.sumPsi;
  

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
    
    
        coef = [0.01 0.8 1.0 1.1 1.2];
    for j = 1:length(coef)

        dt = coef(j) * pos_best.D_T;
        pos.x = pos_best.x + pos_best.vx*dt;
        pos.y = pos_best.y + pos_best.vy*dt;
        pos.D_T = dt;
        
        if j <= length(positions)
            positions(j) = pos;
        else
            positions = [positions pos];
        end
    end
end
tt = toc(t);
disp(tt);
SaveData(x, y);

end
