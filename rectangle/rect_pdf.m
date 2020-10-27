function rect_pdf()
% PDF - probability distribudion function
t=tic;
global D_X D_Y D_T ETTA W H W_2 H_2
N = 100;
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

% "Position" structure. Store an information about coordinates, velocities, sumPsi and maxVelocity
position = {x: x, y: y, vx: vx, vy: vy, sumPsi: 1e100, maxVelocity: 0.0};
% Vector of positions.  
positions = [position];

for iter = 1:10000

    % Find vx, vy, sumPsi, maxVelosity for each position in positions vector
    for j = 1:length(positions)
        pos = positions(j); 
    	[pos.vx, pos.vy, pos.sumPsi, pos.maxVelocity] = FindVelocities(pos.x, pos.y);   
    end

    % Sort all positions by sumPsi
    sort(positions.sumPsi);

    % Best
    pos = positions(1);
    x = pos.x;
    y = pos.y;
    maxVelocity = pos.maxVelocity;
    sumPsi = pos.sumPsi;

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


    dt_array = [0.4 0.8 1.0 1.1 1.2];
    for j = 1:length(dt_array)

        pos = positions(1);
        dt = dt_array(j);
        pos.x = pos.x + pos.vx*dt;
        pos.y = pos.y + pos.vy*dt;
        
        if j < length(positions)
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

function [vx, vy, sumPsi, maxVelocity] = FindVelocities(x, y) 
    % Find velocities (vx, and vy), sumPsi and max velocity for each point position (x,y)
    % input: x, y - arrays of x and y coordinates
    % output: vx, vy - arrays of velocities in x and y directions.
    % sumPsi - sum psi
    % maxVelocity - max velocity as vx*vx+vy*vy
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
end