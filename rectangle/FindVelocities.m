function [vx, vy, sumPsi, maxVelocity] = FindVelocities(x, y) 
    % input: x, y - arrays of x and y coordinates
    % output: vx, vy - arrays of velocities in x and y directions.
    % sumPsi - sum psi
    % maxVelocity - max velocity as vx*vx+vy*vy
    N = length(x);
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