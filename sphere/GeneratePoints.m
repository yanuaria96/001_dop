function [particles, duplicates] = GeneratePoints(N)
for k = 1:N
    theta = pi * rand();
    phi = 2.0 * pi * rand();  
    
    [thetad, phid] = GetMirrored(theta,phi);
    
    particles(k) = struct('theta', theta, 'phi', phi, 'velocity', 0.0, 'isDuplicate', 0);
    duplicates(k) = struct('theta', thetad, 'phi', phid, 'velocity', 0.0, 'isDuplicate', 1);
end  
end