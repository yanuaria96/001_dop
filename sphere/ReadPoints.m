function [particles, duplicates,N] = ReadPoints()

fname = 'D:\YanaOLD\Documents\01work\08-osnPDF\osnPDF\001_spherical_points.txt';

% arrays
theta = []; phi = [];

% Read fortran generated file
f_id = fopen(fname,'r');

% fgetl - get line from the file
line = fgetl(f_id); % skip the first line (head)
line = fgetl(f_id);

while ischar(line)
    numbers = str2num(line);
    theta = [theta   numbers(1)];
    phi   = [phi     numbers(2)];
    line = fgetl(f_id);
end
fclose(f_id);

N = length(theta);

% Check if numbers count is odd!
assert(mod(N,2) == 0, 'Check if odd error!');

for k = 1:N/2
    particles(k) = struct('theta', theta(k), 'phi', phi(k), 'velocity', 0.0, 'isDuplicate', 0);
    duplicates(k) = struct('theta', theta(k+N/2), 'phi', phi(k+N/2), 'velocity', 0.0, 'isDuplicate', 1);
    
    % check 
    [thd, phd] = GetMirrored(theta(k), phi(k));
    
    
end  
N =N/2;
end