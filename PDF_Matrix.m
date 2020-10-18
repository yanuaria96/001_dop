function f=PDF_Matrix(x)
% x - vector on the unit sphere
ax = 2.8936;
ay = 0.0634;
az = 0.8;
m = 2.5883;
A = [[ax 0 0]; [0 ay 0]; [0 0 az]];
f = (x' * A * x)^m;
end