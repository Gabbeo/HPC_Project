close all; clear all; 

fileID = fopen('result.bin','r');
dim = fread(fileID, 2, 'uint');
T = fread(fileID, [dim(1), dim(2)], 'double');

y = linspace(0, 1, size(T, 1));
x = linspace(0, 1, size(T, 2));

p = pcolor(y, x, T);
set(p, 'EdgeColor', 'none');
colorbar;