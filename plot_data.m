close all; clear all; 

T = csvread("result.csv");
x = linspace(0, 1, size(T, 2));
y = linspace(0, 1, size(T, 1));

p = pcolor(x, y, T);
set(p, 'EdgeColor', 'none');
colorbar;