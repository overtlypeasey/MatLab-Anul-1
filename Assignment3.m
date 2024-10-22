clear
clc 
close all

x = linspace(-pi, pi);
f1 = x.*sin(3.*x.^2);
f2 = 1./x .* cos(5.*sqrt(x));
plot(x, f1, 'g')
hold on
plot(x, f2, '*')
hold off
title("Laboratory 3.1")

