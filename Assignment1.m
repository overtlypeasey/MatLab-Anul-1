clear
clc
close all

x = linspace(-2*pi, 2*pi);
y = cos(x);
figure
plot(x, y);
title("Cosinus")

x = linspace(-4, 1);
z = x.^2+3*x-1;
plot(x, z, "G")

x = linspace(0, 2*pi);
w = (sin(x))./(x.^2+1);
figure
plot(x, w, "*")
title("W function")
xlabel("x")
ylabel("w")

x = linspace(0, pi);
t = (x.^2 + 3*x - 1) .^ sin(x);
figure
plot(x, t, "--R")