x = linspace(-5, 5);
y = cos(x);
%plot(x, y, "Y")

%x = linspace(-4, 1);
z = x.^2+3*x-1;
%plot(x, z, "G")

%x = linspace(0, 2*pi);
w = (sin(x))./(x.^2+1);
plot(x, w, "*")
title("W function")
xlabel("x")
ylabel("w")

%x = linspace(0, pi);
t = (x.^2 + 3*x - 1) .^ sin(x);
%plot(x, t, "--R")

%bar(matlabassignment)
%title("Assignment bar")
%xlabel("Position")
%ylabel("Value")