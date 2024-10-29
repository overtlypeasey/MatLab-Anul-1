clear
clc 
close all

x = linspace(-pi, pi);
f1 = x.*sin(3.*x.^2);
f2 = 1./x .* cos(5.*sqrt(abs(x)));
f3 = f1 + f2;
f4 = f1 - f2;
plot(x, f1, 'k')
hold on
plot(x, f2, 'b')
plot(x, f3, '--r', "linewidth", 3)
plot(x, f4, '--g', "LineWidth", 3)
hold off
legend ('f1(x)', 'f2(x)', 'f1(x)+f2(x)', 'f1(x)-f2(x)')
axis([-pi pi -4 4])
xlabel('x')
ylabel('y')
title("Laboratory 3.1")

figure
v1 = rand(1, 10).*100;
v2 = randn(1, 10).*100;
v3 = randi([50, 200], 1, 10);
scatter(v2, v1, v3)
title("Laboratory 3.2")
axis("equal")
xlabel('v2')
ylabel('v1')
legend('v2')

figure
m = [v1; v2; v3];
bar(m', 'grouped')
%hold on
%bar(m(11:20), "", "r")
title("Laboratory 3.3")
xlabel("m")
ylabel("values")
legend({'v1', 'v2', 'v3'})
