clear
clc 
close all

x = linspace(-pi, pi, 10000);
f1 = x .* sin(3 .* x.^2);
f2 = x./pi;
f = f1-f2;
sign_changes = find(diff(sign(f)));
solutions = x(sign_changes);
plot(x, f1)
hold on
plot(x, f2, 'r')
scatter(solutions, f(sign_changes), 'dk', 'filled')
hold off
xlabel('x')
ylabel('y')
title('Laboratory 4.1')
axis('equal')
legend('f1', 'f2', 'solutions')
grid on


v = randi([1, 20], 1, 100);
v1 = v(find(v<7));
v2 = v(find(v>14));
m = v1' * v2;

figure
[row, col] = find(sqrt(m)==floor(sqrt(m)) & m ~= 0);
scatter(row, col, "filled");
axis('equal')
xlabel('row')
ylabel('col')
title('Laboratory 4.2')
legend('Perfect Square')