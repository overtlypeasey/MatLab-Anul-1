N = 500;

X = -10 + (10 - (-10)) * rand(N, 1);
Y = -10 + (10 - (-10)) * rand(N, 1);

distances = sqrt(X.^2 + Y.^2);

inside_circle = distances <= 10;
X = X(inside_circle);
Y = Y(inside_circle);

distances = sqrt(X.^2 + Y.^2);

outside_inner_circle = distances >= 7;
X = X(outside_inner_circle);
Y = Y(outside_inner_circle);

K = convhull(X, Y);

figure;
scatter(X, Y, 'filled');
hold on;
plot(X(K), Y(K), 'r', 'LineWidth', 2);
title('Path');
xlabel('X');
ylabel('Y');
xlim([-10 10]);
ylim([-10 10]);
grid on;
hold off;