clear
clc
close all

x_coordinate =  1 + 9.*rand(1, 10);
y_coordinate = 1 + 9.*rand(1, 10);

figure
plot(x_coordinate, y_coordinate, 'rd')
xlabel('x')
ylabel('y')
title("Graph")
axis equal
grid on
for i = 1:10
    text(x_coordinate(1, i), y_coordinate(1, i), num2str(i))
end
hold on
path = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1];
for i = 1:length(path)-1
    plot(x_coordinate([path(i), path(i+1)]), y_coordinate([path(i), path(i+1)]), 'k')
end
plot(x_coordinate([10, 1]), y_coordinate([10, 1]), 'k')
hold off

lin = @(x1, y1, x2, y2, xp, yp) (xp-x1)./(x2-x1) - (yp-y1)./(y2-y1);
lin(x_coordinate(1), y_coordinate(1), x_coordinate(2), y_coordinate(2), ...
    x_coordinate(8), y_coordinate(8))
inter = @(x1, y1, x2, y2, x3, y3, x4, y4) lin(x1, y1, x2, y2, x3, y3) .* lin(x1, ...
    y1, x2, y2, x4, y4)<0;