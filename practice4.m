clear
clc
close all

x = randi([1, 9], 1, 7);
y = randi([1, 9], 1, 7);
color = ["b", "g", "r", "c", "m", "y", "k"];

[xGrid, yGrid] = meshgrid(0:0.1:10, 0:0.1:10);

% Initialize matrix to hold the index of the nearest master point
nearestIdx = zeros(size(xGrid));

% Calculate the nearest master point for each grid cell
for i = 1:length(x)
    distances = (xGrid - x(i)).^2 + (yGrid - y(i)).^2;
    if i == 1
        minDistances = distances;
        nearestIdx = i * ones(size(xGrid));
    else
        mask = distances < minDistances;
        minDistances(mask) = distances(mask);
        nearestIdx(mask) = i;
    end
end

% Define RGB colors for each master point
colors = [
    0, 0, 1;    % blue
    0, 1, 0;    % green
    1, 0, 0;    % red
    0, 1, 1;    % cyan
    1, 0, 1;    % magenta
    1, 1, 0;    % yellow
    0, 0, 0     % black
];

% Create an RGB image based on nearest master point
colorImage = reshape(colors(nearestIdx, :), [size(xGrid,1), size(xGrid,2), 3]);

% Display the colored field
imagesc([0 10], [0 10], colorImage);

hold on
for i = 1:length(x)
    scatter(x(i), y(i), 36, 'filled');
    text(x(i)+0.01, y(i)+0.01, char(color(i)), 'FontSize', 12);
end
hold off
rectangle('Position', [0, 0, 10, 10], 'EdgeColor', 'k');
axis([0 10 0 10]);
axis('equal')
title('Coordinates')
grid on