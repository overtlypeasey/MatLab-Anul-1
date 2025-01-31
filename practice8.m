clear
clc
close all

N = 40;
x = -100 + 200 * rand(1, N);
y = -100 + 200 * rand(1, N);

% Compute pairwise distances
points = [x', y'];
distanceMatrix = pdist2(points, points);

% Initialize variables
visited = false(N,1);
path = zeros(N,1);
currentPoint = randi(N);
path(1) = currentPoint;
visited(currentPoint) = true;

avgDistance = mean(pdist(points));

% Build path
for i = 2:N
    [~, idx] = min(abs(distanceMatrix(currentPoint, :) - avgDistance) .* ~visited');
    currentPoint = idx;
    path(i) = currentPoint;
    visited(currentPoint) = true;
end

% Extract ordered coordinates
orderedX = x(path);
orderedY = y(path);

% Plot the path
figure;
plot(orderedX, orderedY, '-o');
title('Path Connecting All Points with Equal Segments');
xlabel('X');
ylabel('Y');
grid on;