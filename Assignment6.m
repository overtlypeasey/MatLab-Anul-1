clc
clear
close all

%TETRAHEDRON
vertex = [[1, 1, 1]; [1, -1, -1]; [-1, 1, -1]; [-1, -1, 1]];
face = [1, 2, 3; 1, 2, 4; 1, 3, 4; 2, 3, 4];

figure;
patch('Vertices', vertex, 'Faces', face, ...
    'FaceColor', 'blue', 'EdgeColor', 'black', 'FaceAlpha', 0.7);

axis equal;
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('Tetrahedron');
grid on;
view(3);

%CUBE
vertexCube = [[1, 1, 1]; [-1, 1, 1]; [-1, 1, -1]; [1, 1, -1];
    [1, -1, 1]; [-1, -1, 1]; [-1, -1, -1]; [1, -1, -1]];
faceCube = [1, 2, 3, 4; 1, 2, 6, 5; 1, 4, 8, 5; 5, 6, 7, 8; 
    2, 3, 7, 6; 3, 4, 8, 7];

figure
patch('Vertices', vertexCube, 'Faces', faceCube, ...
    'FaceColor', 'green', 'EdgeColor', 'black', 'FaceAlpha', 0.8)
axis equal;
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('Cube');
grid on;
view(3);

%OCTAHEDRON
vertexOcta = [[1, 0, 0]; [-1, 0, 0]; [0, 1, 0]; 
    [0, -1, 0]; [0, 0, 1]; [0, 0, -1]];
faceOcta = [1, 4, 5, 5; 1, 2, 3, 4; 1, 3, 5, 5; 2, 4, 5, 5; 
    1, 4, 6, 6; 1, 3, 6, 6; 2, 4, 6, 6];

figure
patch('Vertices', vertexOcta, 'Faces', faceOcta, ...
    'FaceColor', 'yellow', 'EdgeColor', 'black', 'FaceAlpha', 0.8)
axis equal;
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('Octahedron');
grid on;
view(3);

%ICOSAHEDRON
phi = (1 + sqrt(5))./2;
vertexIcosa = [-1,  phi,  0; 1,  phi,  0; -1, -phi,  0; 1, -phi,  0; 
    0, -1,  phi; 0,  1,  phi; 0, -1, -phi; 0,  1, -phi;
    phi,  0, -1; phi,  0,  1; -phi,  0, -1; -phi,  0,  1;];
faceIcosa = [5, 6, 12; 3, 5, 12; 3, 5, 4; 5, 6, 10;
    3, 4, 7; 7, 8, 9; 9, 4, 10; 4, 5, 10;
    4, 7, 9; 2, 9, 10; 2, 8, 9; 1, 2, 8;
    3, 7, 11; 3, 11, 12; 1, 8, 11; 1, 6, 12;
    1, 2, 6; 2, 6, 10; 1, 11, 12; 8, 7, 11];

figure
plot3(vertexIcosa(:, 1), vertexIcosa(:, 2), vertexIcosa(:, 3), 'r*')
for i = 1:12
    text(vertexIcosa(i, 1), vertexIcosa(i, 2), vertexIcosa(i, 3), ...
        [' ', num2str(i)])
end
hold on
patch('Vertices', vertexIcosa, 'Faces', faceIcosa, ...
    'FaceColor', 'red', 'EdgeColor', 'black', 'FaceAlpha', 0.8)
axis equal;
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('Icosahedron');
grid on;
view(3);
hold off

%DODECAHEDRON