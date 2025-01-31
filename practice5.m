clc
clear
close all

function [r, phi, z] = cartesianToCylindrical(x, y, z)
    r = sqrt(x^2 + y^2);
    phi = atan2(y, x);
end

function [x, y, z] = cylindricalToCartesian(r, phi, z)
    x = r * cos(phi);
    y = r * sin(phi);
end

% theta = linspace(0, 2*pi, 200);
% z = linspace(-1, 1, 200);
% [Theta, Z] = meshgrid(theta, z);
% R = 1 + 0.3 * sin(3*Theta) .* cos(2*Theta); % Enhanced pattern for logo

% [X, Y, Z_cartesian] = cylindricalToCartesian(R, Theta, Z);

% surf(X, Y, Z_cartesian);
% title('MATLAB Logo in Cylindrical Coordinates');
% xlabel('X');
% ylabel('Y');
% zlabel('Z');
% axis equal;
% shading interp;
% colormap jet;

L = 160*membrane(1,100);

f = figure;
ax = axes;

s = surface(L);
s.EdgeColor = 'none';
view(3)

ax.XLim = [1 201];
ax.YLim = [1 201];
ax.ZLim = [-53.4 160];

ax.CameraPosition = [-145.5 -229.7 283.6];
ax.CameraTarget = [77.4 60.2 63.9];
ax.CameraUpVector = [0 0 1];
ax.CameraViewAngle = 36.7;

ax.Position = [0 0 1 1];
ax.DataAspectRatio = [1 1 .9];

l1 = light;
l1.Position = [160 400 80];
l1.Style = 'local';
l1.Color = [0 0.8 0.8];
 
l2 = light;
l2.Position = [.5 -1 .4];
l2.Color = [0.8 0.8 0];

s.FaceColor = [0.9 0.2 0.2];

s.FaceLighting = 'gouraud';
s.AmbientStrength = 0.3;
s.DiffuseStrength = 0.6; 
s.BackFaceLighting = 'lit';

s.SpecularStrength = 1;
s.SpecularColorReflectance = 1;
s.SpecularExponent = 7;

axis off
f.Color = 'black';