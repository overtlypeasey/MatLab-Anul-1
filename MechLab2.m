%% PROGRAM FOR DETERMINING THE CENTRE OF MASS
%% OF A SYSTEM OF A HOMOGENOUS PLATE

clear all; 
close all; 
clc;

%% Reading the entry data
a = input('> Please write the width of the plate: ');
b = input('> Please write the height of the plate: ');
n = input('> Please write the number of steps: ');

%% Parameter of the parabola
q = b^2 / a;

%% Analytical calculations
Aa = 2 / 3 * a * b; 
xa = 3*a / 5 ; 
ya = 3*b/ 8 ; 

%% Initialising An Sx and Sy

An=0; Sx=0; Sy=0;

%% Numerical calculations

figure(1);
xlim([-0.5, 0.5 + a]);
ylim([-0.5, 0.5 + b]);
hold on;

for i = 0:n-1
 
    x1 = i / n * a;
    y1 = sqrt(q * x1);
    x2 = (i + 1) / n * a;
    y2 = sqrt(q * x2);
    
    A1 = 0.5 * (x2 - x1) * y1;
    xC1 = (2 * x1 + x2) / 3;
    yC1 = y1 / 3;
    
    A2 = 0.5 * (x2 - x1) * y2;
    xC2 = (x1 + 2 * x2) / 3;
    yC2 = (y1 + y2) / 3;
    
    An = An + A1 + A2;
    Sx = Sx + A1 * xC1 + A2 * xC2;
    Sy = Sy + A1 * yC1 + A2 * yC2;
    
    line([x1, x1], [0, y1], 'Color', 'red');
    line([x1, x2], [y1, y2], 'Color', 'green');
    line([x2, x2], [0, y2], 'Color', 'blue');
end


xn = Sx / An;
yn = Sy / An;

plot(xn,yn,'Marker','*','MarkerSize',10,'Color','black')

%% Solution for calculus

fprintf(['Analitical: A=',num2str(Aa),', xC=',num2str(xa),', yC=',num2str(ya),'\n']);
fprintf(['Numerical : A=',num2str(An),', xC=',num2str(xn),', yC=',num2str(yn),'\n']);

