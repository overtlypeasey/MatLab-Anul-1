clear
clc 
close all

% % Define the symbolic variables
syms x y z 

% Define the system of equations
eq1 = x^2 + y^3 + z^4 == 2;
eq2 = x^3 + y^4 + z^2 == 3;
eq3 = x^4 + y^2 + z^3 == 4;

% Solve the system of equations
solutions = solve(eq1, eq2, eq3, x, y, z);

% Extract the solutions
x_solutions = double(solutions.x);
y_solutions = double(solutions.y);
z_solutions = double(solutions.z);

% Plot the solutions in the complex plane
figure;
hold on;
scatter(real(x_solutions), imag(x_solutions), 'b*'); % Blue star for x
scatter(real(y_solutions), imag(y_solutions), 'ro'); % Red circle for y
scatter(real(z_solutions), imag(z_solutions), 'g.'); % Green dot for z
hold off;

% Add labels and title
xlabel('Real Part');
ylabel('Imaginary Part');
title('Solutions in the Complex Plane');
legend('x solutions', 'y solutions', 'z solutions');
grid on;

syms x y z c

eq4 = x^2 + y^3 + z^4 == 90;
eq5 = x^3 + y^4 + z^2 == 26;
eq6 = x^4 + y^2 + z^3 == c;

solutions = solve([eq4, eq5, eq6], [x, y, z, c]);

x_solutions = double(solutions.x);
y_solutions = double(solutions.y);
z_solutions = double(solutions.z);
c_sol = double(solutions.c);

figure;
hold on;
scatter(real(x_solutions), imag(x_solutions), 'b*');
scatter(real(y_solutions), imag(y_solutions), 'ro');
scatter(real(z_solutions), imag(z_solutions), 'g.');
hold off;

xlabel('Real Part');
ylabel('Imaginary Part');
title('Solutions in the Complex Plane');
legend('x solutions', 'y solutions', 'z solutions');
grid on;