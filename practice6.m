clc
clear
close all

theta_deg = 0:5:360;

% Step 2: Convert degrees to radians
theta_rad = deg2rad(theta_deg);

% Step 3: Calculate x and y coordinates
x = 1 + cos(theta_rad);
y = sin(theta_rad);

z = x + 1j * y;
w = z.^2;

plot(real(z), imag(z), 'b.', real(w), imag(w), 'b.')
hold on
for i = 1:length(z)
    % Plot gray lines connecting z and w
    plot([real(z(i)) real(w(i))], [imag(z(i)) imag(w(i))], 'Color', [0.5 0.5 0.5])
        
    % Plot small red circles for current z and w
    plot(real(z(i)), imag(z(i)), 'ro', 'MarkerSize', 5)
    plot(real(w(i)), imag(w(i)), 'ro', 'MarkerSize', 5)

    if i > 1
        % Draw thick red line connecting previous and current z
        plot([real(z(i-1)) real(z(i))], [imag(z(i-1)) imag(z(i))], 'r-', 'LineWidth', 2)
        
        % Draw thick red line connecting previous and current w
        plot([real(w(i-1)) real(w(i))], [imag(w(i-1)) imag(w(i))], 'r-', 'LineWidth', 2)
    end

    drawnow
    pause(0.05) % Adjust the pause duration for animation speed
end
hold off    