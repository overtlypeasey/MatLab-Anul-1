clc;
clear;
close all;

% Parameters for the Weierstrass function
a = 0.5;      % Convergence parameter (0 < a < 1)
b = 3;        % Frequency parameter (b > 1, typically an odd integer)
n_max = 200;  % Number of terms in the series

% Verify the condition ab > 1 + (3/(2*pi)) to ensure non-differentiability
if a * b <= 1 + (3/(2*pi))
    error('Parameters do not satisfy the condition ab > 1 + (3/(2*pi)). Please choose different a and b.');
end

% Define zoom levels: 0 to 8 (total 9 subplots)
num_zoom_levels = 8;
total_plots = num_zoom_levels + 1;

% Precompute a^n and b^n up to n_max to improve efficiency
n = 0:n_max;
a_n = a.^n;
b_n = b.^n;

% Initialize figure
figure('Color', 'w', 'Name', 'Weierstrass Function with Successive Zoom-Ins');

% Define a maximum threshold for b^n * pi * x to prevent numerical overflow
threshold = 1e6;

% Loop over each zoom level
for zoom = 0:num_zoom_levels
    % Calculate the scaling factor for the current zoom level
    scale = 10^zoom;
    
    % Define the range for x by reducing the interval by a factor of 10 each zoom
    x_min = -10/scale;
    x_max = 10/scale;
    
    % Define the number of points in the x-axis
    num_points = 10000;  % You can adjust this for higher/lower resolution
    
    % Create the x vector
    x = linspace(x_min, x_max, num_points);
    
    % Initialize the Weierstrass function values to zero
    W = zeros(size(x));
    
    % Compute the Weierstrass function
    for idx = 1:length(n)
        current_n = n(idx);
        current_a_n = a_n(idx);
        current_b_n = b_n(idx);
        
        % Compute the argument for the cosine function
        argument = current_b_n * pi * x;
        
        % Identify elements where |argument| exceeds the threshold
        valid = abs(argument) <= threshold;
        
        % Compute the term only for valid elements to prevent numerical issues
        term = zeros(size(x));
        term(valid) = current_a_n * cos(argument(valid));
        
        % Add the term to the Weierstrass function
        W = W + term;
        
        % Optional: Display progress every 50 iterations
        if mod(current_n, 50) == 0
            fprintf('Zoom Level %d: n = %d/%d\n', zoom, current_n, n_max);
        end
    end
    
    % Create subplot indices (arranged row-wise)
    subplot_rows = 3;
    subplot_cols = 3;
    subplot_idx = zoom + 1;
    
    % Create subplot
    subplot(subplot_rows, subplot_cols, subplot_idx);
    plot(x, W, 'b', 'LineWidth', 1.5);
    xlabel('x', 'FontSize', 10);
    ylabel('F(x)', 'FontSize', 10);
    
    % Title indicating the zoom level
    title(sprintf('Zoom Level %d: ', zoom), 'FontSize', 12);

    
    % Enable grid and adjust axis
    grid on;
    % axis tight;
end

