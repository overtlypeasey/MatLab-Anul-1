clear
clc
close all

function plot_compliance_ratios_for_fixed_N()
    % Define the fixed value of N
    N = 10;
    % Number of cases for averaging
    num_cases = 3;
    % Initialize results
    ratios = zeros(num_cases, 1);
    
    % Loop through each case
    for case_idx = 1:num_cases
        % Generate N random points
        points = randi([1, 10], N, 2);
        % Compute the ratio of complying paths
        ratios(case_idx) = compute_compliance_ratio(points);
    end
    
    % Compute the average ratio
    avg_ratio = mean(ratios);
    
    % Plot the result using bar()
    figure;
    bar(N, avg_ratio);
    xlabel('Number of Points (N)');
    ylabel('Average Compliance Ratio (\bar{r})', 'Interpreter', 'latex');
    title(['Average Compliance Ratio for N = ', num2str(N)], 'Interpreter', 'latex');
end

function ratio = compute_compliance_ratio(points)
    % Generate all possible permutations of points
    perms_points = perms(1:size(points, 1));
    num_paths = size(perms_points, 1);
    num_compliant_paths = 0;
    
    for i = 1:num_paths
        perm = perms_points(i, :);
        perm_points = points(perm, :);
        path = [perm_points; perm_points(1, :)]; % Close the path
        if ~path_self_intersects(path)
            num_compliant_paths = num_compliant_paths + 1;
        end
    end
    
    % Compute the ratio of complying paths to all possible paths
    ratio = num_compliant_paths / num_paths;
end

function intersects = path_self_intersects(path)
    n = size(path, 1);
    for i = 1:n-1
        for j = i+2:n-1
            if i == 1 && j == n-1
                continue;
            end
            if lines_intersect(path(i, :), path(i+1, :), path(j, :), path(j+1, :))
                intersects = true;
                return;
            end
        end
    end
    intersects = false;
end

function intersect = lines_intersect(p1, p2, p3, p4)
    intersect = ccw(p1, p3, p4) ~= ccw(p2, p3, p4) && ccw(p1, p2, p3) ~= ccw(p1, p2, p4);
end

function result = ccw(A, B, C)
    result = (C(2) - A(2)) * (B(1) - A(1)) > (B(2) - A(2)) * (C(1) - A(1));
end

% Call the main function to execute the program
plot_compliance_ratios_for_fixed_N();