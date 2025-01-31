x = randi([0, 100], 1, 256);

while length(x) > 8
    odd_x = x(1:2:end);
    even_x = x(2:2:end);
    
    % if length(odd_x) > length(even_x)
    %     even_x(end+1:length(odd_x)) = 0;
    % elseif length(even_x) > length(odd_x)
    %     odd_x(end+1:length(even_x)) = 0;
    % end
    
    x = odd_x + even_x;
end

disp(x)