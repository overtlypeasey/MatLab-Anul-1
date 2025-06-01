classdef quadrilateral < Shape
    properties
        Side1
        Side2
        Side3
        Side4
      
    end
    
    methods 
        function obj = quadrilateral(s1, s2, s3, s4)
            obj.Side1 = s1;
            obj.Side2 = s2;
            obj.Side3 = s3;
            obj.Side4 = s4;
           
        end
        
        function area = getArea(obj)
           area=obj.Side1*obj.Side3;
        end
    function plotShape(obj)
    % Define four vertices
    x = [0, obj.Side1, obj.Side2, 0];
    y = [0, 0, obj.Side3, obj.Side4];
    
    % Close the shape
    x(end+1) = x(1);
    y(end+1) = y(1);

    % Plot
    fill(x, y, 'm'); % Magenta color
    axis equal
  
    xlabel('X-axis');
    ylabel('Y-axis');
end

    end
end
