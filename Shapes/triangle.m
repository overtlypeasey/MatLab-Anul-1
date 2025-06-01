classdef triangle < Shape
    properties
        Base
        Height
    end

    methods
        function obj = triangle(base, height)
            obj.Base = base;
            obj.Height = height;
        end

        function area = getArea(obj)
            area = 0.5 * obj.Base * obj.Height;
        end

        function plotShape(obj)
            x = [0, obj.Base, obj.Base/2];
            y = [0, 0, obj.Height];
            fill(x, y, 'm');
            axis equal
            title('Triangle');
        end
        function printArea(shape)
            fprintf('The area is: %.2f\n', shape.getArea());
        end

    end
end
