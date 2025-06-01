classdef trapezoid < Shape
    properties
        Base1
        Base2
        Height
    end

    methods
        function obj = trapezoid(b1, b2, h)
            obj.Base1 = b1;
            obj.Base2 = b2;
            obj.Height = h;
        end

        function area = getArea(obj)
            area = 0.5 * (obj.Base1 + obj.Base2) * obj.Height;
        end

        function plotShape(obj)
            % Rough plot of a trapezoid
            SideLength=0;
            x = [0, obj.Base1, obj.Base2, 0];
            y = [0, 0, obj.Height, obj.Height];
            fill(x, y, 'y'); % Yellow color
            axis equal;
            title('Trapezoid');
        end
    end
end
