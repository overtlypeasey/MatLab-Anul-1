classdef circle<Shape
    properties
        radius
    end
    methods
        function obj=circle(radius)
            obj.radius=radius
        end
        function area=getArea(obj)
            area=pi*obj.radius^2
        end
        function printArea(shape) 
            fprintf('The area is: %f\n', shape.getArea()); 
        end
        function plotShape(obj)
            theta=linspace(0,2*pi,100);
            x=obj.radius*cos(theta);
            y=obj.radius*sin(theta);
            plot(x,y,'r')
            axis equal
            title('circle')
        end
    end
end