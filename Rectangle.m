classdef Rectangle < Quadrilater
    properties
        Length
        Width
    end

    methods
        function obj = Rectangle(length, width)
            obj.Length = length;
            obj.Width = width;
        end

        function area = getArea(obj)
            area = obj.Length * obj.Width;
        end
    end
end
