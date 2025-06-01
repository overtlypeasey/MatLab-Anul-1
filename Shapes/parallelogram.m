classdef parallelogram < trapezoid
    properties
        sideLength
    end
    methods
        function obj = parallelogram(base,sideLength,height)
            obj@trapezoid(base, sideLength, height);
            obj.sideLength=sideLength;
           
        end
    end
end
