classdef rectangle2<quadrilateral
    
    methods
        function obj=rectangle2(length,width)
           obj@quadrilateral(length,length,width,width)
        end
        
    end
end