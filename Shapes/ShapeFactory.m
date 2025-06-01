classdef ShapeFactory
    methods (Static)
        function shape = createShape(shapeType, varargin)
            switch lower(shapeType)
                case 'circle'
                    shape = circle(varargin{:}); % e.g. radius
                case 'rectangle'
                    shape = rectangle2(varargin{:}); % e.g. length, width
                case 'triangle'
                    shape = triangle(varargin{:}); % e.g. base, height
                case 'trapezoid'
                    shape = trapezoid(varargin{:}); % e.g. base1, base2, height
                % add other shapes here
                otherwise
                    error('Unknown shape type: %s', shapeType);
            end
        end
    end
end