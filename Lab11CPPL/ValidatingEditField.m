classdef ValidatingEditField < handle
    properties
        EditField % The underlying uieditfield
        object
        ValidColor = [0.74, 0.94, 0.74]; % Default valid color(light green)
        InvalidColor = [1, 0.6, 0.6]; % Reddish color for
        ValidationFcn % Custom validation function
        handle
    end

    methods
        % Constructor
        function obj = ValidatingEditField(parent, varargin)
            % Create the uieditfield object within the parent container
            obj.EditField = uieditfield(parent, varargin{:});
            % Set the default validation function (simple non-empty check)
            obj.ValidationFcn = @(x) ~isempty(x);
            % Set the ValueChangedFcn callback to the customvalidation function
            obj.EditField.ValueChangedFcn = @(src, event)obj.validateInput();
        end

        % Validation function
        function validateInput(obj)
            % Get the current value of the EditField
            value = obj.EditField.Value;
            % Check validity using the custom validation
            isValid = obj.ValidationFcn(value);
            % Change background color based on validity
            if isValid
                obj.EditField.BackgroundColor = obj.ValidColor;
            else
                obj.EditField.BackgroundColor = obj.InvalidColor;
            end
        end

        function setdigitsValidation(obj)
            obj.ValidationFcn = @(x) ~isempty(x) && all(isstrprop(x, 'digit'));
            % You can set additional properties of the EditField
            obj.EditField.Placeholder = 'Enter digits only';
        end
    end

end

