classdef WeatherSensor < handle
    events
        DataUpdated
    end
    properties  (Access = public)
        Temperature (1,1) double {mustBeFinite, mustBeGreaterThanOrEqual(Temperature, -60), mustBeLessThanOrEqual(Temperature, 60)} = 15
        Humidity    (1,1) double {mustBeFinite, mustBeGreaterThanOrEqual(Humidity, 0), mustBeLessThanOrEqual(Humidity, 100)} = 15
        WindSpeed   (1,1) double {mustBeFinite, mustBeGreaterThanOrEqual(WindSpeed, 0)} = 3
    end

    methods
        function obj=WeatherSensor()

        end

        function update(obj)
            obj.Temperature = -10 + (60+10)*rand();  % Random between -10 and 60 °C
            obj.Humidity    = rand()*100;            % Random between 0 and 100 %%
            obj.WindSpeed   = rand()*30;             % Random between 0 and 30 m/s

            notify(obj, 'DataUpdated')
        end

    end
end