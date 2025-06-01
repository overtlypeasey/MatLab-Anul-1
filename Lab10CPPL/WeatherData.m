classdef WeatherData
    
    properties
        Temperature (1,1) double {mustBeFinite, mustBeGreaterThanOrEqual(Temperature, -60), mustBeLessThanOrEqual(Temperature, 60)}
        Humidity    (1,1) double {mustBeFinite, mustBeGreaterThanOrEqual(Humidity, 0), mustBeLessThanOrEqual(Humidity, 100)} 
        WindSpeed   (1,1) double {mustBeFinite, mustBeGreaterThanOrEqual(WindSpeed, 0)} 
    end

    methods
        function obj = WeatherData(temperature, humidity, windspeed)
            obj.Temperature = temperature;
            obj.Humidity = humidity;
            obj.WindSpeed = windspeed;
        end
    end
end