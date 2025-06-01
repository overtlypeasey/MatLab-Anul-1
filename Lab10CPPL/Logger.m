classdef Logger < handle
    methods
        function obj = Logger(station)
            % Constructor: attach to a WeatherStation
            validateattributes(station, {'WeatherStation'}, {'scalar'});
            % Listen for HighWind events
            addlistener(station, 'HighWind', @(~,~) obj.logHighWind(station));
            % Listen for every aggregated data update event
            addlistener(station, 'DataUpdated', @(~,~) obj.logData(station));
        end
    end

    methods (Access = private)
        function logData(~, station)
            % logData Log aggregated weather readings to console
            data = station.AggregatedData;
            fprintf('%s - Data: T=%.2f°C, H=%.1f%%, W=%.2f m/s\n', ...
                    datestr(now, 'yyyy-mm-dd HH:MM:SS'), ...
                    data.Temperature, data.Humidity, data.WindSpeed);
        end

        function logHighWind(~, station)
            % logHighWind Log high wind notification to console
            avgW = station.AggregatedData.WindSpeed;
            fprintf('%s - ALERT: High wind! Avg Speed = %.2f m/s\n', ...
                    datestr(now, 'yyyy-mm-dd HH:MM:SS'), avgW);
        end
    end
end
