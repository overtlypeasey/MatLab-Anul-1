classdef WeatherStation < handle
    % WeatherStation Handle class aggregating multiple WeatherSensor objects
    %   Maintains aggregated WeatherData and emits events for specific conditions.

    events
        DataUpdated    % Triggered after each aggregation update
        HighWind       % Triggered when aggregated wind speed exceeds fixed threshold
    end

    properties (SetAccess = private)
        Sensors         = WeatherSensor.empty(1,0)  % List of sensors
        AggregatedData  % Instance of WeatherData with current aggregated readings
    end

    methods
        function obj = WeatherStation()
            % Constructor: initialize empty aggregated data
            obj.AggregatedData = WeatherData(0,0,0);
        end

        function addSensor(obj, sensor)
            % ADDSENSOR Add a WeatherSensor and attach listener
            %   Ensures sensor is not already added.
            validateattributes(sensor, {'WeatherSensor'}, {'scalar'});
            % Check for duplicate
            if any(obj.Sensors == sensor)
                warning('WeatherStation:addSensor', 'Sensor is already added.');
                return;
            end
            % Add sensor
            obj.Sensors(end+1) = sensor;
            fprintf('Sensor added. Total sensors: %d\n', numel(obj.Sensors));
            % Listen for data updates, ignore event data
            addlistener(sensor, 'DataUpdated', @(~,~) obj.processUpdate());
        end
    end

    methods (Access = private)
        function processUpdate(obj)
            % PROCESSUPDATE Recompute aggregated data and check conditions
            n = numel(obj.Sensors);
            if n == 0
                return;
            end
            temps = [obj.Sensors.Temperature];
            hums  = [obj.Sensors.Humidity];
            winds = [obj.Sensors.WindSpeed];

            % Compute averages
            avgT = mean(temps);
            avgH = mean(hums);
            avgW = mean(winds);

            % Update aggregated data
            obj.AggregatedData = WeatherData(avgT, avgH, avgW);
            % Notify that data has been updated
            notify(obj, 'DataUpdated');

            % Check high wind condition with fixed threshold = 20 m/s
            if avgW > 20
                notify(obj, 'HighWind');
            end
        end
    end
end
