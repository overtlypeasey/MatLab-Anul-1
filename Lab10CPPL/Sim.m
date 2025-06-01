% runSingleSensor.m

% 1) Clear out old class definitions so MATLAB picks up your latest edits
clear classes

% 2) Create one WeatherSensor
sensor = WeatherSensor();

% 3) Create the WeatherStation and Logger
station = WeatherStation();   % aggregates data and fires events
logger  = Logger(station);    % prints every update & high-wind alerts

% 4) Hook that single sensor into the station
station.addSensor(sensor);    % you should see “Sensor added. Total sensors: 1”

% 5) Trigger one measurement (and watch the console)
sensor.update();              % fires DataUpdated → logData

% (Optional) Trigger repeatedly in a tight loop:
for i = 1:5
    pause(1);                 % 1-second gap
    sensor.update();
end
