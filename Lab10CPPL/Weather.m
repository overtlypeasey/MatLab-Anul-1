sensor = WeatherSensor();

addlistener(sensor, 'DataUpdated', @(src,~) onData(src))
% 
% function onData(src)
%     wd = WeatherData(src.Temperature, src.Humidity, src.WindSpeed);
% end

sensor.update()

sensor.update()