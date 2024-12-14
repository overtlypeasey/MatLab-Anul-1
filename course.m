clear
clc
close all

function time = brachi24(y)
ly = [50 y 0];
time = 0;
v = 0;
g = 9.81;
dx = 10; %xmax = 10*15
for i=2:length(ly)
    dy = ly(i-1)-ly(i);
    ds = sqrt(dx^2+dy^2);
    a = dy/ds*g;
    if abs(a)<0.00001
        dt = ds/v;
    else
        dt = (-v+sqrt(v^2+2*a*ds))/a;
    end
    v = v+a*dt;
    time = time+dt;
end
end
y0 = 90:-10:10;
y = fminsearch(@brachi24, y0);
plot(y, y0)