close all
clc
clear 

xy = [0 0; 3 0; 4 3; 1 4; -1 2; -1.5 0; 0 0];
plot(xy(:, 1), xy(:, 2))

n = length(xy);

xc = sum(xy(:, 1)) / (n-1);
yc = sum(xy(:, 2)) / (n-1);

xmax = max(xy(:, 1));
xmin = min(xy(:, 1));
ymax = max(xy(:, 2));
ymin = min(xy(:, 2));

plot(xy(:, 1), xy(:, 2))
grid on
hold on
plot(xc, yc, 'xr');

nv = 0;
xcm = 0;
ycm = 0;

for k=1:1000
    x0 = xmin+rand()*(xmax-xmin);
    y0 = ymin+rand()*(ymax-ymin);
    

    nintersect = 0;
    for i=1:n-1
        x1 = xy(i, 1);
        y1 = xy(i, 2);
        x2 = xy(i+1, 1);
        y2 = xy(i+1, 2);
        if y0>=ymin && y0<ymax && y2 ~= y1
            xi = (x2-x1)*(y0-y1)/(y2-y1)+x1;
            if x0<xi
                nintersect = nintersect+1;
            end
        end
    end
    in = mod(nintersect, 2) == 0;
    if in 
        % plot(x0, y0, '.g')
        nv = nv+1;
        xcm = xcm+x0;
        ycm = ycm+y0;
    else
        % plot(x0, y0, '.r')
    end
end
xcm = xcm/nv;
ycm = ycm/nv;

plot(xcm, ycm, 'or')

xy(:, 1) = xy(:, 1)-xcm;
xy(:, 2) = xy(:, 2)-ycm;
xamin = [];
yamin = [];
xbmin = []; ybmin = [];
lmin = 1000000000000000000;

for alfa = 0:0.1:180
    if abs(alfa) == 90
        continue
    end
    xa=[];
    ya=[];
    xb=[];
    yb=[];
    for i=1:n-1
        x1 = xy(i, 1);
        y1 = xy(i, 2);
        x2 = xy(i+1, 1);
        y2 = xy(i+1, 2);

        x = (x1*(y2-y1)-y1*(x2-x1))/((y2-y1)-tand(alfa)*(x2-x1));
        if x<max(x1, x2) && x>min(x1, x2)
            if isempty(xa)
                xa = x;
                ya = x*tand(alfa);
            else
                xb = x;
                yb = x*tand(alfa);
            end
        end
    end
    l = sqrt((xa-xb)^2+(ya-yb)^2);
    if l<lmin
        lmin = l;
        xamin = xa;
        yamin = ya;
        xbmin = xb;
        ybmin = yb;
    end
end
plot([xa, xb]+xcm, [ya, yb]+ycm, 'k')
