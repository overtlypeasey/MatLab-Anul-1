% Program for calculating and tracing the trajectory of a point moving in a plane, using Cartesian coordinates
% Se pot defini ecuatii parametrice ale traiectoriei:
%   x(t); y(t) ca functii polinomiale de gradul al doilea 
%
% ====== Setari ale ferestrei grafice ====================
  close all; clc;
  figure(1);
%========Datele problemei ==========================
% Daca miscarea are acceleratie, se dau componentele:
ax=10;          % Acceleration on the Ox direction [m/s^2]
ay=25;           % Acceleration on the Oy direction [m/s^2]
%----------------------------------------------------
x0 = -10;             % Initial position on Ox          [m]
y0 = -20;             % Initial position on Oy           [m]
%---------------------------------------------------
v0x = 30;              % Initial speed on Ox          [m/s]
v0y = 0;              % Initial speed on Oy         [m/s]
%---------------------------------------------------

% Se defineste un set de valori pentru timp: 
 t = 0 : 0.01 : 1.0;       % t =  t0 : time step : final_time
%---------------------------------------------------
% Ecuatiile parametrice pe axe pot fi uniform accelerate.
% Definirea ecuatiilor de miscare:
x = x0 + v0x*t +ax/2*t.^2;              % x(t)  [m]
y = y0 + v0y*t +ay/2*t.^2;              % y(t)  [m]
 
% Calculul vitezelor (vx=dx/dt; vy=dy/dt)
vx=v0x+ax*t;                                  % vx(t)  [m/s]                                         
vy=v0y+ay*t;                                  % vy(t)  [m/s]    
v=sqrt(vx.^2+vy.^2);                      % Modulul vitezei [m/s]
maxv=max(v);                                % Viteza maxima [m/s]

%====== Trasarea traiectoriei punctului si vectori viteza ======
 axis([min(x) max(x) min(y) max(y)]); 
 np=length(t);
 indv=fix(linspace(1,np,11));          % Se vor trasa 11 viteze
 comet (x, y, 0.05);                         % Ca un punct care lasa o urma
 plot(x,y,'Color','g','LineWidth',3); % Retrasare cu o linie mai groasa 
  hold on;
 scv=3/maxv;                                  % Scara trasarii vitezelor
 qp=quiver(x(indv),y(indv),vx(indv)*scv,vy(indv)*scv,0);  % Trasare v
 set(qp,'MaxHeadSize',0.1,'Color','r','LineWidth',3);
  for i=1:11
   text(x(indv(i))+scv/2*vx(indv(i)),y(indv(i))+scv/2*vy(indv(i)), ['v' int2str(i)] ,'fontsize', 20);
  end
 
 % ------ Se scriu titlul, axele, grid, etc. --------
  set(gca, 'linewidth', 2, 'fontsize', 15);
 hold off; axis('equal'); grid; 
 title('Traiectoria & Vectori viteza')
 xlabel('x','FontSize',15); ylabel('y','FontSize',15)

 pause(2)
 %====== Trasarea traiectoriei punctului si vectori acceleratie ======
  figure(2);
  plot(x,y,'Color','g','LineWidth',3); 
  hold on;
  if length(ax)==1                            % Daca ax=const. se genereaza
    ax=ax*ones(size(t));                   %  o matrice cu elementele ax
  end
  if length(ay)==1 
    ay=ay*ones(size(t));
  end 
  maxa=max(sqrt(ax.^2+ay.^2));
  sca=3/maxa;                                  % Scara trasarii acceleratiilor
  qp=quiver(x(indv),y(indv),ax(indv)*scv,ay(indv)*scv,0);  % Trasare a
  set(qp,'MaxHeadSize',0.1,'Color','m','LineWidth',3);
   for i=1:11
   text(x(indv(i))+sca/2*ax(indv(i)),y(indv(i))+sca/2*ay(indv(i)), ['a' int2str(i)] ,'fontsize', 20);
  end
 
 % ------ Se scriu titlul, axele, grid, etc. --------
  set(gca, 'linewidth', 2, 'fontsize', 15);
  hold off; axis('equal'); grid; 
  title('Traiectoria & Vectori acceleratie')
  xlabel('x','FontSize',15); ylabel('y','FontSize',15)
  
%======== Listarea rezultatelor ===========================
 for i=1:np
   fprintf('t=%g [s] x=%g [m]  y=%g [m] vx=%g [m/s] vy=%g [m/s] \n',t(i),x(i),y(i),vx(i),vy(i));
 end