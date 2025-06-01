% Program for calculating and plotting the connecting rod-crank mechanism   
% Prof. M.V. Predoi, Departamentul de Mecanica, UPB.
% Se pot defini geometria mecanismului si rotatia manivelei:
%    theta(t) 
%
% ====== Setari ale ferestri grafice ====================
close all; clc;

%========Datele problemei ==========================
R=0.3;         % Crank radius [m]
L= R+0.4;    % Connecting rod length [m]
if (L<R)
  disp('EROARE: L<R !')
  return
 end
%---------------------------------------------------
% Daca miscarea are acceleratie:
eps=-250;              % Angular acceleration           [rad/s^2]
%----------------------------------------------------
theta0=-pi/2;       % Initial angular position    [rad]
%---------------------------------------------------
omega0=10;    % Initial angular velocity  [rad/s]
%---------------------------------------------------

% Se defineste un set de valori pentru timp: 
 t = linspace(0,0.2*pi/omega0,45);       % t0 : time step :final time
 unu=ones(size(t));
 %---------------------------------------------------
 % Definirea ecuatiilor de miscare:
 theta = theta0 + omega0*t + eps/2*t.^2;      % theta(t)  [rad]
 %---------------------------------------------------
 % Definirea componentelor vitezei:
 omega=omega0 + eps*t;                              % Viteza unghiulara (rad/s)
 vr = 0;                                                              % dr/dt     [m/s]
 vn = omega.*R;                                              % omega*r  [m/s]
  %---------------------------------------------------
 % Definirea componentelor acceleratiei:
 epsilon= eps;
 ar =- omega.^2*R;                                                              % d2r/dt2 - omega^2*r  [m/s^2]
 an = R.*epsilon +2*vr.*omega;                     % epsilon*r +2 dr/dt*omega  [m/s^2]
 %----------------------------------------------------
 % Transformarea de coordonate: polare -> careziene
 xA=R.*cos(theta);
 yA=R.*sin(theta);
 vxA=vr.*cos(theta) - vn.*sin(theta);
 vyA=vr.*sin(theta) + vn.*cos(theta);
 vA=sqrt(vxA.^2 + vyA.^2);                                % Modulul vitezei in A
 axA=ar.*cos(theta) - an.*sin(theta);
 ayA=ar.*sin(theta) + an.*cos(theta);
aA=sqrt(axA.^2 + ayA.^2);                                % Modulul acceleratiei in A
 
 %----------------------------------------------------
 % Calcul pozitie si viteza piston
 % Not: fi = unghi dintre Ox si biela
 sinfi=R/L*sin(theta);  % Teorema sinusurilor
 fi=asin(sinfi);           % Unghiul fi
 omega2=R/L*cos(theta)./cos(fi).*omega;
 xB=R*cos(theta)+L*cos(fi); % Proiectie B pe Ox
 yB=zeros(size(xB));               % B ramane pe axa Ox
 vxB=-R*sin(theta).*(unu+R/L*cos(theta)./cos(fi)).*omega;  % Vx = dx/dt
 vyB=zeros(size(vxB));           % B se misca doar pe Ox
 vB=sqrt(vxB.^2 + vyB.^2);                                % Modulul vitezei in B
 axB=-R*omega.^2.*(cos(theta)+R/L*cos(2*theta)) - R*epsilon.*sin(theta).*(1+R/L*cos(theta));
 ayB=zeros(size(axB));           % B se misca doar pe Ox
 aB=sqrt(axB.^2 + ayB.^2);                                % Modulul acceleratiei in B
 %====== Trasarea miscarii bielei si vectori viteza ======
 figure(1);
 np=length(t);
 indv=fix(linspace(1,np,11));          % Se vor trasa 11 viteze
 origin=zeros(size(xA)); 
 xvecA=[origin xA];
 yvecA=[origin yA];
 xvecB=[xA xB];
 yvecB=[yA yB]; 
  
% plot(x,y,'Color','g','LineWidth',3); % Retrasare cu o linie mai groasa
% hold on; 
 maxv=max(sqrt(vr.^2+vn.^2));
 scv=0.05/maxv;                                  % Scara trasarii vitezelor
 dx=scv;
 dy=0.7*scv;
 axis([-2*R 6*R -4*R 4*R]);
 hold on
 
 for i=1:np
     line('xdata', [0.5*R 1.5*R+L], 'ydata', [0 0]);  % Trasare directie culisa (Ox)
     
     % Trasare manivela
     line ('xdata', [xvecA(i) xvecA(np+i)] , 'ydata', [yvecA(i) yvecA(np+i)],'Color','k','LineWidth',3,'marker','o');
     qp=quiver(xA(i),yA(i),vxA(i)*scv,vyA(i)*scv,0);  % Trasare vA
     set(qp,'autoscalefactor',scv,'MaxHeadSize',0.1,'Color','r','LineWidth',3);
     vAs=sprintf('%4.1f', vA(i)); 
     text(xA(i)+scv/2*vxA(i),yA(i)+scv/2*vyA(i), vAs ,'Color','b','fontsize', 15);  % Marcare v
     plot(0,-0.001,'marker','^','markerfacecolor','k','markersize',12);             % desen articulatie manivela
     
     %Trasare culisa 
      patch([xB(i)-dx,xB(i)+dx,xB(i)+dx,xB(i)-dx],[yB(i)-dy,yB(i)-dy,yB(i)+dy,yB(i)+dy],'g');
    
    % Trasare biela
     line ('xdata', [xvecB(i) xvecB(np+i)] , 'ydata', [yvecB(i) yvecB(np+i)],'Color','B','LineWidth',3,'marker','o');
     qp=quiver(xB(i),yB(i),vxB(i)*scv,vyB(i)*scv,0);  % Trasare vB
     set(qp,'autoscalefactor',scv,'MaxHeadSize',0.1,'Color','r','LineWidth',3);
     vBs=sprintf('%4.1f', vB(i)); 
     text(xB(i)+scv/2*vxB(i),yB(i)+scv/2*vyA(i), vBs ,'Color','b','fontsize', 15);  % Marcare v
    
     % Pregatire grafic noua pozitie
     pause(0.1);            % Pauza intre secvente
     if (i<np)
        cla;                       % Sterge imaginea pentru noul desen
     end
     
 end

 % ------ Se scriu titlul, axele, grid, etc. --------
 set(gca, 'linewidth', 2, 'fontsize', 15);
 hold off;  axis('equal');grid;
 title('Pozitii mecanism & Vectori viteza')
 xlabel('x','FontSize',15); ylabel('y','FontSize',15)
 pause(2);

 %=============================================================================================
 
 %====== Trasarea miscarii bielei si vectori acceleratie ======
 figure(2);
 np=length(t);
 indv=fix(linspace(1,np,11));          % Se vor trasa 11 viteze
 origin=zeros(size(xA)); 
 xvecA=[origin xA];
 yvecA=[origin yA];
 xvecB=[xA xB];
 yvecB=[yA yB]; 
  
% plot(x,y,'Color','g','LineWidth',3); % Retrasare cu o linie mai groasa
% hold on; 
 maxa=max(sqrt(axA.^2+ayA.^2));
 sca=0.05/maxa;                                  % Scara trasarii vitezelor
 dx=scv;
 dy=0.7*scv;
 axis([-2*R 6*R -4*R 4*R]);
 hold on
 
 for i=1:np
     line('xdata', [0.5*R 1.5*R+L], 'ydata', [0 0]);  % Trasare directie culisa (Ox)
     
     % Trasare manivela
     line ('xdata', [xvecA(i) xvecA(np+i)] , 'ydata', [yvecA(i) yvecA(np+i)],'Color','k','LineWidth',3,'marker','o');
     qp=quiver(xA(i),yA(i),axA(i)*sca,ayA(i)*sca,0);  % Trasare vA
     set(qp,'autoscalefactor',sca,'MaxHeadSize',0.1,'Color','r','LineWidth',3);
     aAs=sprintf('%4.1f', aA(i)); 
     text(xA(i)+sca/2*axA(i),yA(i)+sca/2*ayA(i), aAs ,'Color','b','fontsize', 15);  % Marcare a
     plot(0,-0.001,'marker','^','markerfacecolor','k','markersize',12);             % desen articulatie manivela
     
     %Trasare culisa 
      patch([xB(i)-dx,xB(i)+dx,xB(i)+dx,xB(i)-dx],[yB(i)-dy,yB(i)-dy,yB(i)+dy,yB(i)+dy],'g');
    
    % Trasare biela
     line ('xdata', [xvecB(i) xvecB(np+i)] , 'ydata', [yvecB(i) yvecB(np+i)],'Color','B','LineWidth',3,'marker','o');
     qp=quiver(xB(i),yB(i),axB(i)*sca,ayB(i)*sca,0);  % Trasare vB
     set(qp,'autoscalefactor',scv,'MaxHeadSize',0.1,'Color','r','LineWidth',3);
     aBs=sprintf('%4.1f', aB(i)); 
     text(xB(i)+sca/2*axB(i),yB(i)+sca/2*ayA(i), aBs ,'Color','b','fontsize', 15);  % Marcare a
    
     % Pregatire grafic noua pozitie
     pause(0.1);            % Pauza intre secvente
     if (i<np)
        cla;                       % Sterge imaginea pentru noul desen
     end
     
 end

 % ------ Se scriu titlul, axele, grid, etc. --------
 set(gca, 'linewidth', 2, 'fontsize', 15);
 hold off;  axis('equal');grid;
 title('Pozitii mecanism & Vectori acceleratie')
 xlabel('x','FontSize',15); ylabel('y','FontSize',15)
 pause(2);
 
 
 

%======== Listarea rezultatelor ===========================
 for i=1:np
   fprintf('t=%5.2f [s] teta=%4.1f [grade] xA=%5.3f [m] yA=%5.3f [m/s] xB=%5.3f [m] vA=%5.3f [m/s] vB=%5.3f [m/s]  \n',t(i),180/pi*theta(i),xA(i),yA(i),xB(i),vA(i),vB(i));
 end