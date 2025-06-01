% Program for calculating and tracing the trajectory of the quadrilateral mechanism  
% Prof. M.V. Predoi, Departamentul de Mecanica, UPB.
% Se pot defini geometria mecanismului si rotatia manivelei:
%    theta(t) 
% Mecanismul este R-L1-L2 cu proiectie D pe Ox.
% ====== Setari ale ferestri grafice ====================
close all; clc;

%========Datele problemei ==========================
R=0.8;       % Crank radius (bar 1) [m]
L1= 1.7;    % Length Bar 2[m]
L2= 2;    % Length Bar 3 [m]
D= 2.2;       % Distance between fixed joints O1 and O3

if (R+D > L1+L2)
  disp('EROARE: Distanta prea mare intre articulatiile fixe !')
  return
 end
 
%---------------------------------------------------
% Daca miscarea are acceleratie:
eps=35;              % Angular acceleration           [rad/s^2]
%----------------------------------------------------
theta0=pi/7;       % Initial angular position    [rad]
%---------------------------------------------------
omega0=40;    % Initial angular velocity  [rad/s]
%---------------------------------------------------

% Se defineste un set de valori pentru timp: 
 t = linspace(0,0.2*pi/omega0,49);       % t0 : time step : final time
 
  unu=ones(size(t));
  np=length(t);
 %---------------------------------------------------
 % Definirea ecuatiilor de miscare:
 theta = theta0 + omega0*t + eps/2*t.^2 ;     % theta(t)  [rad]
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
 
% ============================================================
% === Rezolvare forma geometrica patrulater =======================
% Notatii: alfa = unghi dintre Ox si L1
%                beta =  unghi dintre Ox si L2
 if( any (D^2 - 2*D*R*cos(theta) + R^2 <0))
   disp('Eroare in datele de intrare')
   return
 end
 
 cosafetbet=(  (D^2 +R^2 -L1^2 -L2^2)*unu - 2*R*D*cos(theta)  )/2/L1/L2;
 % Verificare 
 if (abs(cosafetbet)>1)
   disp('Imposibil geometric; |cos (af+bet)|>1');
   return
 end
   
 fi=real(acos(cosafetbet));
 fideg=fi*180/pi;
 aa=L2+L1*cos(fi);                    %  poate fi: L2-L1*cos(fi)
 bb=L1*sin(fi);                            % poate fi:  -L1*sin(fi)
 cc=R*sin(theta);

 for i=1:np
   
     alfac(1)=atan(real(  (aa(i)*bb(i)  + cc(i)*sqrt(aa(i)^2+bb(i)^2-cc(i)^2)) ./ (aa(i)^2-cc(i)^2)));
     alfac(2)=atan(real(  (aa(i)*bb(i)  - cc(i)*sqrt(aa(i)^2+bb(i)^2-cc(i)^2)) ./ (aa(i)^2-cc(i)^2)));
     alfac(3)=atan(real( - (aa(i)*bb(i)  + cc(i)*sqrt(aa(i)^2+bb(i)^2-cc(i)^2)) ./ (aa(i)^2-cc(i)^2)));
     alfac(4)=atan(real( - (aa(i)*bb(i)  - cc(i)*sqrt(aa(i)^2+bb(i)^2-cc(i)^2)) ./ (aa(i)^2-cc(i)^2)));
     alfac=sort(alfac,2,'descend');
    
   trigangle= rem (theta(i),2*pi);
   if trigangle<pi/2
      alfa(i)=alfac(2);    
   elseif (trigangle>=pi/2 && trigangle<pi)
      alfa(i)=alfac(2); 
   elseif (trigangle>=pi && trigangle<3/2*pi)
      alfa(i)=alfac(1); 
   elseif (trigangle>=3/2*pi && trigangle<=2*pi)
       % Verificare solutie: daca bara L1 trece de verticala:
        alfa(i)=alfac(1);    
        xB=R*cos(theta(i)) + L1*cos(alfa(i));                      % Proiectie B pe Ox
        yB=R*sin(theta(i))  + L1*sin(alfa(i));                        % Proiectie B pe Oy
        L2c=sqrt((xB-D).^2+yB.^2);
        if (abs(L2c-L2)  > 0.001*L2)
           fprintf('theta=%g (deg) & alfa1=%g (deg) : L2= %g  !! \n',theta(i)*180/pi, alfa(i)*180/pi, L2c);  
           alfa(i)=pi+alfac(4);
           xB=R*cos(theta(i)) + L1*cos(alfa(i));                      % Proiectie B pe Ox
           yB=R*sin(theta(i))  + L1*sin(alfa(i));                        % Proiectie B pe Oy
           L2c=sqrt((xB-D).^2+yB.^2);
           fprintf(' ales suplementul: alfa=%g (deg) => L2= %g. \n',alfa(i)*180/pi, L2c);  
        end
   end
 end
 
 beta=alfa-fi;                        % unghiul barei L2 cu orizontala
 
 %----------------------------------------------------
 % Calcul pozitii bare

 xB=R*cos(theta)+L1*cos(alfa);                      % Proiectie B pe Ox
 yB=R*sin(theta)+L1*sin(alfa);                        % Proiectie B pe Oy
 
 OI= D*sin(abs(beta))./sin(theta+abs(beta));
 O3I=OI.*sin(theta)./sin(abs(beta));
 xI=OI.*cos(theta);
 yI=OI.*sin(theta);
 omega2=-vn./(OI-R);
 vnB=omega2.*(O3I-L2);
 vxB=-vnB.*sin(beta);
 vyB=vnB.*cos(beta);
 vB=sqrt(vxB.^2 + vyB.^2);                                % Modulul vitezei in B
       
 
 %====== Trasarea miscarii mecanismului si vectori viteza ======
 figure(1);

 indv=fix(linspace(1,np,11));          % Se vor trasa 11 viteze
 O1=zeros(size(xA)); 
 O3x=D*ones(size(t));
 O3y=zeros(size(xA)); 
 xvec1=[O1 xA];
 yvec1=[O1 yA];
 xvec2=[xA xB];
 yvec2=[yA yB]; 
 xvec3=[xB O3x];
 yvec3=[yB O3y]; 
  
  
 maxv=max(sqrt(vr.^2+vn.^2));
 scv=0.5/maxv;                                  % Scara trasarii vitezelor
 dx=scv;
 dy=0.7*scv;
 axis([-2*R    D+L2   -(D+L2+R)/2    (D+L2+R)/2]);
 hold on
 
 for i=1:np
     line('xdata', [0 D], 'ydata', [0 0]);  % Trasare directie culisa (Ox)
     
     % Trasare bara motoare de lungime R:
     line ('xdata', [xvec1(i) xvec1(np+i)] , 'ydata', [yvec1(i) yvec1(np+i)],'Color','k','LineWidth',3,'marker','o');
     qp=quiver(xA(i),yA(i),vxA(i)*scv,vyA(i)*scv,0);  % Trasare vector vA
     set(qp,'autoscalefactor',scv,'MaxHeadSize',0.1,'Color','r','LineWidth',3);
     
     plot(xI(i),yI(i),'*r');            % Marcare CIR ca * rosie
     text(xI(i),yI(i),'CIR');          % Notatie CIR
     
     vAs=sprintf('%4.1f', vA(i)); 
     text(xA(i)+scv/2*vxA(i),yA(i)+scv/2*vyA(i), vAs ,'Color','b','fontsize', 15);  % Marcare valoare vA
     plot(0,-0.001,'marker','^','markerfacecolor','k','markersize',12);             % desen articulatie O
     text(scv/4,scv/4,'O1','FontSize',16);                                                            % Notatie O1 
     plot(D,-0.001,'marker','^','markerfacecolor','k','markersize',12);             % desen articulatie O3
     text(D+scv/4,scv/4,'O3','FontSize',16);                                                       % Notatie O3
     
     %Trasare bara L1
     line ('xdata', [xvec2(i) xvec2(np+i)] , 'ydata', [yvec2(i) yvec2(np+i)],'Color','k','LineWidth',3,'marker','o');
     qp=quiver(xB(i),yB(i),vxB(i)*scv,vyB(i)*scv,0);  % Trasare vB
     set(qp,'autoscalefactor',scv,'MaxHeadSize',0.1,'Color','r','LineWidth',3);
          
     %Trasare bara L2
     line ('xdata', [xvec3(i) xvec3(np+i)] , 'ydata', [yvec3(i) yvec3(np+i)],'Color','k','LineWidth',3,'marker','o');
     
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

 
%======== Listarea rezultatelor ===========================
 for i=1:np
   fprintf('t=%5.2f [s] teta=%4.1f [%c] A{x,y}={%5.3f,%5.3f} B{x,y}={%5.3f,%5.3f}  \n',t(i),180/pi*theta(i),char(15),xA(i),yA(i),xB(i),yB(i));
   fprintf('           vA=%5.3f [m/s] vB=%5.3f [m/s]  \n',vA(i),vB(i));
 end
%=============================================================================================
 
 