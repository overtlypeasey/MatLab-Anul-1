%% PROGRAM FOR DETERMINING THE RESULTANT FORCE AND RESULTANT TORQUE
%% OF A SYSTEM OF FORCES AND MOMENTS ACTING ON A RIGID BODY

clear all; close all; clc;

%% Reading the edges of the paralelipiped

O1_O2=input('> Please write the length of edge O1_O2: ');
O1_O4=input('> Please write the length of edge O1_O4: ');
O1_O5=input('> Please write the length of edge O1_O5: ');
X= [0,O1_O2,O1_O2,0,0,O1_O2,O1_O2,0];
Y= [0,0,O1_O4,O1_O4,0,0,O1_O4,O1_O4]; 
Z= [0,0,0,0,O1_O5,O1_O5,O1_O5,O1_O5];

%% Reduction set init.

Rx=0;  Ry=0;  Rz=0;

M0x=0; M0y=0; M0z=0;

%% Forces

n=input('> Please write number of forces: ');
for i=1:n
    F = input(['  The magnitude of force ',num2str(i),' is: ']);
    P = input(['  Force ',num2str(i),' starts from point O_: ']);
    Q = input(['  Force ',num2str(i),' ends in point O_: ']);

    B=F/sqrt((X(Q)-X(P))^2+(Y(Q)-Y(P))^2+(Z(Q)-Z(P))^2);
    Fx=B*(X(Q)-X(P));
    Fy=B*(Y(Q)-Y(P));
    Fz=B*(Z(Q)-Z(P));

    Mx=Y(P).*Fz-Z(P).*Fy;
    My=Z(P).*Fx-X(P).*Fz;
    Mz=X(P).*Fy-Y(P).*Fx;

    Rx=Rx+Fx;
    Ry=Ry+Fy;
    Rz=Rz+Fz;

    M0x=M0x+Mx;
    M0y=M0y+My;
    M0z=M0z+Mz;
    
end

%% Moments
m=input('> Please write number of moments: ');
for j=1:m
    M = input(['  The magnitude of torque ',num2str(j),' is: ']);
    P = input(['  Moment ',num2str(j),' starts from point O_: ']);
    Q = input(['  Moment ',num2str(j),' ends in point O_: ']);
    
    C=M/sqrt((X(Q)-X(P))^2+(Y(Q)-Y(P))^2+(Z(Q)-Z(P))^2);
    M2x=C*(X(Q)-X(P));
    M2y=C*(Y(Q)-Y(P));
    M2z=C*(Z(Q)-Z(P));

    M0x=M0x+M2x;
    M0y=M0y+M2y;
    M0z=M0z+M2z;

end



%% Solution for the reduction set

fprintf(['R =(',num2str(Rx),',',num2str(Ry),',',num2str(Rz),');']);
fprintf(['M =(',num2str(M0x),',',num2str(M0y),',',num2str(M0z),');']);3
