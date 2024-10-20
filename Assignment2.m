theta = linspace(0, 2*pi);
zx0 = -0.1;
zxy0 = 0.1;
r = 1.1;
z0 = zx0 + zxy0*i;
z = z0 + r*exp(i*theta);
airfoil = z + 1./z; %Joukowsky transform

figure      %Creates a figure window w/o plot
%plot(real(z), imag(z))
plot(real(airfoil), imag(airfoil))
axis('equal')
grid on
xlabel(x)
ylabel(y)
title("Joukowsky Airfoil")


hold on
plot(-1, 0, 'gd')
plot(1, 0, 'rd')
hold off