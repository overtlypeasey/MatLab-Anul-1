A = 3;
f = 1;
om = 2*pi*f;
t = linspace(0,3);
fi = 0;
zg = randn(size(t)) * 0.5;

signal = A*sin(om*t+fi) + zg;

figure
plot(t, signal)

save datesin