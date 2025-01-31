clc
clear
close all


seq = zeros(1, 1000);
seq(1) = 1;

n = 1;

while true
    s = 1;
    for i=1:n
        s = s+seq(i)^2;
    end
    seq(n+1) = s/n;
    if seq(n + 1) ~= floor(seq(n + 1))
        disp(['First non-integer position: ', num2str(n + 1)])
        break
    end
    n = n+1;
end