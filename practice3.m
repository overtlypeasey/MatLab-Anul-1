clear
clc 
close all

A = 20 * rand(7) - 10;
B = 20 * rand(7) - 10;

A(A<0) = 0;
B(B>0) = 0;

C = A*B - B*A;
C(C<0) = 0;
C(C>0) = 1;
C(1:end, 1) = 1;

% Convert binary columns to ASCII characters without using bi2de
dec_values = C' * 2.^(size(C,1)-1:-1:0)';
text = char(dec_values);
disp(text);