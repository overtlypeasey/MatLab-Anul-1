clear
clc
close all

fileID = fopen("The giant panda.docx");
text = fscanf(fileID, '%c');
fclose(fileID);

up = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
lo = 'abcdefghijklmnopqrstuvwxyz';
sp = '/=-+_)(&^%$#@!,.><;"|\';

for i = 1:length(text)
    for k = 1:length(lo)
        if data(i) == up(k)
            data(i) = lo(k);
        end
    end
end
word = [];
for i = 1:length(text)
    for k = 1:length(lo)
        if data(i) == lo(k)
            word = [word data(i)];
            pause
        end
    end
end