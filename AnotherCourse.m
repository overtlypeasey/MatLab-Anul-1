clc
clear
close all

fhandle = fopen("Thegiantpanda.rtf", 'r');
k = 0;
while ~feof(fhandle)
    line = fgetl(fhandle);
    ii = find (line == ' ' | line == ',' | line == '.');
    found = 0;
    k = 0;
    dict = [];
    if ~isempty(ii)
        j = 1;
        for i=ii
            word = line(j:i);
            j = i+1;
            if ~isempty(word)
                disp(word)
                if k>0
                    found = 0;
                    for kk = 1:k
                        found = strcmp(dict(kk), word);
                    end
                end
                if found<1
                    k = k+1;
                    dict(k) = word;
                end
            else
                j = i+1;
            end
        end
    end
    % disp(line)
    % disp(ii)
    % break
end
fclose(fhandle);
