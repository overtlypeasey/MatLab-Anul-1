clear
clc
close all

fs = 44100;

% eighth note as 0.5 sec
lo_so = key(47-12, 8, fs);
lo_la = key(49-12, 8, fs);
lo_si = key(51-12, 8, fs);
do = key(52-12,8,fs);  
re = key(54-12,8,fs);  
mi = key(56-12,8,fs);  
fa = key(57-12,8,fs);  
so = key(59-12,8,fs);  
la = key(61-12,8,fs);  
si = key(63-12,8,fs);  

% quarter note as 1 sec
lo_so_4 = key(47-12, 4, fs);
lo_la_4 = key(49-12, 4, fs);
lo_si_4 = key(51-12, 4, fs);
do_4 = key(52-12,4,fs);  
re_4 = key(54-12,4,fs);  
mi_4 = key(56-12,4,fs);  
fa_4 = key(57-12,4,fs); 
so_4 = key(59-12,4,fs);  
la_4 = key(61-12,4,fs);  
si_4 = key(63-12,4,fs);  

lo_so_3 = key(47-12, 3, fs);
lo_la_3 = key(49-12, 3, fs);
lo_si_3 = key(51-12, 3, fs);
do_3 = key(52-12,3,fs);  
re_3 = key(54-12,3,fs); 
mi_3 = key(56-12,3,fs); 
fa_3 = key(57-12,3,fs);  
so_3 = key(59-12,3,fs);  
la_3 = key(61-12,3,fs); 
si_3 = key(63-12,3,fs); 

% make a song
line1 = [ re_4 re re_4 re ];
line2 = [ so_4 so so_4 so ];
line3 = [ la_4 la la_4 la ];
line4 = [ si_3 so_3 ];
line5 = [ re_4 re re_4 re];
line6 = [ mi_4 mi mi_4 mi];
line7 = [ re_4 do lo_si_4 lo_la ];
line8 = [ lo_so_3 lo_so_3 ];

song = [line1 line2 line3 line4 line5 line6 line7 line8];
% plot(song)
sound(song,fs,24);

function wave = key(p, n, fs)
    t = 0:1/fs:4/n;
    idx = 440*2^((p-49)/12);
    
%     method 1 - orginal
    % wave = (sin(2*pi*idx*t));

%     method 2 - exponential decreasing 
    tt = 4/n:-1/fs:0;
    wave = (sin(2*pi*idx*t)).*exp(tt);
    wave = wave./max(wave);
    
%     method 3 - triangle decreasing 
    % mid = (t(1)+t(end))/2;
    % tri = -(abs(t-mid)-mid);
    % tri = tri./max(tri);
    % wave = (sin(2*pi*idx*t)).*tri;
end