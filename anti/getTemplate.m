%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Created by Li Yuanman 
%% Jan. 16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A]  = getTemplate(img, win,i,j)
[a, b] = size(img);
%%%%%%%%only for GAP predictor
L      = 7;
w      = win(1):a-win(1)+1;
h      = win(2):b-win(2)+1;
A      = zeros(L, a*b);
i = i-1;
A(1,i*a+j+1) 	  = 1; %e
A(2,i*a+j+2)      = 1; %ee
A(3,(i+1)*a+j-1)  = 1; %sw
A(4,(i+1)*a+j)    = 1; %s
A(5,(i+1)*a+j+1)  = 1; %se
A(6,(i+2)*a+j-1)  = 1; %ssw
A(7,(i+2)*a+j)    = 1; %ss



