%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Created by Li Yuanman 
%% Jan. 10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
function [pre, quan_err, restruction] = GAP_encoder(im_input, tau, winsize)

if nargin < 3 && length(winsize) < 2
    winsize = [3, 3];
end


[a, b]              =       size(im_input);
restruction     =        zeros(a, b);
pre                 =        zeros(a - winsize(1), b - winsize(2));
quan_err                  =        zeros(a - winsize(1), b - winsize(2));

[pre, quan_err, restruction]        =         GAP_coding( im_input, tau,winsize);








