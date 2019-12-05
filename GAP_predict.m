%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Created by Li Yuanman 
%% Jan. 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pre_im, err_im, dir_flags] = GAP_predict( im_comp, win)


if nargin < 2 || length(win) < 2
    win = [3, 3];
end
im_comp       =     double(im_comp);

[a b]               =     size(im_comp);
pre_im           =     zeros(a ,b);         %prediction value
err_im            =     zeros(a, b);         %error value
pre_im           =      double(pre_im);
err_im            =      double(err_im); 
dir_flags         =    zeros(a, b);

for i = win(1) : a-win(1)+1
    for j = win(2): b-win(2)+1
            flag                 = -1;
            [pre_value,flag]     =  get_pre_value(im_comp, i,j);
             pre_im(i,j)         = pre_value;
             dir_flags(i,j)      =  flag;
    end
end
pre_im       =  round(pre_im); 
error_pre   =  im_comp - pre_im;
err_im        = error_pre(win(1):a-win(1), win(2):b-win(2)); 

recover_im          = pre_im(3:b-3, 3:b-3) + err_im;
dif = double(recover_im) - double(im_comp(3:b-3, 3:b-3));


