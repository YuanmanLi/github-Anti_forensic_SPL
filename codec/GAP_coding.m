%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Created by Li Yuanman 
%% Jan. 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pre_im, err_im, im_comp] = GAP_coding( im_comp, tau,win)


if nargin < 2 || length(win) < 2
    win = [3, 3];
end
im_comp       =     double(im_comp);
im_ori            =     im_comp;


[a b]               =     size(im_comp);
pre_im           =     zeros(a ,b);         %prediction value
err_im            =     zeros(a, b);         %error value
pre_im           =      double(pre_im);
err_im            =      double(err_im); 


for i = win(1) : a-win(1)+1
    for j = win(2): b-win(2)+1
            nn               =  im_comp(i-2, j);
            nne             =  im_comp(i-2, j+1);
            nw              =  im_comp(i-1,j-1);
            n                 =  im_comp(i-1, j);
            ne               =  im_comp(i-1, j+1);
            ww             =  im_comp(i, j-2);
            w                =  im_comp(i, j-1);
            
            dh               =  abs(w-ww) + abs(n-nw) + abs(n-ne);
            dv               =  abs(w-nw)  + abs(n-nn)  + abs(ne-nne);
            dv_dh         =  dv - dh;
            pre_value   =  -1;
            
          
            if (dv_dh > 80)
                pre_value = w;
            elseif (dv_dh < -80)
                pre_value = n;
            else
                    hat_value =(w+n)/2 + (ne-nw)/4;
                    if (dv_dh > 32)
                        pre_value = (hat_value + w)/2;
                    elseif (dv_dh > 8)
                        pre_value = (3*hat_value + w)/4;
                    elseif (dv_dh < -32)
                        pre_value = (hat_value + n)/2;
                    elseif (dv_dh < -8)
                        pre_value = (3*hat_value + n)/4;
                    elseif (dv_dh>=-8 && dv_dh <= 8)
                        pre_value = hat_value;     
                    end    
            end            
            
%             if(pre_value == -1)
%                 fprintf('i=%d, j=%d, dv_dh=%d\n', i, j, dv_dh);
%                 error('no change for pre_value! GAP error!'); 
%             end
            pre_value  = round(pre_value);
            pre_im(i,j) =  max(0, pre_value);     % because hat_value may be negtive
            
            curr_err     =   im_comp(i,j) - pre_im(i,j);
            quan_err    =  quan2(curr_err, tau);
            err_im(i,j)  =  quan_err;
            
            curr_reconstruction = quan_err + pre_im(i,j);
          %  curr_reconstruction = min(curr_reconstruction, 255);
         %    curr_reconstruction = max(curr_reconstruction, 0);
            im_comp(i,j)            = curr_reconstruction;
    end
end

psnr = psnrfun(im_comp, im_ori, win);
fprintf('compression success, tau = %d, pnsr = %.3f\n', tau, psnr);

