%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Created by Li Yuanman 
%% Jan. 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pre_value,flag] = get_pre_value(im_comp,i,j)
            flag            =  -1;
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
                flag      = 1;
            elseif (dv_dh < -80)
                pre_value = n;
                flag      = 2;
            else
                    hat_value =(w+n)/2 + (ne-nw)/4;
                    if (dv_dh > 32)
                        pre_value = (hat_value + w)/2;
                        flag      = 3;
                    elseif (dv_dh > 8)
                        pre_value = (3*hat_value + w)/4;
                        flag      = 4;
                    elseif (dv_dh < -32)
                        pre_value = (hat_value + n)/2;
                        flag      = 5;
                    elseif (dv_dh < -8)
                        pre_value = (3*hat_value + n)/4;
                        flag      = 6;
                    elseif (dv_dh>=-8 && dv_dh <= 8)
                        pre_value = hat_value;     
                        flag      = 7;
                    end    
            end
            
%             if(pre_value == -1)
%                 fprintf('i=%d, j=%d, dv_dh=%d\n', i, j, dv_dh);
%                 error('no change for pre_value! GAP error!'); 
%             end      
            pre_value = max(0, pre_value);     % because hat_value may be negtive
end