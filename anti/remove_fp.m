function [im_new,error] =  remove_fp(error_quan,dir_flag,im,pre_im_comp,tau,win)
error=1;
addpath('./codec');
%for test
%lambda1       = 1/3.48;                       %up to change
%lambda2       = 1/3.48;
[w, h]       = size(im);
step       = tau*2+1;
dir_flag_col = dir_flag';
dir_flag_col = dir_flag_col(:);

%%%%%estimate parameter lambda
error_quan1 = error_quan(3:w-4, 4:h-4); %original is error_quan(3:508, 4:508)
error_quan_col   = error_quan1(:);
lambda1      = para_estimate(error_quan_col, tau);
lambda2      = lambda1;

y = exp(-lambda1*step/2);
pro_zero = lambda1/(2*(1-y))*exp(-lambda1*abs(-tau:tau));   %the probability of [-tau tau]
pro_szero = lambda1/(y^(-1)-y)*exp(lambda1*(-tau:tau));

%%%%%revise the lambda and generate the dither variables
error_quan_col   = error_quan(:);
[new_lambda, D_cell] = revisePara(error_quan_col, lambda1, tau);
[v_N,v_q] = hist(error_quan_col,unique(error_quan_col)); %get the number of each qk

% numsamp_zeros    = length(find(error_quan_col==0));
% numsamp_nonzeros = length(error_quan_col)-numsamp_zeros;
% numsamp_b_zeros =  length(find(error_quan_col>0));
% numsamp_s_zeros =  length(find(error_quan_col<0));
% 
% D_zero_value     = rev_ran_generator(lambda1,step,numsamp_zeros,0);   
% D_nonzero_value  = rev_ran_generator(lambda2,step,numsamp_nonzeros,1);
% % save 'D_zero_value.mat'  D_zero_value;
% % save 'D_nonzero_value.mat' D_nonzero_value;
% % load('D_zero_value.mat'); load('D_nonzero_value.mat');
% D_b_zero_value_temp = D_nonzero_value(1:numsamp_b_zeros);
% D_s_zero_value_temp = D_nonzero_value(numsamp_b_zeros+1:length(D_nonzero_value));
% 
% D_b_zero_value = round(D_b_zero_value_temp);
% D_s_zero_value = round(-D_s_zero_value_temp);
% D_zero_value   = round(D_zero_value);


%%%GAP predict
c_bs_tau = 0;
c_cur_quan_err = [];
im_comp       =     double(im);
for i = win(1) : w-win(1)-1      %这里还是考虑到边界的问题。当前像素的L个领域必须都是计算过边界的。行数不用加，但是列数需要
    fprintf('current row is %d\n',i);
    for j = win(2)+1: h-win(2)-1
            A                   = getTemplate(im_comp, win, i,j);
            P_id                = A*dir_flag_col;                 %get the L prediction directions
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%calculate the interval of current dither
            cur_quan_err   = error_quan(i,j);
            old_cur_pre    = pre_im_comp(i,j);
            [new_cur_pre,~]    = get_pre_value(im_comp,i,j);          %get the new current predition value
            c              = round(new_cur_pre) - old_cur_pre;    %the prediction value may be float.c的值可能超过tau
            if c > tau | c < -tau
                c_bs_tau = c_bs_tau+1;
                c_cur_quan_err = [c_cur_quan_err cur_quan_err];
                continue;
            end
            sol      = solve_GAP_theta(im_comp,P_id,cur_quan_err,tau,i,j,pro_zero, pro_szero); %sol_can_1保证不改变方向，而且误差在-tau到tau
            %sol = [-tau tau];
            
%             if(cur_quan_err == 0)                                       %quantized error is zero or nonzero
%                 D_vec = D_zero_value;
%             elseif(cur_quan_err > 0)
%                 D_vec = D_b_zero_value;
%             else 
%                 D_vec = D_s_zero_value;
%             end
            temp_indx = find(v_q == cur_quan_err); 
            D_vec     = D_cell{temp_indx};
            
            
            [d, indx]  =   get_dither(sol,D_vec,cur_quan_err,tau,c);
            
            D_vec(indx) =  [];                                        %remove the used dither
           % D_vec= [D_vec; d];                                         %把它甩到向量末端，本来可以直接删除，但是为了防止直接删除可能导致最后为空，干脆不删
            
           
%             if(cur_quan_err == 0)
%                 D_zero_value    =  D_vec;                           
%             elseif(cur_quan_err > 0)
%                 D_b_zero_value  =  D_vec;
%             else
%                 D_s_zero_value  =  D_vec;
%             end
            D_cell{temp_indx} = D_vec;
           %% 下面这段代码注释掉后 程序应该不会有任何错误，并且分布和压缩后分布一样，用来测试程序。       
            im_comp(i,j) = im_comp(i,j) + d;                          %change current pixel value
          
          
    end
end
            im_new = im_comp;
            c_bs_tau
end

function [d, indx] = get_dither(cand_sol,D_vec,cur_quan_err, tau,c)

    n = length(D_vec);
   %     ran = rand(1,1)/8;
%     len = round(n * ran)+1;
    
    
%     if n == len
%         fprintf('D_vec is empty');
%     end
%len=10;
len= min(10,length(D_vec));
    for indx = 1:len
%         temp = D_vec(indx);
%         if cur_quan_err==0
%             d = round(temp);
%         elseif cur_quan_err>0
%             d = round(temp - 0.5) - tau;
%         elseif cur_quan_err < 0
%             d = round(-temp + 0.5) + tau;
%         end
        d = D_vec(indx)+c;
        if(any(cand_sol == d))
            break;
        end
    end
    
    if indx == len
%         fprintf('no solution found in D_vec! the candiate solution should be %d', cand_sol);
%         error('\n');
        %%only test
        indx = max(1,round(len*rand(1,1)));  %randomly choose a pixel 
        d = D_vec(indx)+c;
       % fprintf('no solution found in D_vec! the candiate solution sholud be %d', cand_sol );
        
    end
end  




