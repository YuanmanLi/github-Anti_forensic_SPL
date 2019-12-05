
function [v_lambda_new, D_cell] = revisePara(error_quan_col, lambda, tau)
    quantized_step = 2*tau + 1;
    errorVector = error_quan_col;
    A = errorVector;
    [v_N1,v_q1] = hist(A,unique(A)); %get the number of each qk
    N = sum(v_N1);
    center_indx = find(v_N1./N > 6*1.0e-3);
    if tau >= 7
         center_indx = find(v_N1./N > 10*1.0e-3);
    end
    
    center_indx_ex = [center_indx(1)-1, center_indx, center_indx(end)+1];
    v_N = v_N1(center_indx_ex);
    v_q = v_q1(center_indx_ex);
    
    v_errortable = [v_q(1)-tau:v_q(end)+tau];
    
   % v_lambda   = repmat(lambda, 1,length(center_indx));      % default lambdas are all estimated one
    
    %%%%%%%%set the lambdas of the boudary quantization bins be low %%
    %%%%%%%%otherwise sometimes, it's unable to revise
    lambda_scale = 0.6;
     if tau >= 7
        lambda_scale = 0.4;
    end
    
    n_q = length(v_N1);
    v_lambda_new = repmat(lambda, 1,n_q);
    v_lambda_new = v_lambda_new * lambda_scale;      
    v_lambda_new (center_indx) =  repmat(lambda, 1,length(center_indx));
    
    
    
    v_lambda1 = repmat(lambda, 1, length(center_indx)+2);
    v_lambda1(1) = v_lambda_new(center_indx(1)-1);   % most left boudary taken in to consideration
    v_lambda1(end) = v_lambda_new(center_indx(end)+1);   % most right boudary taken in to consideration
    
    v_lambda2 =  repmat(v_lambda1, quantized_step, 1);
    v_lambda2 = v_lambda2(:);
    
    f          = @(x,lambda)lambda/2.*exp(-lambda.*abs(x));
    v_probability = f(v_errortable, v_lambda2');
    v_probability = reshape(v_probability, quantized_step, length(v_probability)/quantized_step);
    v_probability = sum(v_probability);   %for test lana256, quantized bin = 0, probability = 0.7286 with the estimated lambda 
    %%%%calculate the U
    U = v_N./(N.*v_probability);        % quantized bin = 0, v_N/N = 0.6848. so the real probability = 0.6848
    
    zero_indx   = find(v_q == 0);
    v_N_neg     = v_N(1:zero_indx);
    v_N_pos     = v_N(zero_indx:end);
    v_q_neg     = v_q(1:zero_indx);
    v_q_pos     =  v_q(zero_indx : end);
    
    v_boundary_pre       = [v_q_neg-tau;v_q_pos+tau];
    v_boundary_pre_pre   = [v_q_neg-tau+1;v_q_pos+tau-1];
    v_boundary_next      = [v_q_neg-tau-1;v_q_pos+tau+1];
     
  
    alpha1 = 0.9;
    alpha2 = 1.3;
    if tau >= 7
        alpha1 = 0.8;
        alpha2 = 1.1;
    end
    %%%%%handle the  case where qk <0
    for i = 2:1:zero_indx
       cur_q        =  v_q(i);
       lambda                 =   v_lambda1(i);      %update lambda
       lambda_pre             =   v_lambda1(i-1);  
       pre          = v_boundary_pre(i);
       prepre       = v_boundary_pre_pre(i);
       next         = v_boundary_next(i);
       boundary_pre_p         = f(pre, lambda) * U(i);
       boundary_pre_pre_p     = f(prepre, lambda)* U(i);
       boundary_next_p        = f(next, lambda_pre)* U(i-1);
       v_dif_pre_p            = boundary_pre_pre_p - boundary_pre_p;
       v_dif_p                = boundary_pre_p - boundary_next_p;
       if(v_dif_p < v_dif_pre_p * alpha1 | v_dif_p > v_dif_pre_p * alpha2)
            opts = optimset('LargeScale','off','display','on','Algorithm','active-set');
            [x,fval,exitflag]=fmincon('obj_fun',lambda,[],[],[],[],lambda*0.01,lambda*100,'con_fun',opts, lambda,pre, prepre,boundary_next_p, quantized_step, v_N(i), N, alpha1, alpha2);
            fprintf('need to revise at pos %d, ori = %f, revised = %f\n', v_q(i), v_lambda1(i), x);
             if (lambda*0.005> x | x> lambda*10 )
                 fprintf('cancell the modification\n');
                 continue;
            end
             v_lambda1(i) = x;
             v_q_cur             =  [cur_q - tau : cur_q + tau];                      %upadate U(i-1)
             proba_cur           =   f(v_q_cur, x); 
             U(i)              =   v_N(i)./(N.*sum(proba_cur));      
       end
    end
  clear lambda_pre;
     %%%%%handle the  case where qk >0
    for i = length(v_lambda1)-1:-1:zero_indx    % lambda of bin 0 has already be revised by the case of negtive. length(v_lambda1) is 1 bigger than that of v_boundary_pre(-5.5)
       cur_q        =  v_q(i);
       lambda_next                 =   v_lambda1(i+1);      %update lambda
       lambda             =   v_lambda1(i);  
       pre          = v_boundary_pre(i+1);
       prepre       = v_boundary_pre_pre(i+1);
       next         = v_boundary_next(i+1);
       boundary_pre_p         = f(pre, lambda) * U(i);
       boundary_pre_pre_p     = f(prepre, lambda)* U(i);
       boundary_next_p        = f(next, lambda_next)* U(i+1);
       v_dif_pre_p            = boundary_pre_pre_p - boundary_pre_p;
       v_dif_p                = boundary_pre_p - boundary_next_p;
       if(v_dif_p < v_dif_pre_p * alpha1 | v_dif_p > v_dif_pre_p * alpha2)
            opts = optimset('LargeScale','off','display','off','Algorithm','active-set');
            [x,fval,exitflag]=fmincon('obj_fun',lambda,[],[],[],[],lambda*0.01,lambda*100,'con_fun',opts, lambda,pre, prepre,boundary_next_p, quantized_step, v_N(i), N, alpha1, alpha2);
            fprintf('need to revise at pos %d, ori = %f, revised = %f\n', cur_q, v_lambda1(i), x);
            if (lambda*0.005 > x | x> lambda*10 )
                 fprintf('cancell the modification\n');
                 continue;
            end
             if(i == zero_indx)
                x = min(v_lambda1(i), x);  %get the smaller lambda at bin0 
             end
             v_lambda1(i)      =   x;   %update lambda
             v_q_cur             =  [cur_q - tau : cur_q + tau];                      %upadate U(i-1)
             proba_cur           =   f(v_q_cur, x); 
             U(i)              =   v_N(i)./(N.*sum(proba_cur));         
       end
    end
   
  
   v_lambda_new(center_indx_ex) = v_lambda1;
   D_cell = cell(1, n_q);
   maxval       = tau*2+1;
    for i = 1 : n_q
        cur_lambda = v_lambda_new(i);
        cur_num    = v_N1(i);
        cur_q      = v_q1(i);
        if(cur_q < 0)
           cur_d        =  rev_ran_generator(cur_lambda,maxval,cur_num,1); 
           D_cell{i}    =  round(-cur_d);
        elseif cur_q ==0
           cur_d        =  rev_ran_generator(cur_lambda,maxval,cur_num,0); 
           D_cell{i}    =  round(cur_d);
        else
           cur_d        =  rev_ran_generator(cur_lambda,maxval,cur_num,1); 
           D_cell{i}    =  round(cur_d);  
        end  
    end
    
%     for i = 1 : n_q
%         cur_q      = v_q1(i);
%         temp_idx   = find(errorVector == cur_q);
%         errorVector(temp_idx) = errorVector(temp_idx) +  D_cell{i};
%     end
%     
%    dis_low = -80; dis_up = 80;
% %draw_error_distribution(err_im_ori, 'ori error', dis_low, dis_up);
% %draw_error_distribution(err_im_comp, 'calic error', dis_low, dis_up);
% draw_error_distribution(errorVector, 'new error', dis_low, dis_up);
%    
% end