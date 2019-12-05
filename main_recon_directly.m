%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Created by Li Yuanman 
%% Jan. 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
dbstop error
addpath('codec');
addpath('anti');
addpath('BM3D');
DataBase_4 = '..\DataBase\'; %28 samples
DataBase_Cur = DataBase_4;
test_n = 28;
for imgidx =29:29 

tolerance = 5;

img_format          = strcat(DataBase_Cur,'Satellite',int2str(imgidx), '.raw'); 

par.win             =  [3, 3];
par.tau             =   tolerance;
par.I               =    double( readraw( img_format ) );

[pre_im_comp1, err_im_comp1, par.im_comp]       =  GAP_coding(par.I, par.tau, par.win);


%%use GAP predict the compressed image
[pre_im_ori, err_im_ori,~]              = GAP_predict(par.I);  %pre_im_ori: prediction value of original image
[pre_im_comp, err_im_comp,dir_flag]     = GAP_predict(par.im_comp);     % the err_im_comp should be equal to err_im_comp1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% only for test. the perfect result, i just add the dither to predition
% error directly
error_quan_col   = err_im_comp1(:);
maxval       = par.tau*2+1;
lambda1      = para_estimate(error_quan_col, par.tau);

[new_lambda, D_cell] = revisePara(error_quan_col, lambda1,  par.tau);
[v_N,v_q] = hist(error_quan_col,unique(error_quan_col)); %get the number of each qk
 errorVector = error_quan_col;
for i = 1 : length(v_N)
    cur_q      = v_q(i);
    temp_idx   = find(errorVector == cur_q);
    errorVector(temp_idx) = errorVector(temp_idx) +  D_cell{i};
end
re_err_mat   = reshape(errorVector',size(par.im_comp)); 

[~,~,recon]  = GAP_Decoding(  par.im_comp,re_err_mat, par.tau,[3,3]);
%[~,~,recon]  = GAP_Decoding(  par.I,err_im_comp1, par.tau,[3,3]);

dis_low = -80; dis_up = 80;
%draw_error_distribution(err_im_ori, 'ori error', dis_low, dis_up);
draw_error_distribution(err_im_comp, 'calic error', dis_low, dis_up);
draw_error_distribution(errorVector, 'new error', dis_low, dis_up);
figure;
imshow(uint8(recon));
%%%%%%%%%%%%%%%%%%%
%only for test.  证明直接对压缩的图像用bm3d去噪后， 是不能抹去压缩痕迹的
% dis_low = -50; dis_up = 50;
% y =  im2double(uint8(par.I));
% z =  im2double(uint8(par.im_comp));
% [~,soft_removed_im] = BM3D(y, z,2.8);
% soft_removed_im  = im2uint8(soft_removed_im);
% [pre_im_removed, err_im_removed,~]      = GAP_predict(soft_removed_im);
% draw_error_distribution(err_im_removed, 'BM3D footprint removed error3.0', dis_low, dis_up);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%remove the footprint
% [par.im_removed_pf,err_im_removed_pf]   =  remove_fp(err_im_comp1,dir_flag,par.im_comp,pre_im_comp,par.tau,par.win);
% %%%%%%%%%%%%%%%%%%%%
% dis_low = -100; dis_up = 100;
% 
% [pre_im_removed, err_im_removed,~]      = GAP_predict(par.im_removed_pf);
% %imshow(uint8(pre_im));
% % draw_error_distribution(err_im_ori, 'prediction error', dis_low, dis_up);
%  draw_error_distribution(err_im_comp, sprintf('quantized prediciton error (step = %d)', 2*tolerance+1), dis_low, dis_up);
% 
% draw_error_distribution(err_im_removed, sprintf('new prediciton error (step = %d)', 2*tolerance+1), dis_low, dis_up);
% 
% psrn_com_re = psnrfun(par.im_comp,par.im_removed_pf,[0 0])
% psrn_re_I   = psnrfun(par.I,par.im_removed_pf,[0 0])
% psrn_com_I  = psnrfun(par.im_comp,par.I,[0 0])
% maxerror    = max(max(abs(par.im_comp-par.im_removed_pf)))
% 
% temp = par.im_removed_pf;
% save im_removed_pf.mat temp
% figure;
% imshow(uint8(par.im_comp));
% 
% y =  im2double(uint8(par.I));
% z =  im2double(uint8(par.im_removed_pf));
% [~,soft_removed_im] = BM3D(y, z,2.8);
% soft_removed_im  = im2uint8(soft_removed_im);
% [pre_im_removed, err_im_removed,~]      = GAP_predict(soft_removed_im);
%draw_error_distribution(err_im_removed, 'revised footprint removed error1 0.3', dis_low, dis_up);

end
