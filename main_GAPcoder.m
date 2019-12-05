%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Created by Li Yuanman 
%% Jan. 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
dbstop error
addpath('codec');
addpath('anti');
DataBase_1 = 'DataBase\'; %
DataBase_Cur = DataBase_1;
test_n = 36;
for imgidx =35:35

tolerance = 7;

img_format          = strcat(DataBase_Cur,'Satellite',int2str(imgidx), '.raw'); 

par.win             =  [3, 3];
par.tau             =   tolerance;
par.I               =    double( readraw( img_format ) );

[pre_im_comp1, err_im_comp1, par.im_comp]       =  GAP_coding(par.I, par.tau, par.win);


%%use GAP predict the compressed image
[pre_im_ori, err_im_ori,~]              = GAP_predict(par.I);  %pre_im_ori: prediction value of original image
[pre_im_comp, err_im_comp,dir_flag]     = GAP_predict(par.im_comp);     % the err_im_comp should be equal to err_im_comp1

dis_low = -100; dis_up = 100;
 draw_error_distribution(err_im_ori, 'ori error', dis_low, dis_up);
 draw_error_distribution(err_im_comp, 'calic error', dis_low, dis_up);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%remove the footprint
[par.im_removed_pf,err_im_removed_pf]   =  remove_fp(err_im_comp1,dir_flag,par.im_comp,pre_im_comp,par.tau,par.win);
%%%%%%%%%%%%%%%%%%%%
dis_low = -100; dis_up = 100;

[pre_im_removed, err_im_removed,~]      = GAP_predict(par.im_removed_pf);
%imshow(uint8(pre_im));
% draw_error_distribution(err_im_ori, 'prediction error', dis_low, dis_up);
 %draw_error_distribution(err_im_comp, sprintf('quantized prediciton error (step = %d)', 2*tolerance+1), dis_low, dis_up);

draw_error_distribution(err_im_removed, sprintf('new prediciton error (step = %d)', 2*tolerance+1), dis_low, dis_up);

psrn_com_re = psnrfun(par.im_comp,par.im_removed_pf,[0 0])
psrn_re_I   = psnrfun(par.I,par.im_removed_pf,[0 0])
psrn_com_I  = psnrfun(par.im_comp,par.I,[0 0])
maxerror    = max(max(abs(par.im_comp-par.im_removed_pf)))

img_save        = strcat(DataBase_Cur,'anti_images\tau',int2str(tolerance),'\Satellite',int2str(imgidx), '.png');
imwrite(uint8(par.im_removed_pf),img_save);

% temp = par.im_removed_pf;
% save im_removed_pf.mat temp
figure;
imshow(uint8(par.im_removed_pf));
end
