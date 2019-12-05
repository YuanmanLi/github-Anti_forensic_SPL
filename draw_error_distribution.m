%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Created by Li Yuanman 
%% Jan. 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [steps, N] = draw_error_distribution(error_mat, label,low, up)
%% label:  for xlabel
%% low, up:  the range of display, for better vision
if(nargin < 2)
    label = 'no label';
elseif nargin < 3
    low = -256;
    up   = 256;
elseif nargin < 4
    up   = 256;
end

error_mat      =       double(round(error_mat));
error_col       =      error_mat(:);
min_value     =        min(error_col);
max_value     =        max(error_col);

min_value     =         max(low, min_value);
max_value     =         min(up, max_value);

steps               =       min_value : max_value;
N                   =        hist(error_col, steps);
N                   =        N./sum(N);
figure;
%plot(steps, N);
bar(steps, N,'b');
%xlabel(label)
%ylabel('probability')
xlabel(label, 'fontsize', 18)
ylabel('probability', 'fontsize',18)
%xlabel('\fontsize{18}Line Number')
%ylabel('\fontsize{18}Mean Value')
set(gca,'FontSize',15)

hold on

%%%%%%%%%%%%%%%%%%%%
%%use MAP to get the laplacian parameter b. 
error_mat_abs      =  abs(error_mat);
[width, height]     =  size(error_mat);
b                           =  sum(sum(error_mat_abs))/(width*height);

%%%%%%%%%%%%%%%%%%%%
%%draw laplacian. 
mu            = 0;
f1              = @(x)0.5/b * exp(-abs(x-mu)*1.0/b);
max_val   = max(max(error_mat));
min_val    = min(min(error_mat));
x               =  min_val : 0.1: max_val;
%plot(x,f1(x),'r-'); 
xlim([low up]);
%legend('prediction error',sprintf('Laplacian: b = %.2f',b));
end