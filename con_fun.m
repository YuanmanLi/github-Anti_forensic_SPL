function [c , ceq]  =  con_fun(x, lambda, pre, prepre, boundary_next_p,quantized_step, Nk, N, alpha1, alpha2)

f          =  @(val,para_lambda)para_lambda/2*exp(-para_lambda.*abs(val));
c_bin = pre:-sign(pre): pre-sign(pre)*(quantized_step-1);    %this should be differen for different sign
Ui  =  @(para_lambda)(Nk/N)/sum(f(c_bin, para_lambda));
pre = f(pre, x)*Ui(x); 
prepre = f(prepre, x)*Ui(x);
c1         =   alpha1*(prepre-pre) - (pre-boundary_next_p);
c2         =   alpha2*(prepre-pre) - (pre-boundary_next_p);
%alpha*prepre + boundary_next_p -  (1+alpha)*pre;
c = [c1, -c2];
%c = -c;
ceq = [];