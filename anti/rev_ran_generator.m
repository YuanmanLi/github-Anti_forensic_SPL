function d= rev_ran_generator(lambda,step,numsamp,flag)
% default 
maxiter= 1e7;

%
d= zeros(numsamp,1);%numsamp=1024

%c= (lambda*maxval)/(1-exp(-lambda*maxval));
%beta=(1/lambda)*(1-exp(-lambda*maxval));

y = exp(-lambda*step/2);
min_v = -step/2;
max_v = step/2;

iter= 1;
sampnum= 1;
if(flag == 0)
    while (sampnum <= numsamp) | (maxiter <= iter) 
        x = rand(1,1)*(max_v-min_v) + min_v;
        u = rand(1,1);
        if u < lambda/(2*(1-y))*exp(-lambda*abs(x));
            d(sampnum)= x;
            sampnum= sampnum+1;
        end
        iter= iter+1;
    end
elseif(flag == 1)
    while (sampnum <= numsamp) | (maxiter <= iter) 
        % generate samples from proposed (uniform) distribution
        x = rand(1,1)*(max_v-min_v) + min_v;
        u = rand(1,1);
        if u < lambda/(y^(-1)-y)*exp(-lambda*x);
            d(sampnum)= x;
            sampnum= sampnum+1;
        end
        iter= iter+1;
    end  
else
    error('flage can only be 0 or 1');
end