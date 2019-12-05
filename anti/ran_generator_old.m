function y= ran_generator(lambda,maxval,numsamp,flag)


% default 
maxiter= 1e6;

%
y= zeros(numsamp,1);%numsamp=1024

%c= (lambda*maxval)/(1-exp(-lambda*maxval));
%beta=(1/lambda)*(1-exp(-lambda*maxval));

iter= 1;
sampnum= 1;
if(flag == 1)
    while (sampnum <= numsamp) | (maxiter <= iter) 
        % generate samples from proposed (uniform) distribution
        x= maxval*rand(1,1);   %maxval就是当前的量化step
                                                %rand(1,1)生成一个0-1的随机数， 均匀分布
        u= rand(1,1);

        if u <= exp(-lambda*x) % %其实案例所已经先求累计分布的反函数。  其实这里用的 accept-reject method.其实前面的系数并不重要。。  因为相当于取到的概率都降低了系数倍。 都降低同样的倍数，其实不影响
            y(sampnum)= x;
            sampnum= sampnum+1;
        end

        iter= iter+1;
    end
elseif(flag == 0)
    while (sampnum <= numsamp) | (maxiter <= iter) 
        % generate samples from proposed (uniform) distribution
        x= maxval*(rand(1,1)-.5);  %这个是处理量化等于0的情况  -Q/2 <=x <= Q/2

        u= rand(1,1);

        if u <= exp(-lambda*abs(x)) %没有看明白 这里是怎么采样的。 当x的绝对值越大，后边的部分越小，那么不等式就更不容易成立，所以一定程度上可能是可行的
            y(sampnum)= x;          %其实案例所已经先求累计分布的反函数。  其实这里用的 accept-reject method.其实前面的系数并不重要。。  因为相当于取到的概率都降低了系数倍。 都降低同样的倍数，其实不影响
            sampnum= sampnum+1;
        end

        iter= iter+1;
    end  
else
    error('flage can only be 0 or 1');
end