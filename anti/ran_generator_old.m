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
        x= maxval*rand(1,1);   %maxval���ǵ�ǰ������step
                                                %rand(1,1)����һ��0-1��������� ���ȷֲ�
        u= rand(1,1);

        if u <= exp(-lambda*x) % %��ʵ�������Ѿ������ۼƷֲ��ķ�������  ��ʵ�����õ� accept-reject method.��ʵǰ���ϵ��������Ҫ����  ��Ϊ�൱��ȡ���ĸ��ʶ�������ϵ������ ������ͬ���ı�������ʵ��Ӱ��
            y(sampnum)= x;
            sampnum= sampnum+1;
        end

        iter= iter+1;
    end
elseif(flag == 0)
    while (sampnum <= numsamp) | (maxiter <= iter) 
        % generate samples from proposed (uniform) distribution
        x= maxval*(rand(1,1)-.5);  %����Ǵ�����������0�����  -Q/2 <=x <= Q/2

        u= rand(1,1);

        if u <= exp(-lambda*abs(x)) %û�п����� ��������ô�����ġ� ��x�ľ���ֵԽ�󣬺�ߵĲ���ԽС����ô����ʽ�͸������׳���������һ���̶��Ͽ����ǿ��е�
            y(sampnum)= x;          %��ʵ�������Ѿ������ۼƷֲ��ķ�������  ��ʵ�����õ� accept-reject method.��ʵǰ���ϵ��������Ҫ����  ��Ϊ�൱��ȡ���ĸ��ʶ�������ϵ������ ������ͬ���ı�������ʵ��Ӱ��
            sampnum= sampnum+1;
        end

        iter= iter+1;
    end  
else
    error('flage can only be 0 or 1');
end