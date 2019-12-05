%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Created by Li Yuanman 
%% Jan. 10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  [quan_val ] = quan2 ( val , tau)

if(tau ~= fix(tau))
    error('tau should be integer!');
end

mask         =    double(2*tau + 1);
quan_val  =   mask * fix(double((val + sign(val)*tau))/ mask);   % fix, only keep the integer part
