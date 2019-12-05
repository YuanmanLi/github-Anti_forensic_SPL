%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Created by Li Yuanman 
%% Jan. 16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [solutions] = solve_GAP_theta(im_comp,flag,quan_err,tau,i,j,pro_zero, pro_szero)  
    solutions = -tau : tau;
   n         = length(flag);   %n should be 7 for calic
   sol_n     = cell(n,1); 
   dh               =  sym('abs(w-ww) + abs(n-nw) + abs(n-ne)');
   dv               =  sym('abs(w-nw)  + abs(n-nn)  + abs(ne-nne)');
   for k = 1: n
       sym r;
       f = dv-dh;
        switch(k)
            case 1;
                i1              =  i; 
                j1               =  j+1;
                nn              =  im_comp(i1-2, j1);
                nne             =  im_comp(i1-2, j1+1);
                nw              =  im_comp(i1-1,j1-1);
                n               =  im_comp(i1-1, j1);
                ne              =  im_comp(i1-1, j1+1);
                ww              =  im_comp(i1, j1-2);
                w               =  im_comp(i1, j1-1);
                f = subs(f,'w','w+r');
            case 2;
                i1              =  i; 
                j1               =  j+2;
                nn              =  im_comp(i1-2, j1);
                nne             =  im_comp(i1-2, j1+1);
                nw              =  im_comp(i1-1,j1-1);
                n               =  im_comp(i1-1, j1);
                ne              =  im_comp(i1-1, j1+1);
                ww              =  im_comp(i1, j1-2);
                w               =  im_comp(i1, j1-1);
                f = subs(f,'ww','ww + r');
            case 3;
                i1              =  i+1; 
                j1               =  j-1;
                nn              =  im_comp(i1-2, j1);
                nne             =  im_comp(i1-2, j1+1);
                nw              =  im_comp(i1-1,j1-1);
                n               =  im_comp(i1-1, j1);
                ne              =  im_comp(i1-1, j1+1);
                ww              =  im_comp(i1, j1-2);
                w               =  im_comp(i1, j1-1);
                f = subs(f,'ne','ne + r');
            case 4;
                i1              =  i+1; 
                j1               =  j;
                nn              =  im_comp(i1-2, j1);
                nne             =  im_comp(i1-2, j1+1);
                nw              =  im_comp(i1-1,j1-1);
                n               =  im_comp(i1-1, j1);
                ne              =  im_comp(i1-1, j1+1);
                ww              =  im_comp(i1, j1-2);
                w               =  im_comp(i1, j1-1);
                f = subs(f,'n','n + r');
            case 5;
                i1              =  i+1; 
                j1               =  j+1;
                nn              =  im_comp(i1-2, j1);
                nne             =  im_comp(i1-2, j1+1);
                nw              =  im_comp(i1-1,j1-1);
                n               =  im_comp(i1-1, j1);
                ne              =  im_comp(i1-1, j1+1);
                ww              =  im_comp(i1, j1-2);
                w               =  im_comp(i1, j1-1);
                f = subs(f,'nw','nw + r');
            case 6;
                i1              =  i+2; 
                j1               =  j-1;
                nn              =  im_comp(i1-2, j1);
                nne             =  im_comp(i1-2, j1+1);
                nw              =  im_comp(i1-1,j1-1);
                n               =  im_comp(i1-1, j1);
                ne              =  im_comp(i1-1, j1+1);
                ww              =  im_comp(i1, j1-2);
                w               =  im_comp(i1, j1-1);
                f = subs(f,'nne','nne + r');
            case 7;
                i1              =  i+2; 
                j1               =  j;
                nn              =  im_comp(i1-2, j1);
                nne             =  im_comp(i1-2, j1+1);
                nw              =  im_comp(i1-1,j1-1);
                n               =  im_comp(i1-1, j1);
                ne              =  im_comp(i1-1, j1+1);
                ww              =  im_comp(i1, j1-2);
                w               =  im_comp(i1, j1-1);
                f = subs(f,'nn','nn + r');
        end
        curr_dir    = flag(k);
        curr_sol    = get_sol(f,nn,nne,nw,n,ne,ww,w,tau,curr_dir); 
        sol_n{k}    = curr_sol;
        solutions   = intersect(solutions, sol_n{k});
        if(isempty(solutions))
            
%             fprintf('no solution for i=%d, j=%d. the %dth neighboor no solution!',i,j,k);
%             error('error!');
%             break;
          %% only for test if no solution, we just random choose a value
%           if quan_err > 0
%              s = -1;
%           else 
%              s = 1;
%           end
            
          temp = rand(1,1);
          indx = 1;
          if quan_err ==0
              while (indx <= tau+1)
                  if (temp < 2 * sum(pro_zero(1:indx)))
                     break; 
                  end
                  indx = indx + 1;
              end
              solutions = sign(rand(1,1)-0.5)*(indx - tau -1);
          else
              s = -1 + 2*(quan_err < 0);
               while (indx <= 2*tau+1)
                  if (temp < sum(pro_szero(1:indx)) | indx == 2*tau+1)
                     break; 
                  end
                  indx = indx + 1;
               end   
              solutions =s*(indx - tau -1);
          end
        
          solutions = round(solutions);
           %fprintf('no solution i=%d, j=%d. the %dth neighboor no solution! random sol=%d\n',i,j,k,solutions);
              
          
        end
   end
end

function [sol] = get_sol(f,nn,nne,nw,n,ne,ww,w,tau,dir_flag)
  r = -tau : tau;
  sol = [];
  switch(dir_flag)
     case 1;
         sol = r(eval(f)>80);
%        for r = -tau : tau 
%            if(eval(f)>80)
%                sol = [sol r];
%            end
%        end
       case 2;
          sol = r(eval(f)<-80);
%        for r = -tau : tau 
%            if(eval(f)<-80)
%                sol = [sol r];
%            end
%        end
       case 3;
          sol = r(eval(f)>32 & eval(f)<=80);
%        for r = -tau : tau 
%            if(eval(f)>32 & eval(f)<=80)
%                sol = [sol r];
%            end
%        end
       case 4;
           sol = r(eval(f)<=32 & eval(f)>8);
%        for r = -tau : tau 
%            if(eval(f)<=32 & eval(f)>8)
%                sol = [sol r];
%            end
%        end
        case 5;
            sol = r(eval(f)<-32 & eval(f)>=-80);
%        for r = -tau : tau 
%            if(eval(f)<-32 & eval(f)>-80)
%                sol = [sol r];
%            end
%        end
       case 6;
            sol = r(eval(f)>=-32 & eval(f)<-8);
%        for r = -tau : tau 
%            if(eval(f)>=-32 & eval(f)<-8)
%                sol = [sol r];
%            end
%        end
        case 7;
            sol = r(eval(f)>=-8 & eval(f)<=8);
%        for r = -tau : tau 
%            if(eval(f)>=-8 & eval(f)<=8)
%                sol = [sol r];
%            end
%        end
  end
   
end