classdef mesh_k
    properties
    k_vec;
    kz_N;
    k_flag;
    end
    
    methods  
        function obj=init(obj,C)

           
            obj.k_flag=0;

            % k_win=3*C.b(1);
%include the case when b=0
            k_win=0.9*C.omega/C.v(1);

            %TaAs
            a=C.wsm_a;
            b=C.wsm_b;
            c=C.wsm_c;
            r=[a,b,c];
            k_end=pi./r;
  
            k_margin=1e7;
           
            % N=[400,400,5];
            % range_max=1;
            % obj=obj.k_linear_direct(k_end,N,k_win,range_max);

         N=6*[10,11,10];%even number to avoid zero
         obj=obj.k_around_cone(k_end,N,C,k_win,k_margin);
         obj.kz_N=length(obj.k_vec{3});


        end

        function obj=k_linear_direct(obj,k_end,N,k_win,range_max)

           for iter=1:3

                obj.k_vec{iter}=linspace(-k_end(iter),k_end(iter),N(iter));
            end
           obj.k_vec{3}=obj.k_vec{3}.*range_max;
        end

        

        function obj=k_around_cone(obj,k_end,N,C,k_win,k_margin)
             %1/3 around the nodes
    
             obj.k_flag=1;
             %obj.k_vec{1}=[linspace(-k_win/2,C.b(1)-k_win,fix(N(1)*0.1)),...
                 % linspace(C.b(1)-k_win+k_margin,C.b(1)+k_win+k_margin,fix(0.8*N(1))),...
                 % linspace(C.b(1)+k_win+2*k_margin,k_end(1),0.1*N(1))];
                 if C.b(1)~=0
                obj.k_vec{1}=[linspace(C.b(1)-k_win,C.b(1)+k_win,fix(0.9*N(1))),...
                 linspace(C.b(1)+k_win+k_margin,k_end(1),0.1*N(1))];
                 else
                obj.k_vec{1}=linspace(-k_win,k_win,N(1));
                 end
             obj.k_vec{2}=linspace(-k_win,k_win,N(2));
             obj.k_vec{3}=linspace(-k_win,k_win,N(3));
          
        end

        function [obj,M]=k_arange(obj,chi,M)
               if obj.k_flag==1
                  if (chi==-1&&abs(obj.k_vec{1}(1))~=obj.k_vec{1}(end))
                    obj.k_vec{1}=-flip(obj.k_vec{1});
                  end
                  M.dkx=diff(obj.k_vec{1});
               end
        end


    
    end 
end