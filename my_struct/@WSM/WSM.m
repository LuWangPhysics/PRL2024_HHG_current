classdef WSM
properties
    stheta;
    ctheta;
    sphi;
    cphi;
    dtheta;
    dphi;
    pi_vec={0,0,0};
    Energy_omega;
    %initial condition 
    Rho;
    Rho0;
    Gamma;
    tau_decay;
    q;M;
    dt;dkx;dky;dkz;
    J;
    gamma_int;
    rho_int;
    switch2D=1;
    A0;
    my_d;
    y_window;
    %y window
    N_left;
    N_right;
    ac;av;
    j_inter;
    j_intra;
    eta_num;

end
methods
    
    function obj=init(obj,K,E,C,A_peak)
        obj.dt=E.dt;
        %charge
        obj.q=-C.e;
        %coefficient of the 2nd case two-nodes connected
        obj.M=3e-10;
        obj.dkx=diff(K.k_vec{1});
        obj.dky=diff(K.k_vec{2});
        obj.dkz=K.k_vec{3}(2)-K.k_vec{3}(1);
        %deal with the singularity
        obj.A0=max(abs(A_peak));
        obj.y_window=10*C.omega/C.c;
        obj.my_d=0*(C.hbar*obj.y_window*C.v(1));
        obj.N_left=fix((length(obj.dky)+1)/2-obj.y_window/obj.dky(1));
        obj.N_right=fix((length(obj.dky)+1)/2+obj.y_window/obj.dky(1));
        obj.tau_decay=10e-15;%C.tau;

       obj.eta_num=0;
          
        %for 2D case 
        if K.kz_N==1
            obj.dkz=2*pi;
            obj.switch2D=0;
        end
        obj.J=[0,0,0];
        obj.j_inter=0;
        obj.j_intra=0;
        
        

    end




   function obj=reset_initial(obj,C,chi,kz_iter,K)
                Fermi=@(E) (1+exp((E-C.Ef)./(C.kb*C.T))).^-1;      
                obj.Gamma=0;
                 
                pi_x=C.hbar.*(K.k_vec{1}-chi*C.b(1)).*C.v(1);      
                pi_y=C.hbar.*(K.k_vec{2}-chi*C.b(2)).*C.v(2);      
                pi_z=C.hbar.*(K.k_vec{3}(kz_iter)-chi*C.b(3)).*C.v(3);
                E=sqrt(pi_x.^2+pi_y'.^2+pi_z.^2);
                obj.Energy_omega=0;
                obj.Rho=Fermi(E)-Fermi(-E);
                obj.Rho0=obj.Rho;
                obj.ac=0.*obj.Rho;
                obj.av=0.*obj.Rho-1;
                obj.eta_num=0;

   end


   
    function [obj]=gamma_rho(obj,chi,C)

            %gamma=acav^*
            %each kz overwrite 
            obj.Energy_omega=obj.Energy_omega+obj.dt.*sqrt(obj.pi_vec{1}.^2+obj.pi_vec{2}'.^2+obj.pi_vec{3}.^2)./C.hbar;
            temp=obj.Gamma.*exp(-2i.*obj.Energy_omega);

            f_rho=-chi*2.*obj.dtheta.*real(temp)+2.*obj.dphi.*obj.stheta.*imag(temp)-(1/obj.tau_decay).*(obj.Rho-obj.Rho0);
            f_gamma= chi*0.5*obj.Rho.*exp(2i.*obj.Energy_omega).*(obj.dtheta-chi*1i.*obj.stheta.*obj.dphi)...
                +chi*1i*obj.ctheta.*obj.dphi.*obj.Gamma-(1/obj.tau_decay).*obj.Gamma;
            
            obj.Rho=f_rho.*obj.dt+obj.Rho;
            % obj.Rho(obj.Rho<-1)=-1;
            %  obj.Rho(obj.Rho>1)=1;
            obj.Gamma=f_gamma*obj.dt+obj.Gamma;
            %obj.Gamma(abs(obj.Gamma)>1)=  obj.Gamma(abs(obj.Gamma)>1)./abs( obj.Gamma(abs(obj.Gamma)>1));



    end

    function [obj]=gamma_rho_frozen_band_num(obj,chi,C)

            %gamma=acav^*
            %each kz overwrite 
            obj.Energy_omega=obj.Energy_omega+obj.dt.*sqrt(obj.pi_vec{1}.^2+obj.pi_vec{2}'.^2+obj.pi_vec{3}.^2)./C.hbar;
            temp=obj.Gamma.*exp(-2i.*obj.Energy_omega);

            f_t_temp=0.5*1i.*chi.*obj.ac.*obj.dphi.*obj.ctheta;
            g_t_temp=0.5.*obj.av.*(1i.*obj.dphi.*obj.stheta-chi.*obj.dtheta).*exp(2i.*obj.Energy_omega);

            obj.ac=(f_t_temp+g_t_temp-0.5.*obj.ac/obj.tau_decay).*obj.dt+obj.ac;

            obj.Rho=abs(obj.ac).^2-1;
            % obj.Rho(obj.Rho<-1)=-1;
            %  obj.Rho(obj.Rho>1)=1;
            obj.Gamma=obj.ac.*conj(obj.av);
            %obj.Gamma(abs(obj.Gamma)>1)=  obj.Gamma(abs(obj.Gamma)>1)./abs( obj.Gamma(abs(obj.Gamma)>1));



    end





    function obj=J_calculate(obj,chi,C,K)
            NN=obj.Rho+1;
            temp=obj.Gamma.*exp(-2i.*obj.Energy_omega);
            coe1={obj.stheta.*obj.cphi,obj.stheta.*obj.sphi,obj.ctheta};
            coe2={-chi*obj.cphi.*obj.ctheta+1i.*obj.sphi,-chi*obj.sphi.*obj.ctheta-1i.*obj.cphi,chi.*obj.stheta};   
            %----------------------------------
            %brute force square mesh integral
            %----------------------------------

               for i_iter=1:3
                   J_kxky=NN.*coe1{i_iter}+2.*real(temp.*coe2{i_iter});
                   if(i_iter==2)
                     obj.j_inter=obj.my_sum(2.*real(temp.*coe2{i_iter}))*C.g*obj.q.*C.v(i_iter);
                     obj.j_intra=obj.my_sum(NN.*coe1{i_iter})*C.g*obj.q.*C.v(i_iter);
                   end
                   obj.J(i_iter)=obj.my_sum(J_kxky)*C.g*obj.q.*C.v(i_iter);
   
               end

            %gamma int rho int
              obj.gamma_int=obj.my_sum(abs(obj.Gamma));
               obj.rho_int=obj.my_sum(obj.Rho);

   

    end



    function out=my_sum(obj,matrix_sum)
                 hold_y=zeros(1,length(obj.dky)+1);
                    for ky_iter=1:length(hold_y)
                         hold_y(ky_iter)=0.5*sum(obj.dkx.*matrix_sum(ky_iter,1:end-1)+obj.dkx.*matrix_sum(ky_iter,2:end));
                    end


                   out=0.5*sum(hold_y(1:end-1).*obj.dky+hold_y(2:end).*obj.dky).*obj.dkz/(2*pi)^3;
    end
   
end 

end
