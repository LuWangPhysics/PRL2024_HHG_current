    function [obj]=gamma_rho_analytic_more_bessel(obj,chi,C,t_iter,E)
          
       %expOmega use the taylor expansion
        t0=E.t_shift(t_iter);

          Pi_xy=sqrt(obj.pi_vec{1}.^2+obj.pi_vec{2}'.^2+obj.my_d.^2)./C.v(1);
          N=-chi*obj.q*E.E0.*exp(-t0^2./C.tau.^2).*obj.cphi.*obj.ctheta./(2*C.omega.*Pi_xy);

        
         Omega_0=4*C.v(1)*abs(obj.q)*E.E0/(C.hbar*pi*C.omega);
         coe_decay=exp(-t0.^2./C.tau.^2).*exp(-0.5.*(E.t_shift(t_iter)-E.t_shift(1))./obj.tau_decay);
     
          coe1={obj.stheta.*obj.cphi,obj.stheta.*obj.sphi,obj.ctheta};
          coe2={-chi*obj.cphi.*obj.ctheta+1i.*obj.sphi,-chi*obj.sphi.*obj.ctheta-1i.*obj.cphi,chi.*obj.stheta};
          nb=6;
         %--------------------------------
         %expand only 1 JA series
         %--------------------------------
          exp_A=exp(1i.*N.*sin(C.omega.*t0));

          int_mat_taylor=(E.E0*C.tau/sqrt(pi)/C.omega).*(1+2*t0./C.tau./sqrt(pi)-0*2.*t0.^3./3./(sqrt(pi)*C.tau^3));
          e_Omega_full_ana=exp(1i*2*C.v(1)*int_mat_taylor.*abs(C.e)./(C.hbar));

            bessel_expnw=bessel_sum(nb,N,t0,C.omega,Omega_0,'coe_all');   
            coe=chi.*obj.stheta./obj.ctheta+1i*obj.stheta.*obj.sphi./obj.cphi;
            a_c=coe.*exp(2i*C.v(1)*abs(obj.q)*E.E0*C.tau/(C.hbar*sqrt(pi)*C.omega))...
            .*bessel_expnw.*conj(exp_A);

           % decay are put in by hand
            a_c=a_c.*coe_decay;
            temp=a_c.*conj(e_Omega_full_ana);


            for i_iter=1:3
               J_kxky=abs(a_c).^2.*coe1{i_iter}+2.*real(temp.*coe2{i_iter});
                obj.J(i_iter)=obj.my_sum(J_kxky)*C.g*obj.q.*C.v(i_iter);

            end
         %--------------------------------
         %expand two JA series in sum
         %--------------------------------
          % c1=obj.stheta.^2.*(obj.cphi.^2+obj.ctheta.^2.*obj.sphi.^2)./obj.cphi.^2./obj.ctheta.^2;
          % J_arr=0;
          % J_nocoe=0;
          % for i_iter=-nb:1:nb
          % J_arr=J_arr+(i_iter/(i_iter+Omega_0/C.omega)).*besselj(i_iter,N)*exp(1i*C.omega*i_iter*t0);
          % J_nocoe=J_nocoe+besselj(i_iter,N)*exp(1i*C.omega*i_iter*t0);
          % end
          % ac_abs2=coe_decay.^2.*c1.*J_arr.*conj(J_arr);
          % c2=chi.*obj.stheta./obj.ctheta+1i*obj.stheta.*obj.sphi./obj.cphi;
          % ac_exp=coe_decay.*c2.*J_arr.*conj(J_nocoe);
          % 
          %      for i_iter=1:3
          %          J_kxky=ac_abs2.*coe1{i_iter}+2.*real(ac_exp.*coe2{i_iter});
          %          obj.J(i_iter)=obj.my_sum(J_kxky)*C.g*obj.q.*C.v(i_iter);
          % 
          %      end


  

    end

