    function [obj]=frozen_band_analytic_withSymmetry(obj,chi,C,t_iter,E)
          Pi_xy=sqrt(obj.pi_vec{1}.^2+obj.pi_vec{2}'.^2+obj.my_d.^2)./C.v(1);
  
          t0=E.t_shift(t_iter);
          N=-E.E0.*exp(-t0^2/C.tau^2)*chi*obj.q.*obj.cphi.*obj.ctheta./(2*C.omega.*Pi_xy);
          Omega_0=4*C.v(1)*abs(obj.q)*E.E0/(C.hbar*pi*C.omega);
          U=Omega_0/C.omega;
          %always even number
          nb=6;
          coe_decay=exp(-t0.^2./C.tau.^2).*exp(-0.5.*(E.t_shift(t_iter)-E.t_shift(1))./obj.tau_decay);
          %---------------------
          %calculate Jx
          %---------------------
          obj.J(1)=0;
          %---------------------
          %calculate Jy
          %---------------------
       
          J_co=0;
          J_ce=0;
          J_o=0;
          J_e=0;
          for i_iter=-nb:2:nb
          J_ce=J_ce+(U/(i_iter+U)).*besselj(i_iter,N)*exp(1i*C.omega*i_iter*t0);
          J_co=J_co+(U/(i_iter+1+U)).*besselj(i_iter+1,N)*exp(1i*C.omega*(i_iter+1)*t0);
          J_o=J_o+besselj(i_iter+1,N)*exp(1i*C.omega*(i_iter+1)*t0);
          J_e=J_e+besselj(i_iter,N)*exp(1i*C.omega*(i_iter)*t0);
          end
          J_o=J_o+besselj(-(nb+1),N)*exp(-1i*C.omega*(nb+1)*t0);
          J_co=J_co+(U/(-nb-1+U)).*besselj(-(nb+1),N)*exp(-1i*C.omega*(nb+1)*t0);
      
          jyc2=(obj.cphi./obj.ctheta+obj.sphi.^2.*obj.ctheta./obj.cphi);
          J_kxky=-coe_decay.*2*chi*obj.stheta.*jyc2.*imag(J_ce.*conj(J_o)+J_co.*conj(J_e));
          obj.J(2)=obj.my_sum(J_kxky)*C.g*obj.q.*C.v(2);


          %---------------------
          %calculate Jz
          %---------------------
  
          jz1=coe_decay.^2.*obj.stheta.^2.*(obj.cphi./obj.ctheta+obj.sphi.^2.*obj.ctheta./obj.cphi);
          jz2=-coe_decay.*2.*(obj.stheta.^2./obj.ctheta);
          J_kxky=2.*jz1.*real(J_ce.*conj(J_co))+jz2.*real((J_ce.*conj(J_o)+J_co.*conj(J_e)));

          obj.J(3)=obj.my_sum(J_kxky)*C.g*obj.q.*C.v(3);

            
 
      

    end