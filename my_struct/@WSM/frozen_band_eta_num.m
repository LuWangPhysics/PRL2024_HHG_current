
    function [obj]=frozen_band_eta_num(obj,chi,C,t_iter,E)
          t0=E.t_shift(t_iter);

          Pi_xy=sqrt(obj.pi_vec{1}.^2+obj.pi_vec{2}'.^2+obj.my_d.^2)./C.v(1);
          Pi_scaler= sqrt(Pi_xy.^2+obj.pi_vec{3}.^2./(C.v(3)^2));
          N=-chi*obj.q*E.E0.*exp(-t0^2./C.tau.^2).*obj.cphi.*obj.ctheta./(2*C.omega.*Pi_xy);
          %exp_A=exp(1i.*N.*sin(C.omega.*t0));
          coe_a=-chi*obj.q.*obj.cphi.*obj.ctheta./(2*C.omega.*Pi_xy);
          exp_A=exp(1i.*coe_a.*E.Ay(t_iter));

          int_mat_ana=(E.E0*C.tau/sqrt(pi)/C.omega).*(erf(t0./C.tau)+1);
          int_mat_num=sum(abs(E.Ay(1:t_iter))).*E.dt;
          e_Omega_full_ana=exp(1i*2*C.v(1)*int_mat_num.*abs(C.e)./(C.hbar));


          coe_angle=(1i*obj.cphi-chi*obj.ctheta.*obj.sphi);
       

          obj.eta_num=obj.eta_num+E.E_vec{2}(t_iter).*exp_A.*e_Omega_full_ana.*E.dt;
          a_c=(-obj.q.*coe_angle./(2.*Pi_scaler)).*conj(exp_A).*obj.eta_num;

          %decay are put in by hand

          coe_decay=exp(-0.5.*(E.t_shift(t_iter)-E.t_shift(1))./obj.tau_decay);
          a_c=a_c.*coe_decay;
          temp=a_c.*conj(e_Omega_full_ana);

            coe1={obj.stheta.*obj.cphi,obj.stheta.*obj.sphi,obj.ctheta};
            coe2={-chi*obj.cphi.*obj.ctheta+1i.*obj.sphi,-chi*obj.sphi.*obj.ctheta-1i.*obj.cphi,chi.*obj.stheta};   
            %----------------------------------
            %brute force square mesh integral
            %----------------------------------
  
               for i_iter=1:3
                   J_kxky=abs(a_c).^2.*coe1{i_iter}+2.*real(temp.*coe2{i_iter});

                   obj.J(i_iter)=obj.my_sum(J_kxky)*C.g*obj.q.*C.v(i_iter);
   
               end


  

    end