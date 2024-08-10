    function obj=theta_phi(obj,chi,E,C,K,t_iter,kz_iter)

      %calculate theta and phi for a given kz;
        for n_iter=1:2
            obj.pi_vec{n_iter}=(C.hbar.*(K.k_vec{n_iter}-chi*C.b(n_iter))-obj.q.*E.A_vec(n_iter)).*C.v(n_iter);      
        end
        obj.pi_vec{3}=(C.hbar.*(K.k_vec{3}(kz_iter)-chi*C.b(3))-obj.q.*E.A_vec(3)).*C.v(3);

        %x horisontal, y vertical
        
        Pi_xy=sqrt(obj.pi_vec{1}.^2+obj.pi_vec{2}'.^2+obj.my_d.^2);

       %remove the sigular part
        % N_b=find(min(abs(obj.pi_vec{1}))==abs(obj.pi_vec{1}));
        % 
        % x_select=N_b+(-fix(obj.y_window/obj.dky(N_b)):(fix(obj.y_window/obj.dky(N_b))));
        % y_select=obj.N_left:obj.N_right;
        % if(min(min(Pi_xy))<obj.my_d)
        % 
        %     Pi_xy(y_select,x_select)=obj.my_d;
        % end
        Pi= sqrt(Pi_xy.^2+obj.pi_vec{3}.^2);
        

        obj.stheta=Pi_xy./Pi;
        obj.ctheta=obj.pi_vec{3}./Pi;
 
        obj.sphi=obj.pi_vec{2}'./Pi_xy;
        obj.cphi=obj.pi_vec{1}./Pi_xy;
        pix_dot=obj.q.*C.v(1).*E.E_vec{1}(t_iter);
        piy_dot=obj.q.*C.v(2).*E.E_vec{2}(t_iter);
        piz_dot=obj.q.*C.v(3).*E.E_vec{3}(t_iter);
        obj.dtheta=-obj.switch2D.*(obj.stheta.*piz_dot-obj.ctheta.*(obj.cphi.*pix_dot+obj.sphi.*piy_dot))./Pi;
        obj.dphi=(piy_dot.*obj.cphi-pix_dot.*obj.sphi)./Pi_xy;


    end

     % 
    % function obj=theta_phi2(obj,chi,E,C,K,t_iter,kz_iter)
    % 
    %   %assume node separation is at x
    % 
    %   %calculate theta and phi for a given kz;
    % 
    % 
    %     k_mod2=(K.k_vec{3}(kz_iter)-obj.q.*E.A_vec(3)/C.hbar).^2+(K.k_vec{2}'-obj.q.*E.A_vec(2)/C.hbar).^2+(K.k_vec{1}-obj.q.*E.A_vec(1)./C.hbar).^2;
    %     obj.pi_vec{1}=C.hbar.*obj.M.*(k_mod2-C.b(1)^2).*C.v(1);
    %     for n_iter=2:3
    %          obj.pi_vec{n_iter}=(C.hbar.*K.k_vec{n_iter}-obj.q.*E.A_vec(n_iter)).*C.v(n_iter);      
    %     end
    % 
    %     %x horisontal, y vertical
    %     Pi= sqrt(obj.pi_vec{1}.^2+obj.pi_vec{2}'.^2+obj.pi_vec{3}.^2);
    %     Pi_xy=sqrt(obj.pi_vec{1}.^2+obj.pi_vec{2}'.^2);
    %     obj.stheta=Pi_xy./Pi;
    %     obj.ctheta=obj.pi_vec{3}./Pi;
    %     obj.sphi=obj.pi_vec{2}'./Pi_xy;
    %     obj.cphi=obj.pi_vec{1}./Pi_xy;
    %     pix_dot=2.*C.v(1)*obj.M*((K.k_vec{1}-obj.q*E.A_vec(1)/C.hbar).*obj.q*E.E_vec{1}(t_iter)...
    %        +(K.k_vec{2}'-obj.q*E.A_vec(2)/C.hbar).*obj.q*E.E_vec{2}(t_iter) ...
    %        +(K.k_vec{3}(kz_iter)-obj.q*E.A_vec(3)/C.hbar).*obj.q*E.E_vec{3}(t_iter));
    %     piy_dot=obj.q.*C.v(2).*E.E_vec{2}(t_iter);
    %     piz_dot=obj.q.*C.v(3).*E.E_vec{3}(t_iter);
    %     obj.dtheta=-obj.switch2D.*(obj.stheta.*piz_dot-obj.ctheta.*(obj.cphi.*pix_dot+obj.sphi.*piy_dot))./Pi;
    %     obj.dphi=(piy_dot.*obj.cphi-pix_dot.*obj.sphi)./Pi_xy;
    % 
    % 
    % end
   
   % function obj=reset_initial2(obj,C,chi,kz_iter,K)
   %              Fermi=@(E) (1+exp((E-C.Ef)./(C.kb*C.T))).^-1;      
   %              obj.Gamma=0;
   %              k_mod2=(K.k_vec{3}(kz_iter)).^2+(K.k_vec{2}').^2+(K.k_vec{1}).^2;
   %              pi_x=(C.hbar.*obj.M.*(k_mod2-C.b(1)^2)).*C.v(1);
   %              pi_y=(C.hbar.*K.k_vec{2}).*C.v(2);      
   %              pi_z=(C.hbar.*K.k_vec{3}).*C.v(3);  
   % 
   % 
   %              E=sqrt(pi_x.^2+pi_y'.^2+pi_z.^2);
   %              obj.Energy_omega=0;
   %              obj.Rho=Fermi(E)-Fermi(-E);
   %              obj.Rho0=obj.Rho;
   % 
   % end
