classdef Field
properties
    tau;
    E_vec;
    A_vec;
    dt;
    test;
    t_shift;
    t;
    A_peak_arr;
    E0;
    Ay;
end
methods
    function obj=init(obj,C,E_0,phi_cep,t)
        obj.tau=C.tau;
        obj.t=t;
       % obj.t=linspace(-55e-15,55e-15,1000);
       obj.E0=E_0;
        obj.dt=obj.t(2)-obj.t(1);
   
        %circularly polarized light phi=pi/2
        phi=0*pi/2;
        %E_vec{Ex,Ey,Ez}
        if obj.t(end)<3*obj.tau
           obj.t_shift=obj.t;
        else
            %shift the incident pulse position so that it is close to t(0)
            obj.t_shift=obj.t+(obj.t(end)-3*obj.tau);
        end
        env=E_0.*exp(-obj.t_shift.^2./obj.tau.^2);
 
        %env=
        obj.E_vec{1}=0.*real(env.*exp(1i*C.omega.*obj.t_shift));
        obj.E_vec{2}=1.*real(env.*exp(1i*C.omega.*obj.t_shift+1i*phi_cep));
        obj.E_vec{3}=0.*real(env.*exp(1i*C.omega.*obj.t_shift+1i*phi));
        obj.Ay=-cumsum(obj.E_vec{2}).*obj.dt;
       
         obj.A_vec=zeros(1,3);
         obj.test=zeros(size(obj.t));
          
         obj.A_peak_arr=zeros(1,3);

    end

    function obj=A(obj,t_iter)
       t_decay=1*obj.tau;
        for n_iter=1:3
           obj.A_vec(n_iter)=-obj.dt*sum(obj.E_vec{n_iter}(1:t_iter));
            %inclde absorption
          %   obj.A_vec(n_iter)=-exp(-obj.t_shift(t_iter)/t_decay).*obj.dt*sum(exp(obj.t_shift(1:t_iter)/t_decay).*obj.E_vec{n_iter}(1:t_iter));
           % obj.A_vec(n_iter)=exp(-obj.dt/t_decay)*obj.A_vec(n_iter)-obj.dt*obj.E_vec{n_iter}(t_iter);
        end
       %obj.test(t_iter)=obj.A_vec(2);
    end

    function obj=A_peak(obj)
     
        for n_iter=1:3
           obj.A_peak_arr(n_iter)=max(abs(-obj.dt*cumsum(obj.E_vec{n_iter})));

        end

    end

end



end
