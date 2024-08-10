classdef CONS 
    properties (Constant)
    c=299792458;
    hbar=1.0545718176461565e-34;
    e=1.602176634e-19;
    kb=1.380649e-23;
    eps=8.854187817620389e-12;
    %degenerate coefficient
    g=2;
    T=300;

    end
    properties
    Ef;
    f_fermi
    v;
    n_fermi;
    k_fermi;
    lambda;
    omega;
    k_range;
    tau;
    b;
    wsm_a;
    wsm_b;
    wsm_c;
    
    end
    methods
        function obj=init(obj)
            obj.omega=obj.c*2*pi/obj.lambda;

            %to couple polarisations, b and k of the field need to be
            %aligned ,bx corresponds to Ey Ez
            %for optical 
            % if obj.lambda<2000e-9
            %     %for optical at 1000nm
            %      obj.b=[3.5*obj.omega/obj.v(1),0,0];
            % else
            %     %for thz at 30um=30x1000nm
            %      obj.b=[30*3.5*obj.omega/obj.v(1),0,0];
            % end
             obj.wsm_a=3.4e-10;
             obj.wsm_b=3.4e-10;
             obj.wsm_c=11.6e-10;
             obj.b= [0.06*pi/obj.wsm_a,0,0];     
        end

    end
end