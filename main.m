clear all;
addpath('my_struct');
addpath('my_function');
addpath('my_plot_scheme');
addpath('my_plot_scheme/print_file');


C=CONS;
C.v=[1,1,1]*1e6;
C.Ef=100e-3*C.e;
C.k_fermi=C.Ef/C.hbar/min(C.v);
%------------------------
%optical pulse parameter
%------------------------
C.lambda=1000e-9;
C.tau=10e-15;
E0=2e8;
%------------------------
%thz pulse parameter
%------------------------
% C.lambda=10*1000e-9;
% C.tau=40e-15;
% E0=2e7;
C=C.init();
%---------------------
%initiate pulse
%---------------------

t=linspace(-5.5*C.tau,5.5*C.tau,2400);

phi_cep=0;
E=Field;
E=E.init(C,E0,phi_cep,t);
E=E.A_peak();

%---------------------
%init k mesh
%---------------------
K=mesh_k;
K=K.init(C);
%init WSM response include theta, phi, Rho, Gamma
M=WSM;
M=M.init(K,E,C,E.A_peak_arr);
J={zeros(size(E.t)),zeros(size(E.t)),zeros(size(E.t))};


tic
for chiral=0:1
    chi=(-1)^chiral;
    [K,M]=K.k_arange(chi,M);
     
    for kz_iter=1:K.kz_N

        for t_iter=1:length(E.t)
                if t_iter==1
                  M=M.reset_initial(C,chi,kz_iter,K);
                end
                E=E.A(t_iter);
                M=M.theta_phi(chi,E,C,K,t_iter,kz_iter);
                %----------------
                %analytical methods
                %----------------
                % M=M.frozen_band_analytic_withSymmetry(chi,C,t_iter,E);
                % M=M.frozen_band_eta_num(chi,C,t_iter,E);
                % M=M.gamma_rho_analytic_more_bessel(chi,C,t_iter,E);
                % M=M.gamma_rho_frozen_band_num(chi,C);
                %----------------------
                %numerical methodss
                %----------------------
                M=M.gamma_rho(chi,C);
                M=J_calculate(M,chi,C,K);
                J{1}(t_iter)=J{1}(t_iter)+M.J(1);
                J{2}(t_iter)=J{2}(t_iter)+M.J(2);
                J{3}(t_iter)=J{3}(t_iter)+M.J(3);


        end
    end
end
toc


plot_et(E,J,C)
[Jf,f_arr]=spec(E.t,J,C,1);

