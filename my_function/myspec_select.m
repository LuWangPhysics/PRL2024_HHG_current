function [Ey,Ez]=myspec_select(Jf,f,t,omega,harmonic_order,title_str,tau,fig_flag)

N_half=fix(length(t)/2)+1;
f_half=f(N_half:end);
if tau<20e-15
t_scale=1;
fwin=1.1;
else
 t_scale=12;
 fwin=1;
end
%--------------------------------------
%get the spectrum of the second harmonic
%--------------------------------------
df=f_half(2)-f_half(1);
if harmonic_order==1
N_peak=fix(harmonic_order*(omega/2/pi)/df);
N_win_l=1;
N_win_r=N_peak+fix(fwin*(omega/2/pi)/df);
else
N_peak=fix(harmonic_order*(omega/2/pi)/df);
N_win_l=N_peak-fix(fwin*(omega/2/pi)/df);
N_win_r=N_peak+fix(fwin*(omega/2/pi)/df);
end

for i_iter=2:3
    sec=Jf{i_iter}(N_win_l:N_win_r);
    if i_iter==2
        E_secondfy=zeros(size(Jf{i_iter}));
        if(max(abs(sec))>fwin*(abs(sec(1))))
        E_secondfy(N_win_l:N_win_r)=sec;
        Ef=[flip(conj(E_secondfy)),0,E_secondfy];
        Et=fftshift(ifft(ifftshift(Ef))).*(length(t)*df);
        Ey=Et(1:end-1);
        else
        Ey=zeros(1,length(Jf{i_iter})*2);
        E_secondfy=sec(1)*ones(size(Jf{i_iter}));
        end
        %add one zero in the middle is important because it makes sure the results
        %are real after FT

    elseif i_iter==3
    E_secondfz=zeros(size(Jf{i_iter}));
        if(max(abs(sec))>fwin*(abs(sec(1))))
        E_secondfz(N_win_l:N_win_r)=sec;
        Ef=[flip(conj(E_secondfz)),0,E_secondfz];
        Et=fftshift(ifft(ifftshift(Ef))).*(length(t)*df);
        Ez=Et(1:end-1);
        else
        Ez=zeros(1,length(Jf{i_iter})*2); 
        E_secondfz=sec(1).*ones(size(Jf{i_iter}));
        end
        %add one zero in the middle is important because it makes sure the results
        %are real after FT

    end


end
t_arr=linspace(t(1),t(end),(length(Et)-1));
if fig_flag==1
figure
plot(t_arr./1e-15,[Ey;Ez]);
title(title_str)
figure
plot(f(((end-length(E_secondfy)+1):end))./(omega/2/pi),[log(abs(E_secondfy));log(abs(E_secondfz))]);
title(title_str)
savefig([title_str '3D_spec.fig'])
%--------------------------------------
%plot the elliptically polarized light
%--------------------------------------
t_peak=find((abs(Ey)+abs(Ez))==max(abs(Ey)+abs(Ez)));
t_center=t_arr(t_peak)/1e-15;
figure
box on
view(45,30)
%plot projection
hold on
plot3(t_arr./1e-15,(max(abs(Ey))*1.3)*ones(size(Ey)), Ez, 'LineWidth', 1,'Color','k'); % project in x-z axis at y=100
hold on
plot3(t_arr./1e-15,Ey, (min(Ez)*1.3)*ones(size(Ez)), 'LineWidth', 1,'Color','k'); % project in y-z axis at z=-2
hold on
plot3(t_arr./1e-15,real(Ey),Ez,"LineWidth",3,'Color','b')
title(title_str)
xlim([t_center-25*t_scale,t_center+t_scale*30])
savefig([title_str '3Dt.fig'])
end
end