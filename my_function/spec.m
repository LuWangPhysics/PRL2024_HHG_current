function [Jf,f_arr]=spec(t,J,C,flag_plot)
dt=t(2)-t(1);
df=1/dt/length(t);
f_limit=length(t).*df/2;
f_arr=linspace(-f_limit,f_limit,length(t));

N_half=fix(length(t)/2)+1;
f_half=f_arr(N_half:end);

Jf={0,0,0};
if flag_plot==1
figure
end
for i_iter=1:3
temp=fftshift(fft(ifftshift(J{i_iter})))./(length(t)*df);
%normalized the energy to be the same as t domain
Jf{i_iter}=sqrt(2).*temp(N_half:end);
plot_f=log(abs(Jf{i_iter}));
if flag_plot==1

plot(f_half./(C.omega/2/pi),plot_f)
end
hold on
end
xlabel('harmonics')
xlim([0,10])
% set(gca, 'YScale', 'log')
% ylabel('log Jf')
legend("x","y","z")


end


