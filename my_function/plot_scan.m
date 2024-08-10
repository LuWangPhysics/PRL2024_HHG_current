function plot_scan(t,cep_loop,C,f,Eyw,Ezw,Eyt,Ezt,folder_path)
if (length(cep_loop)>1)
decay_time=C.tau*1e15;
figure
surf(f./(C.omega/2/pi),cep_loop,log(abs(Eyw)),'LineStyle','none')
xlabel('harmonic order')
ylabel('cep')
zlabel('Jyw')
xlim([0,10])
%savefig([folder_path 'Jyw_decay' num2str(decay_time) 'fs.fig'])
figure
surf(f./(C.omega/2/pi),cep_loop,log(abs(Ezw)),'LineStyle','none')
xlabel('harmonic order')
ylabel('cep')
zlabel('Jzw')
xlim([0,10])
%savefig([folder_path 'Jzw_decay' num2str(decay_time) 'fs.fig'])
figure
surf(t./1e-15,cep_loop,Ezt,'LineStyle','none')
xlabel('t(fs)')
ylabel('cep')
zlabel('Jz')
%savefig([folder_path 'Jz_decay' num2str(decay_time) 'fs.fig'])
figure
surf(t./1e-15,cep_loop,Eyt,'LineStyle','none')
xlabel('t(fs)')
ylabel('cep')
zlabel('Jy')
end
end