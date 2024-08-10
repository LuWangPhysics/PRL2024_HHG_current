function entire_field_t(Ey,Ez,tau,t_arr,title_str)
if tau<20e-15
t_scale=1;
else
 t_scale=12;
end
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
xlim([t_center-30*t_scale,t_center+t_scale*30])
savefig([title_str '3Dt.fig'])
end