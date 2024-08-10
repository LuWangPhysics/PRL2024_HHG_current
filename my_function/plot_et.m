function plot_et(E,J,C)
% figure
% plot(E.t./(1e-15),[jinter;jintra])
% legend('inter','intra')
% title('inter/intra')
figure;
plot(E.t./(1e-15),J{1})
hold on
plot(E.t./(1e-15),J{2})
hold on
plot(E.t./(1e-15),J{3})
legend("x","y","z")
ylabel('J(t)')
xlabel('t (fs)')

%eta0=C.e*C.v(1)*E0*2*pi/C.omega^2/C.hbar;
end