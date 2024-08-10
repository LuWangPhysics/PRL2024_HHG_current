function [test]=bessel_sum(nb,N,t,omega,Omega_0,flag)
%remove the integration at -infinity is important this is the same as
%exluding j0
nb_e=fix(nb/2)*2;
nb_o=nb_e+1;


test=0;
switch flag
    case 'coe_all'
    for n=[-nb:1:nb]
    test=test+besselj(n,N).*exp(1i*(omega.*n+Omega_0).*t).*n./(n+Omega_0/omega);
    end

    case 'plain_odd'
    for n=[-nb_o:2:nb_o]
        test=test+exp(1i*(omega.*n).*t);
    end
    case 'plain_even'
    for n=[-nb_e:2:nb_e]
        test=test+exp(1i*(omega.*n).*t);
    end



end

end