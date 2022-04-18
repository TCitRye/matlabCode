%速度自相关函数
%Tis function implements Eq.(12). Inputs: particle position x [m] and timestep Dt [s]. Outpust: velocity autocorrelation function x [m] and delay time s [s].

function [r, s]=acf(x, Dt)
    v = (x(2:end)-x(1:end-1))/Dt;    %average speed
    r = xcorr(v, ceil(sqrt(length(x))), 'unbiased');
    s = Dt*[0:1:length(r)-1];

    figure
    plot(s, r)
    xlabel('s [s]')
    ylabel('R [m^2]')
end