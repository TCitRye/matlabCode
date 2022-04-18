%Inertial regime
% this function implements Eq.(8). Inputs: number of samples N [1e+5], timestep Dt [1e-8 s], initial positions x1, x2 [0 m, 0 m], paticle radius R [1e-6 m], temperature T [300 K], fluid viscosity eta [0.001 Pa s - water], and density d [2.6e+3 kg/m3], Outputs: particle position x [m] and time t [s].
%对速度自相关数据进行指数函数拟合
% s1 = s(find(s <= 0));
% s2 = s(find(s >= 0));
% fun1 = @(a, s)a(1)*exp(s*a(2)+a(3));
% a = lsqcurvefit(fun1, [0, 0, 0], s1, r(find(s <= 0)));
% fun2 = @(b, s)b(1)*exp(s*b(2)+b(3));
% b = lsqcurvefit(fun2, [0, 0, 0], s2, r(find(s >= 0)));
% s1 = [-8:0.01:0];
% s2 = [0:0.01:8];

% plot(s1, fun1(a, s1), 'r', s2, fun2(b, s2), 'r', [0, 0], [0, 1], 'k--')
function [x, t]=inertial(N, Dt, x1, x2, R, T, eta, d)
    kB = 1.38e-23;      %Boltzmann constant [J/K]
    gamma = 6*pi*R*eta; %friction coefficient
    m = 4/3*pi*R^3*d;   %particle mass
    tau = m/gamma;      %momnentum relaxation time

    x(1)=x1; x2(1)=x2;  %initial conditions
    for i = 3:1:N 
        x(i) = (2+Dt*gamma/m)/(1+Dt*gamma/m)*x(i-1) - 1/(1+Dt*gamma/m)*x(i-2)+sqrt(2*kB*T*gamma)/(m+Dt*gamma)*Dt^(3/2)*randn();
    end
    t = [0:Dt:(N-1)*Dt];

    %figure
    plot(t, x)
    xlabel('time [s]')
    ylabel('x [m]')
end