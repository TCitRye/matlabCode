%Diffusive regime
%This function impolements Eq.(10). Inputs: number of samples N [1e+5], timestep Dt [1e-3 s], initial position x1 [0 m], particle radius R [1e-6 m], temperature T [300 K], and fluid viscosity eta [0.001 pa s - water]. Outputs: particle position x [m] and t [s].

function [x, t]=diffusive(N, Dt, X1, R, T, eta)
    kB = 1.38e-23;      %Boltzmann constant [J/K]
    gamma = 6*pi*R*eta; %friction coefficent
    D = kB*T/gamma;     %diffusion coefficient

    x(1)=x1;            %initial condition
    for i = 2:1:N 
        x(i) = x(i-1) + sqrt(2*D*Dt)*randn();
    end
    t = [0:Dt:(N-1)*Dt];

    %figure
    plot(t, x)
    xlabel('time [s]')
    ylabel('x [m]')
end