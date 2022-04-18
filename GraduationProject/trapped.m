%Trapped particle
%This function implements Eq.(16). Inputs: number of sampples N [1e+5], timestep Dt [1e-3 s], initial position x1, x2, x3 [0 m, 0 m, 0 m], particle radius R [1e-6 m],temperature T [300 K], flufid viscosity eta [0.001 Pa s - water], and trap stiffness kx, ky, kz [1e-6 N/m, 1e-6 N/m, 0.2e-6 N/m]. Outputs: particle position x, y, z [m] and time t [s].

function [x, y, z, t]=trapped(N, Dt, x1, y1, z1, R, T, eta, kx, ky, kz)
    kB = 1.38e-23;      %Boltzmann constant [J/K]
    gamma = 6*pi*R*eta; %friction coefficient
    D = kB*T/gamma;     %diffusion coefficient

    x(1)=x1; y(1)=y1; z(1)=z1;     %initial condition
    for i = 2:1:N 
        %Deterministic step
        x(i) = x(i-1) - kx*Dt/gamma*x(i-1);
        y(i) = y(i-1) - ky*Dt/gamma*y(i-1);
        z(i) = z(i-1) - kz*Dt/gamma*z(i-1);
        %Diffusive step
        x(i) = x(i) + sqrt(2*D*Dt)*randn();
        y(i) = y(i) + sqrt(2*D*Dt)*randn();
        z(i) = z(i) + sqrt(2*D*Dt)*randn();
    end
    t = [0:Dt:(N-1)*Dt];

    figure
    polt3(x, y, z)
    xlabel('x [m]')
    ylabel('y [m]')
    zlabel('z [m]')
    axis qeual
end