# MatLab codes

Giorgio Volpe and Giovanni Volpe

​	Below, wo provide the key MatLab codes to perform the simulations presented in the article. These codes can be modified by the readers themselves to explore more complex cases. In each subsection, we present a function that implements one of the FDE presented above and we suggest some parameters to run the codes. Note that we use thoroughly SI measurements units. 

## Inertial regime

this function implements Eq.(8). Inputs: number of samples N [1e+5], timestep Dt [1e-8], initial positions x1, x2 [0 m, 0 m], paticle radius R [1e-6 m], temperature T [300 K], fluid viscosity eta [0.001 Pa s - water], and density d [2.6e+3 kg/m3], Outputs: particle position x [m] and time t [s].



```matlab
function [x, t]=inertial(N, Dt, x1, x2, R, T, eta, d)
    kB = 1.38e-23;    %Boltzmann constant [J/K]
    gamma = 6*pi*R*eta; %friction coefficient
    m = 4/3*pi*R^3*d;  %particle mass
    tau = m/gamma;    %momnentum relaxation time
    
    x(1)=x1; x2(1)=x2;  %initial conditions
    for i = 3:1:N 
        x(i) = (2+Dt*gamma/m)/(1+Dt*gamma/m)*x(i-1) - 1/(1+Dt*gamma/m)*x(i-2)+sqrt(2*kB*T*gamma)/(m+Dt*gamma)*randn();
    end
    t = [0:Dt:(N-1)*Dt];
    
    figure
    plot(t, x)
    xlabel('time [s]')
    ylabel('x [m]')
end
```



## Diffusive regime

This function implements Eq.(10). Inputs: number of samples N [1e+5], timestep Dt [1e-3 s], initial position x1 [0 m], particle radius R [1e-6 m], temperature T [300 K], and fluid viscosity eta [0.001 pa s - water]. Outputs: particle position x [m] and t [s].

```matlab
function [x, t]=diffusive(N, Dt, X1, R, T, eta)
    kB = 1.38e-23;      %Boltzmann constant [J/K]
    gamma = 6*pi*R*eta; %friction coefficent
    D = kB*T/gamma;     %diffusion coefficient

    x(1)=x1;            %initial condition
    for i = 2:1:N 
        x(i) = x(i-1) + sqrt(2*D*Dt)*randn();
    end
    t = [0:Dt:(N-1)*Dt];

    figure
    plot(t, x)
    xlabel('time [s]')
    ylabel('x [m]')
end
```

## Velocity autocorrelation function

Tis function implements Eq.(12). Inputs: particle position x [m] and timestep Dt [s]. Outpust: velocity autocorrelation function x [m] and delay time s [s].

```matlab
function [r, s]=acf(x, Dt)
    v = (x(2:end)-x(1:end-1))/Dt;    %average speed
    r = xcorr(v, ceil(sqrt(length(x))), 'unbiased');
    s = Dt*[0:1:length(r)-1];

    figure
    plot(s, r)
    xlabel('s [s]')
    ylabel('R [m^2]')
end
```

## Mean square displacement

This function implements Eq.(14). Inputs: position x [m] and timestep Dt [s]. Outputs: mean square displacement med [m2] and delay time s [s]

```matlab
function [msd, s]=MSD(x, Dt)
    for n = 0:1:sqrt(length(x))
        msd(n+1) = mean((x(n+1:end)-x(1:end-n)).^2);
    end
    s = Dt*[0:1:length(med)-1];

    figure
    plot(s, msd)
    xlabel('s [s]')
    ylabel('MSD [m^2]')
end
```

## Trapped particle

This function implements Eq.(16). Inputs: number of sampples N [1e+5], timestep Dt [1e-3 s], initial position x1, x2, x3 [0 m, 0 m, 0 m], particle radius R [1e-6 m],temperature T [300 K], flufid viscosity eta [0.001 Pa s - water], and trap stiffness kx, ky, kz [1e-6 N/m, 1e-6 N/m, 0.2e-6 N/m]. Outputs: particle position x, y, z [m] and time t [s].

```matlab
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
```

