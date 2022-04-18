%Mean square displacement
%This function implements Eq.(14). Inputs: position x [m] and timestep Dt [s]. Outputs: mean square displacement med [m2] and delay time s [s]

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