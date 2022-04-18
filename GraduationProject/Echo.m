clear

t = 30;
Dt = 0.1;
m = 10;

figure('units','inches','position',[0.5 0.5 9 6])
tiledlayout(2, 3, 'tileSpacing', 'compact', 'Padding', 'compact')

[x, t] = tra(Dt, t, m);
r = xcorr(x, 'normalized');
[Max, I] = max(r);
t = [0:Dt:ceil(2*t/Dt)*Dt-Dt]-I*Dt;
nexttile(1)
plot(t(I:end-1), x)
hold on
nexttile(2)
plot(t(1:end-1), r)
hold on




function [x, t]=tra(Dt, t, m)
    for i=1:ceil(t/Dt)
        x1(i) =sin(0.1*i)*randn();
        if i > ceil(m/Dt)
            x(i) =x1(i) + x1(i-ceil(m/Dt));
        else
           x(i) = x1(i);  
        end
    end
end
