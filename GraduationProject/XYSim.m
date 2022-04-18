clear

T = 30;       %模拟时间
m = 10000;     %模拟次数

Dt = [1, 0.5, 0.1];

figure('units','inches','position',[0.5 0.5 9 6])

for i=1:3

    [W, t1]=Wi(Dt(i), T);
    [x, Stdx, t2]=Tra(Dt(i), T, m);

    subplot(2, 3, i)
    scatter(t1, W, 5, 'filled')
    hold on;
    axis([0, 30, -8, 8])
    yticks(-4:4:8)
    box on;

    if i > 1
        set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[]);
    end

    subplot(2, 3, 3+i)
    area(t2, Stdx, 'LineStyle', 'none', 'ShowBaseLine', 'off', 'FaceColor', [0.8 0.8 0.8])
    hold on;
    area(t2, -Stdx, 'LineStyle', 'none', 'ShowBaseLine', 'off', 'FaceColor', [0.8 0.8 0.8])
    hold on;
    plot(t2, x(1, :))
    hold on;
    axis([0, 30, -8, 8])
    yticks(-4:4:8)
    box on;
    
    if i > 1
        set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[]);
    end

end

function [W, t]=Wi(Dt, T)

    for i=1:ceil(T/Dt)
        W(i) = randn()/sqrt(Dt);
    end

    t = [Dt:Dt:ceil(T/Dt)*Dt];
end

function [x, Stdx, t]=Tra(Dt, T, m)

    x(1:m,1) = 0;
    for j=1:m
        for i=1:ceil(T/Dt)
            x(j,i+1) = x(j,i) + sqrt(Dt)*randn();
        end
    end

    for i=1:ceil(T/Dt)+1
        Stdx(i) = std(x(:, i));
    end


    t = [0:Dt:ceil(T/Dt)*Dt];
end


