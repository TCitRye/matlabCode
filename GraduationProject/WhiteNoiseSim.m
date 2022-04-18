clear

t = 30;         %模拟时间 [s]
m = 10000;      %模拟次数
SN = [{'(a)'}, {'(b)'}, {'(c)'}, {'(d)'}, {'(e)'}, {'(f)'}];

Dt = [1, 0.5, 0.1];   %时间步长 [s]

figure('units','inches','position',[0.5 0.5 9 6])
tiledlayout(2, 3, 'tileSpacing', 'compact', 'Padding', 'compact')

for i=1:3

    [W, t1]=Wi(Dt(i), t);               %白噪声模拟：白噪声 W
    [x, Stdx, t2]=tra(Dt(i), t, m);     %粒子轨迹模拟：粒子轨迹 x，粒子位子标准差 Stdx
    %绘制白噪声图像
    nexttile(i)
    scatter(t1, W, 5, 'filled')
    hold on;
    axis([0, 30, -8, 8])
    yticks(-8:4:8)
    text(1, 7, SN(i), 'FontSize', 14)
    text(1, -7, ['\Delta t = ', num2str(Dt(i))], 'FontSize', 14)
    set(gca, 'xtick', [], 'xticklabel', [])
    box on;
    if i == 1
        ylabel('W_i', 'FontSize', 16)
    else
        set(gca, 'ytick', [], 'yticklabel', [])
    end
    %绘制粒子轨迹图像
    nexttile(i+3)
    area(t2, Stdx(:, 1), 'LineStyle', 'none', 'ShowBaseLine', 'off', 'FaceColor', [0.8 0.8 0.8])
    hold on;
    area(t2, Stdx(:, 2), 'LineStyle', 'none', 'ShowBaseLine', 'off', 'FaceColor', [0.8 0.8 0.8])
    hold on;
    plot(t2, x(1, :), 'b')
    hold on;
    axis([0, 30, -8, 8])
    yticks(-8:4:8)
    xlabel('t', 'FontSize', 16)
    text(1, 7, SN(i+3), 'FontSize', 14)
    text(1, -7, ['\Delta t = ', num2str(Dt(i))], 'FontSize', 14)
    box on;
    if i == 1
        ylabel('W_i', 'FontSize', 16)
    else
        set(gca, 'ytick', [], 'yticklabel', [])
    end
end
%白噪声模拟函数
function [W, t]=Wi(Dt, t)
    for i=1:ceil(t/Dt)
        W(i) = randn()/sqrt(Dt);
    end

    t = [Dt:Dt:ceil(t/Dt)*Dt];
end
%粒子轨迹模拟函数
function [x, Stdx, t]=tra(Dt, t, m)
    x(1:m,1) = 0;
    for j=1:m
        for i=1:ceil(t/Dt)
            x(j,i+1) = x(j,i) + sqrt(Dt)*randn();
        end
    end
    %分别计算 xi 大于零和小于零的粒子对 xbar=0 的标准差
    for i=1:ceil(t/Dt)+1
        a = 0;
        b = 0;
        xp = [];
        xn = [];
        a = find(x(:, i) >= 0);
        xp = x(a,i);
        b = find(x(:, i) <= 0);
        xn = x(b,i);
        Stdx(i, 1) = sqrt(sum(xp.^2)/prod(size(xp)));
        Stdx(i, 2) = -sqrt(sum(xn.^2)/prod(size(xn)));
    end

    t = [0:Dt:ceil(t/Dt)*Dt];
end


