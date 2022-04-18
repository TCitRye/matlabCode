clear

N1 = 1e+6;           %样本数量
N2 = 1e+7;
Dt = 1e-4;          %时间步长 [s]
x1 = 0;             %初始位置 [m]
R = 1e-6;           %粒子半径 [m]
T = 300;            %温度 [K]
eta = 0.001;        %流体粘滞性 [Ns/m^2]
kx = 1e-6;          %阱刚度 [N/m]
ky = 1e-6;
kz = 0.2e-6;
Dp = 2;             %距离步长 [nm]
MaxX = 200;         %x 轴最大值 [nm]
MaxY = 200;         %y 轴最大值 [nm]
MaxZ = 400;         %z 轴最大值 [nm]
Kxy = [0.2 1  5];
SN = [{'(a)'}, {'(b)'}, {'(c)'}];


figure('units','inches','position',[0.5 0.5 6 3])
tiledlayout(1, 2, 'tileSpacing', 'compact', 'Padding', 'compact')

[xr]=trapped(N1, Dt, x1, R, T, eta, Kxy(3)*kx);
xr = 1e+9.*xr;
    [r]=acf(xr);
[Max, I] = max(r);
s1 = Dt*1e+3*([0:1:length(r)-1]-I);
s2 = Dt*1e+3*[1:100];


nexttile(1)
box on
plot(s1, r, 'b')
hold on
axis([-500, 500, 0, 1.1])
xticks(-500:250:500)
yticks(0.01:1.1)
xlabel('t [ms]', 'FontSize', 16)
ylabel('C_v(t) [a.u.]', 'FontSize', 16)
text(0.01, 0.97, SN(1), 'FontSize', 14, 'Unit', 'normalized')

nexttile(2)
box on
plot(s2, xr(1:100), 'b')
hold on
% xticks(1e-2*power(10, 0:4))
% yticks(1e+0*power(10, 0:5))
xlabel('t [ms]', 'FontSize', 16)
ylabel('x(t)^2 [nm^2]', 'FontSize', 16)
text(0.01, 0.97, SN(2), 'FontSize', 14, 'Unit', 'normalized')
legend('k = 0.2 fN/nm', 'Location','southeast')

% 光阱中粒子的布朗运动
function [x]=trapped(N, Dt, x1, R, T, eta, kx)
    kB = 1.38e-23;      %波尔兹曼常数 [J/K]
    gamma = 6*pi*R*eta; %摩擦系数
    D = kB*T/gamma;     %扩散系数
    x(1)=x1;   %初始条件
    for i = 2:1:N 
        %Deterministic step
        x(i) = x(i-1) - kx*Dt/gamma*x(i-1) + sqrt(2*D*Dt)*randn();
    end
end
%速度自相关函数
function [r, s]=acf(x)
    r = xcorr(x, ceil(sqrt(length(x))), 'normalized');
    % [Max, I] = max(r);
    % s = ([0:1:length(r)-1]-I);
end