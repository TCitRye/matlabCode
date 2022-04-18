clear

N = 1e+5;           %样本数量
Dt = 1e-3;          %时间步长 [s]
x1 = 0;             %初始位置 [m]
y1 = 0;
z1 = 0;
R = 1e-6;           %粒子半径 [m]
T = 300;            %温度 [K]
eta = 0.001;        %流体粘滞性 [Ns/m^2]
kx = 1e-6;          %阱刚度 [N/m]
ky = 1e-6;
kz = 0.2e-6;
Dp = 10;             %距离步长 [nm]
MaxX = 200;         %x 轴最大值 [nm]
MaxY = 200;         %y 轴最大值 [nm]
MaxZ = 400;         %z 轴最大值 [nm]

figure('units','inches','position',[0.5 0.5 6 6])
tiledlayout(3, 3, 'tileSpacing', 'compact', 'Padding', 'compact')

nexttile([3 2])
box on
[x, y, z, t]=trapped(N, Dt, x1, y1, z1, R, T, eta, kx, ky, kz);
x = 1e+9*x;
y = 1e+9*y;
z = 1e+9*z;
plot3(x, y, z)
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
% view(90,0)
axis equal

hold on
nexttile([2 1])
box on

for i = -MaxX:Dp:MaxX
    x = ceil(x/Dp)*Dp;
    y = ceil(y/Dp)*Dp;
    indexy = find(ceil(x) == i);
    for j = -MaxY:Dp:MaxY
        indexz = find(ceil(y(indexy)) == j);
        Z(i+MaxX+1:i+MaxX+Dp, j+MaxY+1:j+MaxY+Dp) = length(z(indexz));
    end
end

imagesc([-MaxX MaxX], [-MaxY MaxY], Z);
colormap(flipud(gray))
axis equal
hold on






% 光阱中粒子的布朗运动
function [x, y, z, t]=trapped(N, Dt, x1, y1, z1, R, T, eta, kx, ky, kz)
    kB = 1.38e-23;      %波尔兹曼常数 [J/K]
    gamma = 6*pi*R*eta; %摩擦系数
    D = kB*T/gamma;     %扩散系数

    x(1)=x1; y(1)=y1; z(1)=z1;     %初始条件
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
end


