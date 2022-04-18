% clear

N = 1e+7;           %样本数量
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
Dp = 2;             %距离步长 [nm]
MaxX = 200;         %x 轴最大值 [nm]
MaxY = 200;         %y 轴最大值 [nm]
MaxZ = 400;         %z 轴最大值 [nm]
SN = [{'(a)'}, {'(b)'}, {'(c)'}];

figure('units','inches','position',[0.5 0.5 6 6])
tiledlayout(13, 13, 'tileSpacing', 'compact', 'Padding', 'compact')


% [x, y, z, t]=trapped(N, Dt, x1, y1, z1, R, T, eta, kx, ky, kz);
% x = 1e+9*x;
% y = 1e+9*y;
% z = 1e+9*z;
% [xzPlane] = ProbabilityDensity(x, z, y, Dp, MaxX, MaxZ);
% [xyPlane] = ProbabilityDensity(x, y, z, Dp, MaxX, MaxY);
% sigmax = std(x);
% sigmay = std(y);
% sigmaz = std(z);
% [X, Y, Z] = sphere;
% X = 3*sigmax*X;
% Y = 3*sigmay*Y;
% Z = 3*sigmaz*Z;
% C(:, :, 1) = 0.5*ones(size(X));
% C(:, :, 3) = 0.5*ones(size(Y));
% C(:, :, 2) = 0.5*ones(size(Z));

nexttile([13 8]);
surf(X, Y, Z, C, 'EdgeColor', 'none', 'FaceAlpha', 0.3)
hold on
plot3(x(1:1e3), y(1:1e3), z(1:1e3), 'LineWidth', 1.5)
hold on
axis([-200, 200, -200, 200, -500, 500])
xticks(-200:200:200)
yticks(-200:200:200)
zticks(-400:200:400)
xlabel('x [nm]', 'FontSize', 16, 'Position', [0 -300 -500])
ylabel('y [nm]', 'FontSize', 16, 'Position', [300 0 -500])
zlabel('z [nm]', 'FontSize', 16)
set(gca,'FontSize', 12)
view(60, 15)
grid on
box on
ax = gca;
ax.BoxStyle = 'full';
text(0.01, 0.97, SN(1), 'FontSize', 14, 'Unit', 'normalized')


nexttile([9 5])
box on
imagesc([-MaxX MaxX], [MaxZ -MaxZ], xzPlane);
colormap(flipud(gray))
yticks(-400:200:400)
ylabel('z [nm]', 'FontSize', 16)
set(gca, 'xticklabel', [], 'FontSize', 12, 'TickLength', [0 0])
set(gca, 'yaxislocation', 'right')
axis xy
hold on
axis equal
text(0.01, 0.97, SN(2), 'FontSize', 14, 'Unit', 'normalized')

nexttile([4 5])
box on
imagesc([-MaxX MaxX], [-MaxY MaxY], xyPlane);
colormap(flipud(gray))
axis([-200, 200, -200, 200])
yticks(-200:200:200)
xticks(-200:200:200)
xlabel('x [nm]', 'FontSize', 16)
ylabel('y [nm]', 'FontSize', 16)
set(gca, 'yaxislocation', 'right', 'FontSize', 12, 'TickLength', [0 0])
axis xy
pbaspect([1 1 1]) 
hold on
text(0.01, 0.93, SN(3), 'FontSize', 14, 'Unit', 'normalized')





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

%粒子在某一平面的概率密度投影
function [xyPlane]=ProbabilityDensity(x, y, z, Dp, MaxX, MaxY)
    for i = -MaxX:Dp:MaxX
        x = ceil(x/Dp)*Dp;
        y = ceil(y/Dp)*Dp;
        indexy = find(ceil(x) == i);
        for j = -MaxY:Dp:MaxY
            indexz = find(ceil(y(indexy)) == j);
            xyPlane(i+MaxX+1:i+MaxX+Dp, j+MaxY+1:j+MaxY+Dp) = length(z(indexz));
        end
    end
end