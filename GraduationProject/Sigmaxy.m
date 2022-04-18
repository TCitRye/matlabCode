clear

N = 1e+7;           %样本数量
Dt = 1e-3;          %时间步长 [s]
x1 = 0;             %初始位置 [m]
y1 = 0;
z1 = 0;
R = 1e-6;           %粒子半径 [m]
T = 300;            %温度 [K]
eta = 0.001;        %流体粘滞�? [Ns/m^2]
kx = 1e-6;          %阱刚�? [N/m]
ky = 1e-6;
kz = 0.2e-6;
Dp = 2;             %距离步长 [nm]
MaxX = 300;         %x 轴最大�?? [nm]
MaxY = 300;         %y 轴最大�?? [nm]
Kxy = [0.2 0.3 0.4 0.5 1 2 3 4 5 6 7 8 9];
SN = [{'(a)'}, {'(b)'}, {'(c)'}, {'(D'}];

figure('units','inches','position',[0.5 0.5 6 6])
tiledlayout(3, 3, 'tileSpacing', 'compact', 'Padding', 'compact')

for i = 1:size(Kxy)
    [x(:,i), y(:,i), z(:,i), t(i)]=trapped(N, Dt, x1, y1, z1, R, T, eta, Kxy(i)*kx, Kxy(i)*ky, kz);
end
x = 1e+9.*x;
y = 1e+9.*y;
z = 1e+9.*z;
sigmax = std(x);
sigmay = std(y);
fun = @(a, xdata)a(1)/xdata
a = lsqcurvefit(fun, [1], Kxy, sigmax.*sigmay);
[xyPlane(:,:,1)] = ProbabilityDensity(x(:,1), y(:,1), z(:,1), Dp, MaxX, MaxY);
[xyPlane(:,:,2)] = ProbabilityDensity(x(:,5), y(:,5), z(:,5), Dp, MaxX, MaxY);
[xyPlane(:,:,3)] = ProbabilityDensity(x(:,9), y(:,9), z(:,9), Dp, MaxX, MaxY);

nexttile([2 3]);
plot(Kxy, sigmax.*sigmay, 'ko', [0:0.1:10],fun(x,times),'b-')
hold on
axis([0, 10, 0, 5])

nexttile([1 1]);
box on
imagesc([-MaxX MaxX], [-MaxY MaxY], xyPlane(:,:,1));
colormap(flipud(gray))
axis([-300, 300, -300, 300])
yticks(-200:200:200)
xticks(-200:200:200)
xlabel('x [nm]', 'FontSize', 16)
ylabel('y [nm]', 'FontSize', 16)
set(gca, 'yaxislocation', 'right', 'FontSize', 12, 'TickLength', [0 0])
axis xy
pbaspect([1 1 1]) 
hold on
text(0.01, 0.93, SN(2), 'FontSize', 14, 'Unit', 'normalized')

nexttile([1 1]);
box on
imagesc([-MaxX MaxX], [-MaxY MaxY], xyPlane(:,:,2));
colormap(flipud(gray))
axis([-300, 300, -300, 300])
yticks(-200:200:200)
xticks(-200:200:200)
xlabel('x [nm]', 'FontSize', 16)
ylabel('y [nm]', 'FontSize', 16)
set(gca, 'yaxislocation', 'right', 'FontSize', 12, 'TickLength', [0 0])
axis xy
pbaspect([1 1 1]) 
hold on
text(0.01, 0.93, SN(3), 'FontSize', 14, 'Unit', 'normalized')

nexttile([1 1]);
box on
imagesc([-MaxX MaxX], [-MaxY MaxY], xyPlane(:,:,3));
colormap(flipud(gray))
axis([-300, 300, -300, 300])
yticks(-200:200:200)
xticks(-200:200:200)
xlabel('x [nm]', 'FontSize', 16)
ylabel('y [nm]', 'FontSize', 16)
set(gca, 'yaxislocation', 'right', 'FontSize', 12, 'TickLength', [0 0])
axis xy
pbaspect([1 1 1]) 
hold on
text(0.01, 0.93, SN(4), 'FontSize', 14, 'Unit', 'normalized')





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

%粒子在某�?平面的概率密度投�?
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