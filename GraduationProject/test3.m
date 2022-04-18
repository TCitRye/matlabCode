clear

N = 1e+5;           %��������
Dt = 1e-3;          %ʱ�䲽�� [s]
x1 = 0;             %��ʼλ�� [m]
y1 = 0;
z1 = 0;
R = 1e-6;           %���Ӱ뾶 [m]
T = 300;            %�¶� [K]
eta = 0.001;        %����ճ���� [Ns/m^2]
kx = 1e-6;          %��ն� [N/m]
ky = 1e-6;
kz = 0.2e-6;
Dp = 10;             %���벽�� [nm]
MaxX = 200;         %x �����ֵ [nm]
MaxY = 200;         %y �����ֵ [nm]
MaxZ = 400;         %z �����ֵ [nm]

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






% ���������ӵĲ����˶�
function [x, y, z, t]=trapped(N, Dt, x1, y1, z1, R, T, eta, kx, ky, kz)
    kB = 1.38e-23;      %������������ [J/K]
    gamma = 6*pi*R*eta; %Ħ��ϵ��
    D = kB*T/gamma;     %��ɢϵ��

    x(1)=x1; y(1)=y1; z(1)=z1;     %��ʼ����
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


