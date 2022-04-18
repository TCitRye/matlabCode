clear
clc

N = 1e+7;        %��������
Dt = 1e-8;       %���� [s]
x1 = 0;          %��ʼλ�� [m]
x2 = 0;          
R = 1e-6;        %���Ӱ뾶 [m]
T = 300;         %�¶� [K]
eta = 0.001;     %����ճ���� [Ns/m^2]
d = 2.6e+3;      %�ܶ� [kg/m^3]
SN = [{'(a)'}, {'(b)'}, {'(c)'}, {'(d)'}];


figure('units','inches','position',[0.5 0.5 6 6])
tiledlayout(2, 2, 'tileSpacing', 'compact', 'Padding', 'compact')

[xi, xn, t]=inertial(N, Dt, x1, x2, R, T, eta, d);    %��ʵ���ӵĲ����˶�:�й����� xi���޹����� xn
[r, s]=acf(xi, Dt, R, eta, d);      %�ٶ������
[msdi, si]=MSD(xi, Dt, R, eta, d);  %����λ���й����� msdi���޹����� msdn
[msdn, sn]=MSD(xn, Dt, R, eta, d);


nexttile(1)
box on
plot(t, xi*1e9, 'r', t, xn*1e9, 'k--')
hold on
axis([0, 1, -0.5, 0.5])
xticks(0:0.5:1)
yticks(-0.5:0.5:0.5)
xlabel('t/\tau', 'FontSize', 16)
ylabel('x [nm]', 'FontSize', 16)
text(0.01, 0.97, SN(1), 'FontSize', 14, 'Unit', 'normalized')


nexttile(2)
box on
plot(t, xi*1e9, 'r', t, xn*1e9, 'k--')
hold on
axis([0, 100, -5, 5])
xticks(0:50:100)
yticks(-5:5:5)
xlabel('t/\tau', 'FontSize', 16)
ylabel('x [nm]', 'FontSize', 16)
text(0.01, 0.97, SN(2), 'FontSize', 14, 'Unit', 'normalized')


nexttile(3)
box on
plot(s, r, 'r', [0, 0], [0, 1], 'k--')
hold on
axis([-8, 8, 0, 1.1])
xticks(-6:3:6)
yticks(0:1)
xlabel('t/\tau', 'FontSize', 16)
ylabel('C_v(t) [a.u.]', 'FontSize', 16)
text(0.01, 0.97, SN(3), 'FontSize', 14, 'Unit', 'normalized')


nexttile(4)
box on
loglog(si, msdi*1e18, 'r', sn, msdn*1e18, 'k--')
hold on
axis([-5e-7, 50, 1e-5, 100])
xticks(0.1*power(10, 0:2))
yticks(1e-6*power(10, 0:2:8))
% get(gca,'xtick')  % �õ������ʵ�ʴ�С
set(gca,'xticklabel',get(gca,'xtick')) % �� x ��ʾ���ַ��滻Ϊʵ�ʴ�С
xlabel('t/\tau', 'FontSize', 16)
ylabel('x(t)^2 [nm^2]', 'FontSize', 16)
text(0.01, 0.97, SN(4), 'FontSize', 14, 'Unit', 'normalized')
legend('inertial','non-inertial', 'Location','southeast')



%��ʵ���ӵĲ����˶�
function [xi, xn, t]=inertial(N, Dt, x1, x2, R, T, eta, d)
    kB = 1.38e-23;      %������������ [J/K]
    gamma = 6*pi*R*eta; %Ħ��ϵ��
    m = 4/3*pi*R^3*d;   %���ӵ�����
    tau = m/gamma;      %������ԥʱ��

    xi(1) = x1; xi(2) = x2;  %��ʼ����
    xn(1) = x1; xn(2) = x2;
    for i = 3:1:N 
        w(i) = randn();
        %�й��������
        xi(i) = (2+Dt*gamma/m)/(1+Dt*gamma/m)*xi(i-1) - 1/(1+Dt*gamma/m)*xi(i-2) + sqrt(2*kB*T*gamma)/(m+Dt*gamma)*Dt^(3/2)*w(i);
        %�޹��������
        xn(i) = xn(i-1) + sqrt(2*kB*T/gamma*Dt)*w(i);
    end

    t = [0:Dt:(N-1)*Dt]./tau;
end
%�ٶ�����غ���
function [r, s]=acf(x, Dt, R, eta, d)
    gamma = 6*pi*R*eta; %Ħ��ϵ��
    m = 4/3*pi*R^3*d;   %���ӵ�����
    tau = m/gamma;      %������ԥʱ��

    v = (x(2:end)-x(1:end-1))/Dt;    %ƽ���ٶ�
    r = xcorr(v, ceil(sqrt(length(x))), 'normalized');
    [Max, I] = max(r);
    s = Dt*([0:1:length(r)-1]-I)./tau;
end
%����λ��
function [msd, s]=MSD(x, Dt, R, eta, d)
    gamma = 6*pi*R*eta; %Ħ��ϵ��
    m = 4/3*pi*R^3*d;   %���ӵ�����
    tau = m/gamma;      %������ԥʱ��
    for n = 0:1:sqrt(length(x))
        msd(n+1) = mean((x(n+1:end)-x(1:end-n)).^2);
    end
    s = Dt*[0:1:length(msd)-1]./tau;
end