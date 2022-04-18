% clear

N1 = 1e+6;           %��������
N2 = 1e+5;
Dt1 = 1e-3;          %ʱ�䲽�� [s]
Dt2 = 1e-4;
x1 = 0;             %��ʼλ�� [m]
R = 1e-6;           %���Ӱ뾶 [m]
T = 300;            %�¶� [K]
eta = 0.001;        %����ճ���� [Ns/m^2]
kx = 1e-6;          %��ն� [N/m]         %���벽�� [nm]
MaxX = 200;         %x �����ֵ [nm]
Kxy = [0.2 1  5];
SN = [{'(a)'}, {'(b)'}, {'(c)'}];

figure('units','inches','position',[0.5 0.5 6 3])
tiledlayout(1, 2, 'tileSpacing', 'compact', 'Padding', 'compact')

parfor i = 1:length(Kxy)
    [xr(:,i)]=trapped(N1, Dt1, x1, R, T, eta, Kxy(i)*kx);
end
xr = 1e+9.*xr;
parfor i = 1:3
    [r(:,i)]=acf(xr(:,i));
end
[Max, I] = max(r(:,1));
s1 = Dt1*1e+3*([0:1:length(r)-1]-I);

% parfor i = 1:length(Kxy)
%     [xm(:,i)]=trapped(N2, Dt2, x1, R, T, eta, Kxy(i)*kx);
% end
% xm = 1e+9.*xm;
% parfor i = 1:3
%     [msd(:,i)]=MSD(xm(:,i));
% end
% s2 = Dt2*1e+3*[1:1:length(msd)];

nexttile(1)
box on
plot(s1, r(:,1), 'm', s1, r(:,2), 'r', s1, r(:,3), 'k', 'LineWidth', 2)
hold on
axis([-500, 500, -0.05, 1.1])
xticks(-500:250:500)
yticks(0:1)
xlabel('t [ms]', 'FontSize', 16)
ylabel('C_x(t) [a.u.]', 'FontSize', 16)
text(0.01, 0.97, SN(1), 'FontSize', 14, 'Unit', 'normalized')

nexttile(2)
box on
for i = 1:10
    parfor i = 1:length(Kxy)
        [xm(:,i)]=trapped(N2, Dt2, x1, R, T, eta, Kxy(i)*kx);
    end
    xm = 1e+9.*xm;
    parfor i = 1:3
        [msd(:,i)]=MSD(xm(:,i));
    end
    s2 = Dt2*1e+3*[0:1:length(msd)-1];

    plot(s2, msd(:,1), 'm', s2, msd(:,2), 'r', s2, msd(:,3), 'k', 'LineWidth', 2)
    hold on
    if i == 10
        legend('k = 0.2 fN/nm','k = 1 fN/nm', 'k = 5 fN/nm', 'Location','southeast')
    end
end



% plot(s2, msd(:,1), 'm', s2, msd(:,2), 'r', s2, msd(:,3), 'k', 'LineWidth', 1.5)
% hold on
set(gca, 'Xscale', 'log', 'YScale', 'log')
axis([1e-2, 1e+2, 1e+0, 1e+5])
xticks(1e-2*power(10, 0:4))
yticks(1e+0*power(10, 0:5))
xlabel('t [ms]', 'FontSize', 16)
ylabel('x(t)^2 [nm^2]', 'FontSize', 16)
text(0.01, 0.97, SN(2), 'FontSize', 14, 'Unit', 'normalized')
% legend('k = 0.2 fN/nm','k = 1 fN/nm', 'k = 5 fN/nm', 'Location','southeast')


% ���������ӵĲ����˶�
function [x, y, z]=trapped(N, Dt, x1, R, T, eta, kx)
    kB = 1.38e-23;      %������������ [J/K]
    gamma = 6*pi*R*eta; %Ħ��ϵ��
    D = kB*T/gamma;     %��ɢϵ��

    x(1)=x1;   %��ʼ����
    for i = 2:1:N 
        x(i) = x(i-1) - kx*Dt/gamma*x(i-1) + sqrt(2*D*Dt)*randn();
    end
end
%�ٶ�����غ���
function [r]=acf(x)
    r = xcorr(x, ceil(sqrt(length(x))), 'normalized');
end
%����λ��
function [msd]=MSD(x)
    for n = 0:1:ceil(length(x))/50
        msd(n+1) = mean((x(n+1:end)-x(1:end-n)).^2);
    end
end