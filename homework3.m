% x = 0:10; y = sin(x);
% xi = 0:.15:10;
% yi1 = interp1(x, y, xi, 'nearest'); %最邻近插值
% yi2 = interp1(x, y, xi, 'linear'); %线性插值
% yi3 = interp1(x, y, xi, 'spline'); %分段三次样条插值
% yi4 = interp1(x, y, xi, 'cubic'); %保型分段三次插值
% x = 0:0.1:10; y = sin(x);
% subplot(2, 2, 1); plot(x, y, x, y, 'b', xi, yi1, xi, yi1, 'r:.');
% subplot(2, 2, 2); plot(x, y, x, y, 'b', xi, yi2, xi, yi2, 'r:.');
% subplot(2, 2, 3); plot(x, y, x, y, 'b', xi, yi3, xi, yi3, 'r:.');
% subplot(2, 2, 4); plot(x, y, x, y, 'b', xi, yi4, xi, yi4, 'r:.');

x = 0:2 * pi; y = 1 - sin(x);
xi = 0:pi / 20:2 * pi;
yi1 = interp1(x, y, xi, 'nearest'); %最邻近插值
yi2 = interp1(x, y, xi, 'linear'); %线性插值
yi3 = interp1(x, y, xi, 'spline'); %分段三次样条插值
yi4 = interp1(x, y, xi, 'pchip'); %保型分段三次插值
x = 0:pi / 30:2 * pi; y = 1 - sin(x);
subplot(2, 2, 1); polarplot(x, y, x, y, 'b', xi, yi1, xi, yi1, 'r:.');
subplot(2, 2, 2); polarplot(x, y, x, y, 'b', xi, yi2, xi, yi2, 'r:.');
subplot(2, 2, 3); polarplot(x, y, x, y, 'b', xi, yi3, xi, yi3, 'r:.');
subplot(2, 2, 4); polarplot(x, y, x, y, 'b', xi, yi4, xi, yi4, 'r:.');
