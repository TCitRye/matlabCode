clc
clear
pic_num = 1;
t = 0:0.01:2 * pi;
x = cos(t);
y = sin(t);

plot(x, y);
axis equal;
hold on;

lineX1 = [0, 0];
lineY1 = [0, 1];
lineX2 = [0, 0];
lineY2 = [0, 1];
h1 = plot(lineX1, lineY1);
h2 = plot(lineX2, lineY2);

theta = 0;

for i =     1:382 * pi
    theta = theta + pi / 10;
    lineX1(2) = sin(theta / 60);
    lineY1(2) = cos(theta / 60);
    lineX2(2) = sin(theta);
    lineY2(2) = cos(theta);
    set(h1, 'XData', lineX1, 'Ydata', lineY1);
    set(h2, 'XData', lineX2, 'Ydata', lineY2);
    drawnow;

    F = getframe(gcf);
    I = frame2im(F);
    [I, map] = rgb2ind(I, 256);

    if pic_num == 1
        imwrite(I, map, 'test.gif', 'gif', 'Loopcount', inf, 'DelayTime', 0.2);
    else
        imwrite(I, map, 'test.gif', 'gif', 'WriteMode', 'append', 'DelayTime', 0.2);
    end

    pic_num = pic_num + 1;
end
