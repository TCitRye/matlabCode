x = [1:0.01:1000];

y = sqrt(x).*sin(x.*pi);

figure

plot(x, y, 'LineWidth', 2.5)
axis([0, 1000, -1000, 1000])