clear
sigma  = [1 2];

fun = @(sigma, x, y)1/(sigma(1)*sigma(2)*2*pi)*exp(-1/2*(x.^2+y.^2)/(sigma(1)*sigma(2)));

[x, y] = meshgrid(-10:0.1:10);


surf(x, y, fun(sigma, x, y), 'EdgeColor', 'none')