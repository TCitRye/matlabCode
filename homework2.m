% fzero 求某点附近或某区间上的零点
% sym(x)
% y = @sin
% x0 = 2
% x = fzero(y, x0)
% x1 = [2, 10]
% x = fzero(y, x1)

% roots 以向量形式返回以向量表示的多项式的根
% y = [1 -4 -6 -16 4]
% r = roots(y)

% % fsolve ，求解方程组给出初始值最近的零点
% y=@sin
% x0= [1 4]
% x =fsolve(y,x0)

% % solve 求解方程的根
% clear all
% syms x y
% q=4*x+3*y==4
% w=solve(q,x)

% linsolve 求解线性方程组AX = B
% A=[1 2 3; 2 3 4; 3 4 5]
% B=[3 2 1; 3 2 1; 3 2 1]
% X=linsolve(A,B)

% fminbnd 求解函数定区间上最小值
% clear all
% double x
% y = @(x)x^4 - 4 * x^3 - 6 * x^2 - 16 * x + 4
% x=fminbnd(y,-1,4)

% clear all
% x = -1:1:1;
% y = x.^4 - 4 * x.^3 - 6 * x.^2 - 16 * x + 4;
% plot(x, y)

% clear all
% sym x;
% y = @(x)x^4 - 4 * x^3 - 6 * x^2 - 16 * x + 4;
% x0=[-1 4];
% fzero(y,x0)
% x = fminbnd(y, -1, 4)
% y=x^4 - 4 * x^3 - 6 * x^2 - 16 * x + 4

% 结果：
% >> homework2

% ans =

%     0.2278


% x =

%     4.0000


% y =

%  -156.0000