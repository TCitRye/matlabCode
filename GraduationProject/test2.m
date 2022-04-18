clear 
Dt = 1;
T = 30;
m = 1000;
x(1:m,1) = 0;
for j=1:m
    for i=1:ceil(T/Dt)
        x(j,i+1) = x(j,i) + sqrt(Dt)*randn();
    end
end

for i=1:ceil(T/Dt)+1
    a = 0;
    b = 0;
    xp = [];
    xn = [];
    a = find(x(:, i) >= 0);
    xp = x(a,i);
    b = find(x(:, i) <= 0);
    xn = x(b,i);
    Stdx(i, 1) = std(xp);
    Stdx(i, 2) = -std(xn);


end
t = [0:Dt:ceil(T/Dt)*Dt];