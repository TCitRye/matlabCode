clear
Dt = 1;
N = 30;
m = 10000;

[x, Stdx]=Tra(Dt, N, m);
Sum=Test(10)

function [x, Stdx]=Tra(Dt, N, m)
    x(1:m,1) = 0;
    for j=1:m
        for i=1:ceil(N/Dt)
            x(j,i+1) = x(j,i) + sqrt(Dt)*randn();
        end
    end

    t2 = [0:Dt:ceil(N/Dt)*Dt];

    plot(t2, x(1,:))
    hold on;

    for i=1:ceil(N/Dt)+1
        Stdx(i) = std(x(:, i));
    end

    plot(t2, Stdx)
    hold on;
end

function Sum=Test(n)
    Sum = 0;
    for i=1:n
        Sum = Sum + i;
    end
end

