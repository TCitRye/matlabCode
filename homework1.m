a = [1, 2, 3; 4, 5, 6; 7, 8, 9]
b = zeros(3, 3)%生成0矩阵
c = ones(3, 3)%生成1矩阵
d = rank(a)%矩阵求秩
e = mtimes(a, c)%矩阵相乘
f = mpower(a, 2)%矩阵的幂
g = hadamard(4)%hadamard矩阵
h = hankel(1:3, 3:5)%hankel矩阵
k = hilb(3)%hilbert矩阵
l = invhilb(3)%hilbert矩阵逆矩阵
m = magic(3)%幻方矩阵
n = pascal(3)%帕斯卡矩阵
o = rosser%rosser矩阵
p = toeplitz(1:3)%托普利茨矩阵
q = vander(1:3)%vandermond矩阵
r = inv(f)%矩阵求逆
s = eig(a)%矩阵特征向量
[L,U,P] = lu(a)%矩阵LU分解
t = transpose(a)%矩形转置
u = sqrtm(f)%矩阵平方根
v = expm(a)%矩阵的指数
w = logm(a)%矩阵的对数
x = cross(a,c)%矩阵叉积
y = dot(a,b)%矩阵点积
z = tril(a)%矩阵下三角部分
A = triu(a)%矩阵上三角
B = det(a)%矩阵的行列式
C = null(a)%矩阵零空间
D = orth(a)%适用于矩阵范围的标准正交基

for i = 1:3

    for j = 1:3

        if mod(a(i, j), 2) == 1
            b(i, j) = 1;
        else
            b(i, j) = 0;
        end

    end

end

fprintf('b =\n');
disp(b);
