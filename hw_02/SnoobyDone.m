clc
clear

%% 曲线坐标
% 坐标
x1 = [1, 2, 5, 6, 7, 8, 10, 13, 17];
y1 = [3.0, 3.7, 3.9, 4.2, 5.7, 6.6, 7.1, 6.7, 4.5];
g1_0 = 1.0;
g1_1 = -0.67;

x2 = [17, 20, 23, 24, 25, 27, 27.7];
y2 = [4.5, 7.0, 6.1, 5.6, 5.8, 5.2, 4.1];
g2_0 = 3.0;
g2_1 = -4.0;

x3 = [27.7, 28, 29, 30];
y3 = [4.1, 4.3, 4.1, 3.0];
g3_0 = 0.33;
g3_1 = -1.5;

[outx1, outy1] = cubicSpline(x1, y1, g1_0, g1_1);
[outx1, outy1] = cubicSpline(x2, y2, g2_0, g2_1);
[outx1, outy1] = cubicSpline(x3, y3, g3_0, g3_1);
% plot(x1, y1, 'b', x2, y2, 'r', x3, y3, 'g');

%% 原始函数
function [outx, outy] = cubicSpline(x, y, fp0, fpn)
n = length(x);
y1_deri = fp0;
yn_deri = fpn;

for i = 1:n - 1
    h(i) = x(i+1) - x(i);
end
% fprintf('计算 h 结果为:\n');
% h
for i = 2:n - 1
    u(i-1) = h(i-1) / (h(i-1) + h(i));
    lamda(i) = h(i) / (h(i-1) + h(i));
end

u(n-1) = 1;
lamda(1) = 1;
% fprintf('计算 μ 结果为: \n');
% u
% fprintf('计算 λ 结果为：\n');
% lamda

for i = 2:n - 1
    d(i) = 6 * ((y(i+1) - y(i)) / (x(i+1) - x(i)) - (y(i) - y(i-1)) / (x(i) - x(i-1))) / (h(i-1) + h(i));
end

d(1) = 6 / h(1) * ((y(2) - y(1)) / (x(2) - x(1)) - y1_deri);
d(n) = 6 / h(n-1) * (yn_deri - ((y(n) - y(n-1)) / (x(n) - x(n-1))));
% fprintf('计算 d 的结果：\n');
% d

matrix1 = zeros(n, n);
for i = 1:n - 1
    matrix1(i, i) = 2;
    matrix1(i, i+1) = lamda(i);
    matrix1(i+1, i) = u(i);
end
matrix1(n, n) = 2;
matrix1;
% fprintf('求得 M 结果:\n');
M = matrix1 \ d';
for i = 1:n - 1
    clear S
    syms t
    k = x(i):0.001:x(i+1);
    
    S = M(i) * (x(i+1) - t)^3 / (6 * h(i)) + M(i+1) * (t - x(i))^3 / (6 * h(i)) + (y(i) - M(i) * h(i)^2 / 6) * (x(i+1) - t) / h(i) + (y(i+1) - M(i+1) * h(i)^2 / 6) * (t - x(i)) / h(i);
    s = M(i) * (x(i+1) - k).^3 / (6 * h(i)) + M(i+1) * (k - x(i)).^3 / (6 * h(i)) + (y(i) - M(i) * h(i)^2 / 6) * (x(i+1) - k) / h(i) + (y(i+1) - M(i+1) * h(i)^2 / 6) * (k - x(i)) / h(i);
    
    hold on;
    outx = k;
    outy = s;
    fprintf('绘图区间为[ %.3f : %.3f]\n', x(i), x(i+1));
    plot(k, s,'*');
    
end
end