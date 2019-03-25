% 史努比狗狗
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

% 运行
[outx1, outy1] = cubicSpline(x1, y1, g1_0, g1_1, 8);
[outx2, outy2] = cubicSpline(x2, y2, g2_0, g2_1, 6);
[outx3, outy3] = cubicSpline(x3, y3, g3_0, g3_1, 3);

% 画图
plot(x1, y1, 'b', x2, y2, 'r', x3, y3, 'g', outx1, outy1, 'p', outx2, outy2, 'p', outx3, outy3, 'p');

%% 原始函数
function [outx, outy] = cubicSpline(x, a, g0, g1, n)
% 0
h = zeros(n+1);
l = zeros(n+1);
u = zeros(n+1);
z = zeros(n+1);
c = zeros(n+1);
alpha = zeros(n+1);
b = zeros(n+1);
d = zeros(n+1);
outx = [];
outy = [];

% 1
for i = 1:n
    h(i) = x(i+1) - x(i);
end

% 2
alpha_0 = 3 * (a(2) - a(1)) / h(1) - 3 * g0;
alpha_n = 3 * g1 - 3 * (a(n) - a(n-1)) / h(n-1);

% 3
alpha(1) = alpha_0;
for i = 2:n
    alpha(i) = 3 / h(i) * (a(i+1) - a(i)) - 3 / h(i-1) * (a(i) - a(i-1));
end
alpha(n+1) = alpha_n;

% for i = 1:n+1
%     alpha(i)
% end

% 4
l(1) = 2 * h(1);
u(1) = 0.5;
z(1) = alpha(1) / l(1);

% 5
for i = 2:n
    l(i) = 2 * (x(i+1) - x(i-1)) - h(i-1) * u(i-1);
    u(i) = h(i) / l(i);
    z(i) = (alpha(i) - h(i-1) * z(i-1)) / l(i);
end

% 6
l(n+1) = h(n) * (2 - u(n));
z(n+1) = (alpha(n+1) - h(n) * z(n)) / l(n+1);
c(n+1) = z(n+1);


% 7
for j = n:-1:1
    c(j) = z(j) - u(j) * c(j+1);
    b(j) = (alpha(j+1) - alpha(j)) / h(j) - h(j) * (c(j+1) + 2 * c(j)) / 3;
    d(j) = (c(j+1) - c(j)) / (3 * h(j));
end

% 8 : plot
allx = [];
ally = [];
for j = 1:n
    xx = x(j):0.01:x(j+1);
    ak = a(j) + b(j) * (xx - x(j)) + c(j) * (xx - x(j)).^2 + d(j) * (xx - x(j)).^3;
    allx = [allx, xx];
    ally = [ally, ak];
    
end

outx = [outx, allx];
outy = [outy, ally];
end