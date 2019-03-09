%% 补充作业：画出牛顿法求解方程 z^4-1=0 的收敛域。
%
% % % % % % % % % % % % % % % % %
% 解方程
%
% z*-1=0
%
% 在复平面上有 4 个根：
%
% @.2=+1,, 03,4=土 i
%
% 对每--个根 w，我们要确定集合 Dw，使得对任意
%
% Zo∈D%,
%
% 都有 limzn = 0. 这里，{z, } 为牛顿法迭代：
%
% n-》∞0
%
% z4-1
%
% Zn+1 = zn 一
%
% 4 zin
%
% 'n
% % % % % % % % % % % % % % % % %

%% set para
d = 6;
tol = 1e-5;
maxIter = 100;
r = -2:0.01:2; %实部虚部的范围
[x, y] = meshgrid(r); %产生实部虚部二维网格
Z = x + 1i * y; %Z对应网格的虚平面

%% Define fuction
f = @(x, d) (x.^d) - 1;
fprime = @(x, d) d * (x.^(d - 1));

%% Perform Newton iterations
for k = 1:maxIter
    Z = Z - (f(Z, d) ./ fprime(Z, d));
end

%% Find d roots of unity, and the  mask
renderMat = 0;
for j = 1:d
    root = exp(2*pi*1i/d)^j; % the jth root
    Mj = abs(Z-root); %  distance  Z中每点都这个根的距离
    % Each root gets a unique number in [1,d]
    mask = (Mj <= tol) * j; %Mj<=tol返回满足误差的逻辑矩阵
    %满足误差部分为1*j，不满足部分为0
    renderMat = renderMat + mask;
    %加起来之后renderMat中收敛于第j个根的区域数据都是j
    %那么收敛于第j个根的区域都是同一种颜色
end
colormap(hsv(d+1)); % Set the color map
imagesc(r, r, renderMat) % Render the fractal
xlabel('Re(Z)');
ylabel('Im(Z)');
h = colorbar;
set(h, 'ytick', (2 * (0:d) + 1)*d/(d + 1)/2);
str = arrayfun(@(x)num2str(x, '%.2f'), exp(2*pi*1i/d).^(1:d), 'uniformoutput', false);
set(h, 'yticklabel', [{'未收敛'}, str]);