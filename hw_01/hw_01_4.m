clc;
clear;

%% 题目:4、
% 应用Newton法求f(x）的零点，e=10^-6
% 这里f(x)=x-sinx。
% 再用 Steffensen's method 加速其收敛。

%% 0.作图获取粗取值范围:
x1 = -pi * 10:0.1:10 * pi;
y2 = formula(x1);
plot(x1, y2);
% % % % % % % % % % % % % % % % % % % % % % %

% fprintf(num2str(formula(-1.752325e-08)));

%% 1.牛顿法
p0 = -0.0002;
tol = 1e-12;

fprintf('\n牛顿法开始:%%%%%%%%%%%%%%%%%%%%\n');
res_nt = funNewton(p0, 100, tol);
fprintf('res:\t%d\n', res_nt);
fprintf('结果约为0\n');
fprintf('牛顿法结束%%%%%%%%%%%%%%%%%%%%\n\n');


%求y=x-sinx
%y'=1-cosx
%% 2.用Steffensen’s method使其加速收敛
p0 = 1;
p1 = p0 - (p0 - sin(p0)) / (1 - cos(p0));
N = 50;
tol = 1e-6;
n = 0;
p(1) = p0;
while n <= N
    fprintf('循环中:%d\n', n);
    p(2) = p(1) - (p(1) - sin(p(1))) / (1 - cos(p(1)));
    p(3) = p(2) - (p(2) - sin(p(2))) / (1 - cos(p(2)));
    p1 = p(1) - (p(2) - p(1))^2 / (p(3) - 2 * p(2) + p(1));
    f0 = p1 - sin(p1);
    if abs(f0) < tol
        break
    end
    n = n + 1;
    p(1) = p1;
end
fprintf('最后结果:%d\n', p1);
fprintf('循环次数:%d\n', n+1);

% % % % % % % % % % % % % % % % % % % % % % %
% a.牛顿法
function xk = funNewton(x0, max_steps, tol)
syms x
dif_f = matlabFunction(diff(formula(x)));
clear x
x = x0;
for k = 1:max_steps
    xk = x;
    %     disp(['the ', num2str(k), ' time is ', num2str(x)])
    x = x - formula(x) / dif_f(x);
    if (abs(xk-x) < tol)
        break;
    end
end
end


% % % % % % % % % % % % % % % % % % % % % % %
% 函数:原式
function y = formula(x)
y = x - sin(x);
end


%% 3.分析
