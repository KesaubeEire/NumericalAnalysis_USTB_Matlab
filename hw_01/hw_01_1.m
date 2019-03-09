%% 题目:1、
% 分别用不动点迭代与 Newton 法求解方程 [x-exp((3*x))+5=0] 的正根与负根。
% % % % % % % % % % % % % % % % % % %
format long
set(gca, 'XGrid', 'on');
set(gca, 'YGrid', 'on');

%% 1.获取粗取值范围: -5 和 (0.5,0.6)
% x1 = -6:0.01:1.9;
% y2 = x1 - exp((3*x1))+5;
% plot(x1,y2);

% y1 = formula(0.5726218);
% fprintf('%f\n',y1);

%% 正式运算
% 2.1 参数
% % % % % % % % % % % % % % % % % % %
n0 = 20;
p0 = 0.2;
p1 = -4;
tol = 1e-4;

% 2.2 不动点迭代:总是在-5处收敛
fprintf('\n\n不动点迭代\n');
fprintf('p0:\n');
for x = 0.5
    fpi(x, tol);
    
end

% 2.3 newton迭代:
fprintf('\n\nnewton迭代\n');
fprintf('p0:\n');
funNewton(p0, 100, tol);
fprintf('p1:\n');
funNewton(p1, 100, tol);
% % % % % % % % % % % % % % % % % % %
% 定点法
function xc = fpi(x0, tol)
x(1) = x0;
i = 1;
k = 0;
while i < 10
    x(i + 1) = formula_g(x(i));
    disp(['the ', num2str(i), ' time .'])
    if (abs(x(i+1)-x(i)) < tol)
        k = 1;
        break
    end
    i = i + 1;
end

if (k == 0)
    fprintf('%.2f迭代失败\n', x0);
    return
end

xc = x(i+1);
fprintf('%.2f:\t%.4d\n', x0, x(i+1));
end
% % % % % % % % % % % % % % % % % % %

% 牛顿法
function xk = funNewton(x0, max_steps, tol)
syms x
dif_f = matlabFunction(diff(formula(x)));
clear x
x = x0;
for k = 1:max_steps
    xk = x;
    disp(['the ', num2str(k), ' time is ', num2str(x)])
    %xk to save the last time value of x
    x = x - formula(x) / dif_f(x);
    %newton solve
    if (abs(xk-x) < tol)
        %decide whether to break out
        break;
    end
end
end

% 函数:原式
function y = formula(x)
y = x - exp((3 * x)) + 5;
end
% 函数:不动点式
function y = formula_g(x)
y = exp((3 * x)) - 5;
end

%% 分析
% 从结果来看,不动点迭代经过了4次
% 牛顿法在不同的初始点迭代次数是不同的
% 使用牛顿法还是要注意初始点的选择