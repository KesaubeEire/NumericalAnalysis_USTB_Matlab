% 题目:4、
% 应用Newton法求f(x）的零点，e=10^-6
% 这里f(x)=x-sinx。
% 再用 Steffensen's method 加速其收敛。

% 1.作图获取粗取值范围:
x1 = -pi*10:0.1:10*pi;
y2 = formula(x1);
plot(x1, y2);
% % % % % % % % % % % % % % % % % % % % % % %

% fprintf(num2str(formula(-1.752325e-08)));

p0 = -0.0002;
tol = 1e-12;

fprintf('\n牛顿法开始:%%%%%%%%%%%%%%%%%%%%\n');
res_nt = funNewton(p0, 100, tol);
fprintf('res:\t%d\n', res_nt);
fprintf('结果约为0\n');
fprintf('牛顿法结束%%%%%%%%%%%%%%%%%%%%\n\n');


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

