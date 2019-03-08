% 题目:3、
% 应用Newton法求f(x）的零点，
% e=10^-6这里f(x)=x-sinx。
% 再用求重根的两种方法求 f(x）的零点。
%
format long;
% 1.作图获取粗取值范围:
x1 = -pi * 10:0.1:10 * pi;
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


fprintf('第一种求重根方法:\n');
fx = x - sin(x);
dfx = diff(fx, x);
fprintf('fx =\t');
disp([fx]);
fprintf('dfx =\t');
disp([dfx]);
fprintf('u(x) = fx/dfx =\t');
ux = fx / dfx;
disp([ux]);
fprintf('du(x) =\t');
dux = diff(ux, x);
disp([dux]);
fprintf('\n');
fprintf('g(x) = x - u(x)/du(x) =\t');
fprintf('\n');
g = x - ux / dux;
disp([g]);

e = 1;
fprintf('e:\t%d\n', e);
for i = 1:4
    e = vpa(subs(g, e));
    fprintf('e:\t%d\n', e);
end
fprintf('第一种求重根方法求得为0.');
fprintf('\n\n第二种求重根方法:\n');
fprintf('牛顿法结果无限接近于0,已知根为0,需要确认是几重根\n以下为formula的导数:\n');

syms x;
fx = x - sin(x);
fxx = fx;
for i = 1:10
    
    if (i <= 3)
        fxx = diff(fxx, x);
        disp([fxx]);
    else
        break;
    end
end
fprintf('当x=0时,二阶导数为0,三阶导数不为0,可知0为3重根\nm=3\n');
fprintf('由公式: g(x) = x - m*f(x)/f`(x) \n');
m = 3;
g2 = x - m * fx / dfx;

e = 1;
fprintf('e:\t%d\n', e);
for i = 1:4
    e = vpa(subs(g2, e));
    fprintf('e:\t%d\n', e);
end
fprintf('第二种求重根方法求得为0.\n');

% % % % % % % % % % % % % % % % % % % % % % %
% a.牛顿法
function xk = funNewton(x0, max_steps, tol)
syms x
dif_f = matlabFunction(diff(formula(x)));
clear x
x = x0;
for k = 1:max_steps
    xk = x;
    %     这里注意,需要
    %     disp(['the ', num2str(k), ' time is ', num2str(x)])
    %xk to save the last time value of x
    x = x - formula(x) / dif_f(x);
    %newton solve
    if (abs(xk-x) < tol)
        %decide whether to break out
        break;
    end
end
end


% % % % % % % % % % % % % % % % % % % % % % %
% 函数:原式
function y = formula(x)
y = x - sin(x);
end