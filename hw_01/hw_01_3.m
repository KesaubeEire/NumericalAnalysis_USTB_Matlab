clc;clear;
% ��Ŀ:3��
% Ӧ��Newton����f(x������㣬
% e=10^-6����f(x)=x-sinx��
% �������ظ������ַ����� f(x������㡣
%
format long;
% 1.��ͼ��ȡ��ȡֵ��Χ:
x1 = -pi * 10:0.1:10 * pi;
y2 = formula(x1);
plot(x1, y2);
% % % % % % % % % % % % % % % % % % % % % % %

fprintf(num2str(formula(-1.752325e-08)));

p0 = -0.0002;
tol = 1e-12;

fprintf('\nţ�ٷ���ʼ:%%%%%%%%%%%%%%%%%%%%\n');
res_nt = funNewton(p0, 100, tol);
fprintf('res:\t%d\n', res_nt);
fprintf('���ԼΪ0\n');
fprintf('ţ�ٷ�����%%%%%%%%%%%%%%%%%%%%\n\n');


fprintf('��һ�����ظ�����:\n');
syms x;
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
fprintf('��һ�����ظ��������Ϊ0.');
fprintf('\n\n�ڶ������ظ�����:\n');
fprintf('ţ�ٷ�������޽ӽ���0,��֪��Ϊ0,��Ҫȷ���Ǽ��ظ�\n����Ϊformula�ĵ���:\n');

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
fprintf('��x=0ʱ,���׵���Ϊ0,���׵�����Ϊ0,��֪0Ϊ3�ظ�\nm=3\n');
fprintf('�ɹ�ʽ: g(x) = x - m*f(x)/f`(x) \n');
m = 3;
g2 = x - m * fx / dfx;

e = 1;
fprintf('e:\t%d\n', e);
for i = 1:4
    e = vpa(subs(g2, e));
    fprintf('e:\t%d\n', e);
end
fprintf('�ڶ������ظ��������Ϊ0.\n');

% % % % % % % % % % % % % % % % % % % % % % %
% a.ţ�ٷ�
function xk = funNewton(x0, max_steps, tol)
syms x
dif_f = matlabFunction(diff(formula(x)));
clear x
x = x0;
for k = 1:max_steps
    xk = x;
    %     ����ע��,��Ҫ
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
% ����:ԭʽ
function y = formula(x)
y = x - sin(x);
end