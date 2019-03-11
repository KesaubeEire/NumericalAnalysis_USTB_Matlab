clc;
clear;

%% ��Ŀ:4��
% Ӧ��Newton����f(x������㣬e=10^-6
% ����f(x)=x-sinx��
% ���� Steffensen's method ������������

%% 0.��ͼ��ȡ��ȡֵ��Χ:
x1 = -pi * 10:0.1:10 * pi;
y2 = formula(x1);
plot(x1, y2);
% % % % % % % % % % % % % % % % % % % % % % %

% fprintf(num2str(formula(-1.752325e-08)));

%% 1.ţ�ٷ�
p0 = -0.0002;
tol = 1e-12;

fprintf('\nţ�ٷ���ʼ:%%%%%%%%%%%%%%%%%%%%\n');
res_nt = funNewton(p0, 100, tol);
fprintf('res:\t%d\n', res_nt);
fprintf('���ԼΪ0\n');
fprintf('ţ�ٷ�����%%%%%%%%%%%%%%%%%%%%\n\n');


%��y=x-sinx
%y'=1-cosx
%% 2.��Steffensen��s methodʹ���������
p0 = 1;
p1 = p0 - (p0 - sin(p0)) / (1 - cos(p0));
N = 50;
tol = 1e-6;
n = 0;
p(1) = p0;
while n <= N
    fprintf('ѭ����:%d\n', n);
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
fprintf('�����:%d\n', p1);
fprintf('ѭ������:%d\n', n+1);

% % % % % % % % % % % % % % % % % % % % % % %
% a.ţ�ٷ�
function xk = funNewton(x0, max_steps, tol)
syms x
dif_f = matlabFunction(diff(formula(x)));
clear x
x = x0;
for k = 1:max_steps
    xk = x;
        disp(['the ', num2str(k), ' time is ', num2str(x)])
    x = x - formula(x) / dif_f(x);
    if (abs(xk-x) < tol)
        break;
    end
end
end


% % % % % % % % % % % % % % % % % % % % % % %
% ����:ԭʽ
function y = formula(x)
y = x - sin(x);
end


%% 3.����
% ����ѭ�����ﵽ����
% ������10�౶
