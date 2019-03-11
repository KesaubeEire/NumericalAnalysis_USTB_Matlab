
%% ��Ŀ:1��
% �ֱ��ò���������� Newton ����ⷽ�� [x-exp((3*x))+5=0] �������븺����
% % % % % % % % % % % % % % % % % % %
format long
set(gca, 'XGrid', 'on');
set(gca, 'YGrid', 'on');

%% 1.��ȡ��ȡֵ��Χ: -5 �� (0.5,0.6)
% x1 = -6:0.01:1.9;
% y2 = x1 - exp((3*x1))+5;
% plot(x1,y2);

% y1 = formula(0.5726218);
% fprintf('%f\n',y1);

%% ��ʽ����
% 2.1 ����
% % % % % % % % % % % % % % % % % % %
n0 = 20;
p0 = 0.2;
p1 = -4;
tol = 1e-4;

% 2.2 ���������:������-5������
fprintf('\n\n���������\n');
fprintf('p0:\n');
for x = 0.5
    fpi(x, tol);
end

% 2.3 newton����:
fprintf('\n\nnewton����\n');
fprintf('p0:\n');
funNewton(p0, 100, tol);
fprintf('p1:\n');
funNewton(p1, 100, tol);
% % % % % % % % % % % % % % % % % % %
% ���㷨
function xc = fpi(x0, tol)
x(1) = x0;
i = 1;
k = 0;
fprintf('���������㷽��:y = exp((3 * x)) - 5;\n');
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
    fprintf('%.2f����ʧ��\n', x0);
    return
end

xc = x(i+1);
fprintf('%.2f:\t%.4d\n', x0, x(i+1));

fprintf('\n');

x(1) = x0;
i = 1;
k = 0;
fprintf('���������㷽��:y = 3 * x - ln(5);\n');
while i < 100
    x(i + 1) = formula_g2(x(i));
    disp(['the ', num2str(i), ' time .', num2str(x(i + 1))])
    if (abs(x(i+1)-x(i)) < tol)
        k = 1;
        break
    end
    i = i + 1;
end

if (k == 0)
    fprintf('%.2f����ʧ��\n', x0);
    return
end

xc = x(i+1);
fprintf('%.2f:\t%.4d\n', x0, x(i+1));
end
% % % % % % % % % % % % % % % % % % %

% ţ�ٷ�
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

% ����:ԭʽ
function y = formula(x)
y = x - exp((3 * x)) + 5;
end
% ����:������ʽ_1_����
function y = formula_g(x)
y = exp((3 * x)) - 5;
end
% ����:������ʽ_1_����
function y = formula_g2(x)
y = (log(x+5)) / 3;
end

%% ����
% �ӽ������,���������������4��
% ţ�ٷ��ڲ�ͬ�ĳ�ʼ����������ǲ�ͬ��
% ʹ��ţ�ٷ�����Ҫע���ʼ���ѡ��