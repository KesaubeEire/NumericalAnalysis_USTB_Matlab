clc;
clear;

%% ��Ŀ:2��
% Use each of the following methods to find a solution in [0.1,1]
% accurate to within 1e-4 for [600 * x.^4 - 550 * x.^3 + 200 * x.^2 - 20 * x - 1]
%
% a.bisection
% b.newton's
% c.secant
% d.false position
% e.muller
%

%% 1.��ȡ��ȡֵ��Χ:
x1 = 0.1:0.001:0.3;
y2 = formula(x1);
y3 = 0;
plot(x1, y2,x1,y3,'*');
%%%%%%%%%%%%%%%%%%%%%

%% 2.��ʽ����:
p0 = 0.1;
p1 = 1;
tol = 1e-7;

a = formula(0.232354);
disp(['��׼��: ', num2str(a)]);

index = 1;
fprintf('%d:', index);
index = index + 1;
fprintf('\n���ַ���ʼ:%%%%%%%%%%%%%%%%%%%%\n');
res_bi = bisection(p0, p1, tol);
fprintf('res:\t%d\n', res_bi);
fprintf('���ַ�����%%%%%%%%%%%%%%%%%%%%\n\n');

fprintf('%d:', index);
index = index + 1;
fprintf('\nţ�ٷ���ʼ:%%%%%%%%%%%%%%%%%%%%\n');
res_nt = funNewton(p0, 100, tol);
fprintf('res:\t%d\n', res_nt);
fprintf('ţ�ٷ�����%%%%%%%%%%%%%%%%%%%%\n\n');

fprintf('%d:', index);
index = index + 1;
fprintf('\n���߷���ʼ:%%%%%%%%%%%%%%%%%%%%\n');
res_se = mysecant(p0, p1, 100, tol);
fprintf('res:\t%d\n', res_se);
fprintf('���߷�����%%%%%%%%%%%%%%%%%%%%\n\n');

fprintf('%d:', index);
index = index + 1;
fprintf('\n��λ����ʼ:%%%%%%%%%%%%%%%%%%%%\n');
res_fp = falseposition(p0, p1, 1000, tol);
fprintf('res:\t%d\n', res_fp);
fprintf('��λ������%%%%%%%%%%%%%%%%%%%%\n\n');


p2 = 0.3;

fprintf('%d:', index);
index = index + 1;
fprintf('\nMuller����ʼ:%%%%%%%%%%%%%%%%%%%%\n');
res_mu = muller_method(p0, p1, p2, 1000, tol);
fprintf('res:\t%d\n', res_mu);
fprintf('Muller������%%%%%%%%%%%%%%%%%%%%\n\n');

%% 3.��������
% a.���ַ�
function res = bisection(p0, p1, tol)
% assign
p = (p0 + p1) / 2;
k = 0;
% ����
while (abs(p0-p1) > tol)
    k = k + 1;
    p = (p0 + p1) / 2;
    disp(['the ', num2str(k), ' time is ', num2str(p)])
    if (formula(p0) * formula(p) <= 0)
        if (p == 0)
            res = p;
            return
        end
        p1 = p;
    else
        p0 = p;
    end
end
res = p;
end

% b.ţ�ٷ�
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

% c.secant method
function x = mysecant(x0, x1, max_steps, tol)
% Solves f(x) = 0 by doing n steps of the secant method
% starting with x0 and x1.
% x0 -- starting guess , a number
% x1 -- second starting geuss
% n -- the number of steps to do
% Output : x -- the approximate solution
y0 = formula(x0);
y1 = formula(x1);
for i = 1:max_steps % Do n times
    x = x1 - (x1 - x0) * y1 / (y1 - y0); % secant formula .
    y = formula(x); % y value at the new approximate solution .
    % Move numbers to get ready for the next step
    x0 = x1;
    y0 = y1;
    x1 = x;
    y1 = y;
    disp(['the ', num2str(i), ' time is ', num2str(x)])
    if (abs(x0-x1) < tol)
        break;
    end
end
end

% d.false position method_1
function res = falseposition(x0, x1, max_steps, tol)
f = @(x) 600 * x^4 - 550 * x^3 + 200 * x^2 - 20 * x - 1;
for i = 0:max_steps
    x2 = x1 - (f(x1) * (x1 - x0) / (f(x1) - f(x0)));
    %     disp(['the ', num2str(i), ' time is ', num2str(x2)]);
    c = f(x2);
    absolute_c = abs(c);
    
    if absolute_c < tol
        res = x2;
        disp(['the ', num2str(i), ' time is ', num2str(x2)]);
        break
    end
    
    if f(x0) * c < 0
        x1 = x2;
        continue
    else
        x0 = x2;
        continue
    end
    res = x2;
end
end

% d.false position method_2_���Լ�д��,Ч�ʸ���,��
% function res = falseposition(x0, x1, max_steps, tol)
% y0 = formula(x0);
% y1 = formula(x1);
% for i = 1:max_steps % Do n times
%     x = x1 - (x1 - x0) * y1 / (y1 - y0); % secant formula .
%     y = formula(x); % y value at the new approximate solution .
%     % Move numbers to get ready for the next step
%
%     if (y * y0 <= 0)
%         x1 = x;
%         y1 = y;
%     else
%         x0 = x;
%         y0 = y;
%     end
%
%     disp(['the ', num2str(i), ' time is ', num2str(x)])
%
%     if (abs(formula(x)) < tol)
%         res = x;
%         break;
%     end
%     res = x;
% end
% end

% d.muller method
function res = muller_method(p0, p1, p2, max_steps, TOL)
f = @ (x) 600 * x^4 - 550 * x^3 + 200 * x^2 - 20 * x - 1;
h1 = p1 - p0;
h2 = p2 - p1;
DELTA1 = (f(p1) - f(p0)) / h1;
DELTA2 = (f(p2) - f(p1)) / h2;
d = (DELTA2 - DELTA1) / (h2 + h1);
i = 3;

while i <= max_steps
    
    b = DELTA2 + h2 * d;
    D = (b^2 - 4 * f(p2) * d)^(1 / 2);
    if abs(b-D) < abs(b+D)
        E = b + D;
    else
        E = b - D;
    end
    
    h = -2 * f(p2) / E;
    p = p2 + h;
    
    
    if abs(h) < TOL
        %         p;
        fprintf('Bingo:------------\n');
        disp(['the ', num2str(i-2), ' time is ', num2str(p)])
        fprintf('Bingo:------------\n\n');
        break
    end
    
    disp(['the ', num2str(i-2), ' time is ', num2str(p)])
    
    p0 = p1;
    p1 = p2;
    p2 = p;
    h1 = p1 - p0;
    h2 = p2 - p1;
    DELTA1 = (f(p1) - f(p0)) / h1;
    DELTA2 = (f(p2) - f(p1)) / h2;
    d = (DELTA2 - DELTA1) / (h2 + h1);
    i = i + 1;
    
    res = p;
end

if i > max_steps
    formatSpec = string('The method failed after N0 iterations,N0= %d \n');
    fprintf(formatSpec, max_steps);
end

end

% % % % % % % % % % % % % % % % % % % % % % %
% ����:ԭʽ
function y = formula(x)
y = 600 * x.^4 - 550 * x.^3 + 200 * x.^2 - 20 * x - 1;
end

%% ����:
% �Ӽ��ֽⷨ��ѭ����������
% Muller����ţ�ٷ�ѭ����������Ϊ6��
% ���߷���֮8�� ���ַ�24��
% ��λ��ѭ���˼��ٴ�
% �õ�һ���ȽϷ������۵���ֵ