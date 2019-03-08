% ��Ŀ:4��
% Ӧ��Newton����f(x������㣬e=10^-6
% ����f(x)=x-sinx��
% ���� Steffensen's method ������������

% 1.��ͼ��ȡ��ȡֵ��Χ:
x1 = -pi*10:0.1:10*pi;
y2 = formula(x1);
plot(x1, y2);
% % % % % % % % % % % % % % % % % % % % % % %

% fprintf(num2str(formula(-1.752325e-08)));

p0 = -0.0002;
tol = 1e-12;

fprintf('\nţ�ٷ���ʼ:%%%%%%%%%%%%%%%%%%%%\n');
res_nt = funNewton(p0, 100, tol);
fprintf('res:\t%d\n', res_nt);
fprintf('���ԼΪ0\n');
fprintf('ţ�ٷ�����%%%%%%%%%%%%%%%%%%%%\n\n');


% % % % % % % % % % % % % % % % % % % % % % %
% a.ţ�ٷ�
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
% ����:ԭʽ
function y = formula(x)
y = x - sin(x);
end

