% 1.��ȡ��ȡֵ��Χ:
x1 = -pi:0.1:pi;
y2 = formula(x1);
plot(x1, y2);
% % % % % % % % % % % % % % % % % % % % % % %

% % % % % % % % % % % % % % % % % % % % % % %
% ����:ԭʽ
function y = formula(x)
y = x - sin(x);
end