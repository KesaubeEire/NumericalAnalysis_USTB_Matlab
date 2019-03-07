% 1.获取粗取值范围:
x1 = -pi:0.1:pi;
y2 = formula(x1);
plot(x1, y2);
% % % % % % % % % % % % % % % % % % % % % % %

% % % % % % % % % % % % % % % % % % % % % % %
% 函数:原式
function y = formula(x)
y = x - sin(x);
end