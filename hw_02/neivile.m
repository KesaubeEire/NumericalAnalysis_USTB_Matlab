clc
clear

%% Neville's
clc;
clear;
n = 12; % 在这里更改
k = 100;

% function
x0 = linspace(-1, 1, k);
re = [];
x = linspace(-1, 1, n);
for h = 1:k
    q = zeros(n);
    for i = 1:n
        q(i, 1) = 1 / (1 + 25 * x(i)^2);
    end
    for i = 2:n
        for j = 2:i
            q(i, j) = ((x0(h) - x(i-j+1)) * q(i, j-1) - (x0(h) - x(i)) * q(i-1, j-1)) / (x(i) - x(i-j+1));
        end
    end
    re(h) = q(n, n);
end
y = 1 ./ (1 + 25 * x0.^2);
y1 = 1 ./ (1 + 25 * x.^2);
plot(x0, y, 'r');
hold on;
plot(x0, re, 'k', x, y1, 'k+');