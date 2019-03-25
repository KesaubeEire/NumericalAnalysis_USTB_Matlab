clc
clear

%% ����Romberg����
clear;

% fun = @(x)sin(x+eps) / (x + eps)
fun = @(x)sqrt(1+cos(x).^2)
a = 0
b = 48
epsilon = 0.000001
Romberg(fun, a, b, epsilon)

%% Main Function
%Romberg.m
function [int_value] = Romberg(fun, a, b, epsilon)
%Romberg: ����������ʽ������������Ķ�����ֵ�����ػ���ֵ�ļ�����������ӡ������̾���
% inputs:
%   -fun�����ֺ������
%   -a/b������������
%   -epsilon: �����ľ���
% Outputs:
%   -int_value: ���ֽ��
% Example:
%   fun=@(x)sin(x + eps) / (x + eps)
%   value = Romberg(fun, 0, 1, 0.000001)

T = [];
delta = intmax;
n = 1;
while (delta > epsilon)
    T(n, 1) = Tn(fun, a, b, n-1);
    sig = 2;
    if (n >= sig)
        T(n, sig) = T(n, sig-1) + 1 / (4.^(sig - 1) - 1) * (T(n, sig-1) - T(n-1, sig-1));
        sig = sig + 1;
        if (n >= sig)
            T(n, sig) = T(n, sig-1) + 1 / (4.^(sig - 1) - 1) * (T(n, sig-1) - T(n-1, sig-1));
            sig = sig + 1;
            if (n >= sig)
                T(n, sig) = T(n, sig-1) + 1 / (4.^(sig - 1) - 1) * (T(n, sig-1) - T(n-1, sig-1));
                sig = sig + 1;
                if (n >= sig)
                    T(n, sig) = T(n, sig-1) + 1 / (4.^(sig - 1) - 1) * (T(n, sig-1) - T(n-1, sig-1));
                    sig = sig + 1;
                    if (n >= sig)
                        T(n, sig) = T(n, sig-1) + 1 / (4.^(sig - 1) - 1) * (T(n, sig-1) - T(n-1, sig-1));
                        sig = sig + 1;
                        if (n >= sig)
                            T(n, sig) = T(n, sig-1) + 1 / (4.^(sig - 1) - 1) * (T(n, sig-1) - T(n-1, sig-1));
                            sig = sig + 1;
                            if (n >= sig)
                                T(n, sig) = T(n, sig-1) + 1 / (4.^(sig - 1) - 1) * (T(n, sig-1) - T(n-1, sig-1));
                                sig = sig + 1;
                                if (n >= sig)
                                    T(n, sig) = T(n, sig-1) + 1 / (4.^(sig - 1) - 1) * (T(n, sig-1) - T(n-1, sig-1));
                                    sig = sig + 1;
                                    if (n >= sig)
                                        T(n, sig) = T(n, sig-1) + 1 / (4.^(sig - 1) - 1) * (T(n, sig-1) - T(n-1, sig-1));
                                        sig = sig + 1;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    if (n >= 9)
        delta = abs(T(n, 4)-T(n-1, 4));
    end
    n = n + 1;
end
int_value = T(n-1, 4);
T
end

%% Romberg Function
%Tn.m
function [T1] = Tn(fun, a, b, k)
%Tn: ���ø������ι�ʽ���㸴���������Ϊ(2^k)ʱ�Ļ���ֵ
% inputs:
%   -fun�����ֺ������
%   -a/b������������
%   -k�������������Ϊ2^k��
% Outputs:
%   -T1; �����������Ϊ(2^k)ʱ�Ļ���ֵ
% Example
%   fun=@(x)4./(1+x^2);
%   T = Tn(fun,0,1,0)

T1 = (b - a) / 2 * (fun(a) + fun(b));
n = 1;
cnt = 1;
del = intmax;
while (n < 2^k)
    sum = 0;
    for i = 1:n
        sum = sum + fun(a+(2 * i - 1)*(b - a)/(2 * n));
    end
    T2 = 1 / 2 * (T1 + (b - a) / n * sum);
    del = T2 - T1;
    T1 = T2;
    n = n * 2;
    cnt = cnt + 1;
end
T1;

end
