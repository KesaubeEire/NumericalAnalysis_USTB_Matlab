%% ������ҵ������ţ�ٷ���ⷽ�� z^4-1=0 ��������
%
% % % % % % % % % % % % % % % % %
% �ⷽ��
%
% z*-1=0
%
% �ڸ�ƽ������ 4 ������
%
% @.2=+1,, 03,4=�� i
%
% ��ÿ--���� w������Ҫȷ������ Dw��ʹ�ö�����
%
% Zo��D%,
%
% ���� limzn = 0. ���{z, } Ϊţ�ٷ�������
%
% n-����0
%
% z4-1
%
% Zn+1 = zn һ
%
% 4 zin
%
% 'n
% % % % % % % % % % % % % % % % %

%% set para
d = 6;
tol = 1e-5;
maxIter = 100;
r = -2:0.01:2; %ʵ���鲿�ķ�Χ
[x, y] = meshgrid(r); %����ʵ���鲿��ά����
Z = x + 1i * y; %Z��Ӧ�������ƽ��

%% Define fuction
f = @(x, d) (x.^d) - 1;
fprime = @(x, d) d * (x.^(d - 1));

%% Perform Newton iterations
for k = 1:maxIter
    Z = Z - (f(Z, d) ./ fprime(Z, d));
end

%% Find d roots of unity, and the  mask
renderMat = 0;
for j = 1:d
    root = exp(2*pi*1i/d)^j; % the jth root
    Mj = abs(Z-root); %  distance  Z��ÿ�㶼������ľ���
    % Each root gets a unique number in [1,d]
    mask = (Mj <= tol) * j; %Mj<=tol�������������߼�����
    %��������Ϊ1*j�������㲿��Ϊ0
    renderMat = renderMat + mask;
    %������֮��renderMat�������ڵ�j�������������ݶ���j
    %��ô�����ڵ�j������������ͬһ����ɫ
end
colormap(hsv(d+1)); % Set the color map
imagesc(r, r, renderMat) % Render the fractal
xlabel('Re(Z)');
ylabel('Im(Z)');
h = colorbar;
set(h, 'ytick', (2 * (0:d) + 1)*d/(d + 1)/2);
str = arrayfun(@(x)num2str(x, '%.2f'), exp(2*pi*1i/d).^(1:d), 'uniformoutput', false);
set(h, 'yticklabel', [{'δ����'}, str]);