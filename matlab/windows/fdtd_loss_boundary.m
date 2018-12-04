%% ���������� FDTD. ������ 1.6
% ������� ������� ������ - ����������.
% �� ������ ������� ������� ���������� ���������� � ��������.
clear

% �������� ������������� ���������� ������������
W0 = 120 * pi;

% ������ � �����. loss = sigma * dt / (2 * eps * eps0)
loss = 0.02;

% ������ ���������������� ����
layer_x = 100;

% ��� ���������� ����������� ����������
layer_loss_x = 160;

% ����� ������� � ��������
maxTime = 750;

% ������ ������� ������������� � ��������
maxSize = 200;

% ��������� ���������
sourcePos = 50;

% ��������� �������, ��������������� ����
probePos = 60;

Ez = zeros (1, maxSize);
Hy = zeros (1, maxSize - 1);

eps = ones (size (Ez));
eps(layer_x: end) = 9.0;

% ������������ ��� ������� ���� E
ceze = ones (1, maxSize);
ceze(layer_loss_x: end) = (1 - loss) / (1 + loss);

cezh = (ones (1, maxSize) * W0 ./ eps);
cezh(layer_loss_x: end) = cezh(layer_loss_x: end) / (1 + loss);

% ������������ ��� ������� ���� H
chyh = ones (1, maxSize - 1);
chyh(layer_loss_x: end) = (1 - loss) / (1 + loss);

chye = ones (1, maxSize - 1) / W0;
chye(layer_loss_x: end) = chye(layer_loss_x: end) / (1 + loss);

% ����, ������������������ � ������� � ����������� �� �������
probeTimeEz = zeros (1, maxTime);

figure

for t = 1: maxTime
    % ������ ���������� ���� H
    for m = 1: maxSize - 1
        % �� ���� ������ Hy(n) ������ �������� ���������� Hy
        % �� ���������� ������ �������
        Hy(m) = chyh(m) * Hy(m) +...
            chye(m) * (Ez(m + 1) - Ez(m));
    end
    
    Hy(sourcePos - 1) = Hy(sourcePos - 1) -...
        exp (-(t - 30.0) ^ 2 / 100.0) / W0;
    
    % ������ ���������� ���� E
    Ez(1) = Ez(2);
    
    for m = 2: maxSize - 1
        % �� ���� ������ Ez(n) ������ �������� ���������� EzS
        % �� ���������� ������ �������
        % ������ W0 / eps ������������ cezh
        Ez(m) = ceze(m) * Ez(m) + ...
            cezh(m) * (Hy(m) - Hy(m - 1));
    end

    % �������� �����������
    Ez(sourcePos) = Ez(sourcePos) +...
        exp (-(t + 0.5 - (-0.5) - 30.0) ^ 2 / 100.0);
    
    % ����������� ���� � �����
    probeTimeEz(t) = Ez(probePos);
    
    plot (Ez);
    xlim ([1, maxSize]);
    ylim ([-1.1, 1.1]);
    xlabel ('x, ������')
    ylabel ('Ez, �/�')
    line ([layer_x, layer_x], [-1.1, 1.1], ...
        'Color',[0.0, 0.0, 0.0]);
    line ([layer_loss_x, layer_loss_x], [-1.1, 1.1], ...
        'Color',[0.0, 0.0, 0.0]);
    grid on
    hold on
    plot (probePos, 0, 'xk');
    plot (sourcePos, 0, '*r');
    hold off
    pause (0.03)
end

figure
plot (probeTimeEz)
xlabel ('t, ������')
ylabel ('Ez, �/�')
grid on