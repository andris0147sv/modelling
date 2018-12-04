%% ���������� FDTD. ������ 1.3
% �������������� ������� ������� ���������������� � ���� ������� (TFSF boundary)
clear

% �������� ������������� ���������� ������������
W0 = 120 * pi;

% ����� ������� � ��������
maxTime = 250;

% ������ ������� ������������� � ��������
maxSize = 200;

% ��������� �������, ��������������� ����
probePos = 60;

% ��������� ��������� �����������
sourcePos = 50;

Nwg = 20;
Ndg = 50;
N_lambda = 8;


Ez = zeros (1, maxSize);
Hy = zeros (1, maxSize);

% ����, ������������������ � ������� � ����������� �� �������
probeTimeEz = zeros (1, maxTime);

figure

for t = 1: maxTime
    % ������ ���������� ���� H
    Hy(maxSize) = Hy(maxSize - 1);
    for m = 1: maxSize - 1
        % �� ���� ������ Hy(n) ������ �������� ���������� Hy
        % �� ���������� ������ �������
        Hy(m) = Hy(m) + (Ez(m + 1) - Ez(m)) / W0;
    end
    
    Hy(sourcePos - 1) = Hy(sourcePos - 1) -...
        sin (2 * pi * t / N_lambda) *...
        exp (- ((t - Ndg) / Nwg) ^ 2) / W0;
    
    % ������ ���������� ���� E
    Ez(1) = Ez(2);
    for m = 2: maxSize
        % �� ���� ������ Ez(n) ������ �������� ���������� EzS
        % �� ���������� ������ �������
        Ez(m) = Ez(m) + (Hy(m) - Hy(m - 1)) * W0;
    end

    % �������� �����������
    Ez(sourcePos) = Ez(sourcePos) +...
        sin (2 * pi * (t + 0.5 - (-0.5)) / N_lambda) *...
        exp (- ((t + 0.5 - (-0.5) - Ndg) / Nwg) ^ 2);
    
    % ����������� ���� � �����
    probeTimeEz(t) = Ez(probePos);
    
    plot (Ez);
    xlim ([1, maxSize]);
    ylim ([-1.1, 1.1]);
    xlabel ('x, ������')
    ylabel ('Ez, �/�')
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
ylim ([-1.1, 1.1]);
grid on