%% ���������� FDTD. ������ 1.3.1
% ������� ������� ���������������� � ���� ������� (TFSF boundary)
% �������� ��������� � �����������
clear

% �������� ������������� ���������� ������������
W0 = 120 * pi;

% ����� �������
Sc = 1;

% ����� ������� � ��������
maxTime = 350;

% ������ ������� ������������� � ��������
maxSize = 200;

% ��������� �������, ��������������� ����
probePos = 60;

% ��������� ��������� �����������
sourcePos = 50;

Ez = zeros (1, maxSize);
Hy = zeros (1, maxSize);

eps = ones (size (Ez));
eps(1: end) = 9.0;

mu = ones (size (Ez));
%mu(1: end) = 4.0;

% ����, ������������������ � ������� � ����������� �� �������
probeTimeEz = zeros (1, maxTime);

figure

for t = 1: maxTime
    % ������ ���������� ���� H
    for m = 1: maxSize - 1
        % �� ���� ������ Hy(n) ������ �������� ���������� Hy
        % �� ���������� ������ �������
        Hy(m) = Hy(m) + (Ez(m + 1) - Ez(m)) / (W0 * mu(m)) * Sc;
    end
    
    % �������� ����������� � �������������� ������ 
    % Total Field / Scattered Field
    Hy(sourcePos - 1) = Hy(sourcePos - 1) - ...
        Sc / (W0 * mu(sourcePos - 1)) * ...
        exp (-(t - 30.0) ^ 2 / 100.0);
    
    % ������ ���������� ���� E
    for m = 2: maxSize
        % �� ���� ������ Ez(n) ������ �������� ���������� EzS
        % �� ���������� ������ �������
        Ez(m) = Ez(m) + (Hy(m) - Hy(m - 1)) * (W0 / eps(m)) * Sc;
    end

    % �������� ����������� � �������������� ������ 
    % Total Field / Scattered Field
    Ez(sourcePos) = Ez(sourcePos) + ...
      Sc / (sqrt(eps(sourcePos) * mu(sourcePos))) *...
      exp(-(t + 0.5 - (-0.5 * sqrt(eps(sourcePos) * mu(sourcePos)) / Sc)...
      - 30.0) ^ 2 / 100.0);
    
    % ����������� ���� � �����
    probeTimeEz(t) = Ez(probePos);
    
    plot (Ez);
    xlim ([1, maxSize]);
    ylim ([-1.1, 1.1]);
    xlabel ('x, ������')
    ylabel ('Ez, �/�')
    title (sprintf('%d', t))
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