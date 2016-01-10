%% Одномерный FDTD. Версия 1.8
% Граница раздела. Используются граничные условия ABC первого порядка.
clear

% Волновое сопротивление свободного пространства
W0 = 120 * pi;

% Число Куранта
Sc = 1;

% Время расчета в отсчетах
maxTime = 550;

% Размер области моделирования в отсчетах
maxSize = 200;

% Положение датчика, регистрирующего поля
probePos = 50;

layer_x = 100;

Ez = zeros (1, maxSize);
Hy = zeros (1, maxSize - 1);

eps = ones (size (Ez));
eps(layer_x: end) = 9.0;

mu = ones (size (Ez));

% !!!
% Ez(2) в предыдущий момент времени
oldEzLeft = 0;

% !!!
% Ez(end-1) в предыдущий момент времени
oldEzRight = 0;

% !!!
% Расчет коэффициентов для граничных условий
tempLeft = Sc / sqrt (mu(1) * eps(1));
koeffABCLeft = (tempLeft - 1) / (tempLeft + 1);

tempRight = Sc / sqrt (mu(end) * eps(end));
koeffABCRight = (tempRight - 1) / (tempRight + 1);

% Поле, зарегистрированное в датчике в зависимости от времени
probeTimeEz = zeros (1, maxTime);

figure

for t = 1: maxTime
    % Расчет компоненты поля H
    for m = 1: maxSize - 1
        % До этой строки Hy(n) хранит значение компоненты Hy
        % за предыдущий момент времени
        Hy(m) = Hy(m) + (Ez(m + 1) - Ez(m)) / W0 / mu(m);
    end
    
    Hy(49) = Hy(49) - exp (-(t - 30.0) ^ 2 / 100.0) / W0;
    
    % Расчет компоненты поля E
  
    for m = 2: maxSize - 1
        % До этой строки Ez(n) хранит значение компоненты EzS
        % за предыдущий момент времени
        Ez(m) = Ez(m) + (Hy(m) - Hy(m - 1)) * W0 / eps (m);
    end

    % Источник возбуждения
    Ez(50) = Ez(50) + exp (-(t + 0.5 - (-0.5) - 30.0) ^ 2 / 100.0);
    
    % Граничные условия ABC первой степени
    Ez(1) = oldEzLeft + koeffABCLeft * (Ez(2) - Ez(1));
    oldEzLeft = Ez(2);
    
    Ez(end) = oldEzRight + koeffABCRight * (Ez(end-1) - Ez(end));
    oldEzRight = Ez(end-1);
    
    % Регистрация поля в точке
    probeTimeEz(t) = Ez(probePos);
    
    plot (Ez);
    xlim ([1, maxSize]);
    ylim ([-1.1, 1.1]);
    xlabel ('x, отсчет')
    ylabel ('Ez, В/м')
    line ([layer_x, layer_x], [-1.1, 1.1], ...
        'Color',[0.0, 0.0, 0.0]);
    pause (0.01)
end

figure
plot (probeTimeEz)
xlabel ('t, отсчет')
ylabel ('Ez, В/м')
grid on