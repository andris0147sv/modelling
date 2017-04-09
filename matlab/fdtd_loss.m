%% Одномерный FDTD. Версия 1.5
% Среда с потерями.
clear

% Волновое сопротивление свободного пространства
W0 = 120 * pi;

% Потери в среде. loss = sigma * dt / (2 * eps * eps0)
loss = 0.01;

% Время расчета в отсчетах
maxTime = 1000;

% Размер области моделирования в отсчетах
maxSize = 200;

% Положение датчика для регистрации поля
probePos = 60;

% Начало слоя с потерями
layer_x = 100;

Ez = zeros (1, maxSize);
Hy = zeros (1, maxSize - 1);

eps = ones (size (Ez));
eps(layer_x: end) = 9.0;

ceze = ones (1, maxSize);
ceze(layer_x: end) = (1 - loss) / (1 + loss);

cezh = ones (1, maxSize) * W0 ./ eps;
cezh(layer_x: end) = cezh(layer_x: end) / (1 + loss);

% Поле, зарегистрированное в датчике в зависимости от времени
probeTimeEz = zeros (1, maxTime);

figure

for t = 1: maxTime
    % Расчет компоненты поля H
    for m = 1: maxSize - 1
        % До этой строки Hy(n) хранит значение компоненты Hy
        % за предыдущий момент времени
        Hy(m) = Hy(m) + (Ez(m + 1) - Ez(m)) / W0;
    end
    
    Hy(49) = Hy(49) - exp (-(t - 30.0) ^ 2 / 100.0) / W0;
    
    % Расчет компоненты поля E
    % Т.е. пока справа простейшее граничное условие и так не работает,
    % правое граничное условие убрано.
    Ez(1) = Ez(2);
    
    for m = 2: maxSize - 1
        % До этой строки Ez(n) хранит значение компоненты EzS
        % за предыдущий момент времени
        % Вместо W0 / eps используется cezh
        Ez(m) = ceze(m) * Ez(m) + ...
            cezh(m) * (Hy(m) - Hy(m - 1));
    end

    % Источник возбуждения
    Ez(50) = Ez(50) + exp (-(t + 0.5 - (-0.5) - 30.0) ^ 2 / 100.0);
    
    % Регистрация поля в точке
    probeTimeEz(t) = Ez(probePos);
    
    plot (Ez);
    xlim ([1, maxSize]);
    ylim ([-1.1, 1.1]);
    xlabel ('x, отсчет')
    ylabel ('Ez, В/м')
    line ([layer_x, layer_x], [-1.1, 1.1], ...
        'Color',[0.0, 0.0, 0.0]);
    grid on
    pause (0.03)
end

figure
plot (probeTimeEz)
xlabel ('t, отсчет')
ylabel ('Ez, В/м')
grid on