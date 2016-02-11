%% Двумерный FDTD. Версия 1.0
% Не используются граничные условия
clear
clc

%% Параметры области моделирования
% Шаг сетки (d = dx = dy)
d = 1e-3;

% Время моделирования в секундах
maxTime_sec = 1e-9;

% Размер области моделирования в метрах
sizeX_m = 0.2;
sizeY_m = 0.2;

% Положение точечного источника в метрах
port_x_m = 0.1;
port_y_m = 0.1;

% Параметры гауссова сигнала
gauss_width_sec = 2e-11;
gauss_delay_sec = 2.5 * gauss_width_sec;

%% Физические константы
% Магнитная постоянная
mu0 = pi * 4e-7;

% Электрическая постоянная
eps0 = 8.854187817e-12;

% Скорость света в вакууме
c = 1.0 / sqrt (mu0 * eps0);


%% Расчет "дискретных" параметров моделирования
% "Одномерный" аналог числа Куранта для случая 2D
Cdtds = 1.0 / sqrt (2.0);

dt = d / (c) * Cdtds;

% Волновое сопротивление свободного пространства
W0 = 120 * pi;

% Время расчета в отсчетах
maxTime = ceil (maxTime_sec / dt);

% Размер области моделирования в отсчетах
sizeX = ceil (sizeX_m / d);
sizeY = ceil (sizeY_m / d);

% Положение точки возбуждения
port_x = ceil (port_x_m / d);
port_y = ceil (port_y_m / d);

gauss_width = gauss_width_sec / dt;
gauss_delay = gauss_delay_sec / dt;


%% Параметры среды
% Диэлектрическая проницаемость среды
eps = ones (sizeX, sizeY);

% Магнитная проницаемость среды
mu = ones (sizeX, sizeY);

% Проводимость среды
sigma = zeros (sizeX, sizeY);

% "Магнитная проводимость" среды
sigma_m = zeros (sizeX, sizeY);

%% Коэффициенты для конечно-разностной схемы
Chxh = (1 - sigma_m .* dt ./ (2 * mu * mu0)) ./ ...
    (1 + sigma_m .* dt ./ (2 * mu * mu0));

Chxe = 1 ./ (1 + (sigma_m .* dt ./ (2 * mu * mu0))) .*...
    dt ./ (mu * mu0 * d);

Chyh = Chxh;

Chye = Chxe;

Ceze = (1 - sigma .* dt ./ (2 * eps * eps0)) ./ ...
    (1 + sigma .* dt ./ (2 * eps * eps0));

Cezh = 1 ./ (1 + (sigma .* dt ./ (2 * eps * eps0))) .*...
    dt ./ (eps * eps0 * d);

%% Компоненты поля
Hx = zeros (sizeX, sizeY);
Hy = zeros (sizeX, sizeY);

Ez = zeros (sizeX, sizeY);

figure;
[x, y] = meshgrid (1:sizeX, 1:sizeY);
x = x';
y = y';

%% Конечно-разностная схема
for t = 1: maxTime        
    for m = 1:sizeX
        for n = 1:sizeY - 1
            Hx(m, n) = Chxh(m, n) * Hx(m, n) -...
                Chxe(m, n) * (Ez(m, n + 1) - Ez(m, n));
        end
    end

    for m = 1:sizeX - 1
        for n = 1:sizeY
            Hy(m, n) = Chyh(m, n) * Hy(m, n) +...
                Chye(m, n) * (Ez(m + 1, n) - Ez(m, n));
        end
    end

    for m = 2:sizeX - 1
        for n = 2:sizeY - 1
            Ez(m, n) = Ceze(m, n) * Ez(m, n) +...
                Cezh(m, n) * ((Hy(m, n) - Hy(m - 1, n)) -...
                (Hx(m, n) - Hx(m, n-1)));
        end
    end
    
    Ez(port_x, port_y) = exp (-(t - gauss_delay) ^ 2 / (gauss_width ^ 2));
    
    %surfl(x, y, Ez);
    %shading interp;

    imagesc(Ez, [-1, 1]);
    colormap jet;
    
    zlim([-1.1, 1.1])
    pause (0.01);
end

