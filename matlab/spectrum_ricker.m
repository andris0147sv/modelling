%% Расчет спектра гауссова импульса
clear

% Размер массива
size = 256;

% шаг по времени
dt = 1e-10;

% Шаг по частоте
df = 1.0 / (size * dt);

fp = 1e9;
Md = 2;
dr = Md / fp

% Вейвлет Рикера
time = (1:size) * dt;
impulse = (1 - 2 * (pi * fp * (time - dr)) .^ 2) .*...
    exp (-(pi * fp * (time - dr)) .^ 2);

% Расчет спектра
spectrum = fft(impulse);
spectrum = fftshift (spectrum);

% Расчет частоты
freq = (-size / 2:size / 2 - 1) * df;

% Отображение импульса
subplot (2, 1, 1)
plot (time, impulse)
grid on
xlabel ('Время, с')
ylabel ('Ez')

% Отображение спектра
subplot (2, 1, 2)
plot (freq, abs (spectrum))
grid on
xlabel ('Частота, Гц')
ylabel ('|P|')
set(gca,'XTick',[-5e9: 1e9: 5e9])