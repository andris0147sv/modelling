%% Расчет спектра гауссова импульса

% Размер массива
size = 256;

% шаг по времени
dt = 1e-10;

A_0 = 100;
A_max = 100;
F_max = 1e9;

% Шаг по частоте
df = 1.0 / (size * dt);

w_g = sqrt(log(A_max)) / (pi * F_max);
d_g = w_g * log (A_0);

% Гауссов импульс
time = (1:size) * dt;
gauss = exp (-((time - d_g) / w_g) .^ 2);

% Расчет спектра
spectrum = fft(gauss);
spectrum = fftshift (spectrum);

% Расчет частоты
freq = (-size / 2:size / 2 - 1) * df;

% Отображение импульса
subplot (2, 1, 1)
plot (time, gauss)
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