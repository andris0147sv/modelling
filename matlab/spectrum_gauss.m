%% Расчет спектра гауссова импульса

% Размер массива
size = 256;

% шаг по времени
dt = 1e-10;

% Шаг по частоте
df = 1.0 / (size * dt);

% Гауссов импульс
time = 1:size;
gauss = exp (-(time - 30.0) .^ 2 / 100.0);

% Расчет спектра
spectrum = fft(gauss);
spectrum = fftshift (spectrum);

% Расчет частоты
freq = (-size / 2:size / 2 - 1) * df;

% Отображение импульса
subplot (2, 1, 1)
plot (time * dt, gauss)
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