# -*- coding: utf-8 -*-
'''
Двумерный FDTD. Версия 1.0
Поляризация TMz. Граничные условия - PEC
'''

import pylab
import numpy

if __name__ == '__main__':
    # Физические константы
    # Магнитная постоянная
    mu0 = numpy.pi * 4e-7

    # Электрическая постоянная
    eps0 = 8.854187817e-12

    # Скорость света в вакууме
    c = 1.0 / numpy.sqrt(mu0 * eps0)

    # Шаг сетки (d = dx = dy)
    d = 1e-3

    # Время моделирования в секундах
    maxTime_sec = 1.1e-9

    # Размер области моделирования в метрах
    sizeX_m = 0.3
    sizeY_m = 0.2

    # Положение точечного источника в метрах
    port_x_m = 0.15
    port_y_m = 0.1

    # Положение пробника в метрах
    probe_x_m = 0.12
    probe_y_m = 0.08

    # Частота источника
    f = 10e9
    phi_0 = 0
    wavelength_m = c / f
    wavelength = wavelength_m / d

    # Расчет "дискретных" параметров моделирования
    # "Одномерный" аналог числа Куранта для случая 2D
    Cdtds = 1.0 / numpy.sqrt(2.0)

    dt = d / c * Cdtds

    # Волновое сопротивление свободного пространства
    W0 = 120 * numpy.pi

    # Время расчета в отсчетах
    maxTime = int(numpy.ceil(maxTime_sec / dt))

    # Размер области моделирования в отсчетах
    sizeX = int(numpy.ceil(sizeX_m / d))
    sizeY = int(numpy.ceil(sizeY_m / d))

    # Положение точки возбуждения
    port_x = int(numpy.ceil(port_x_m / d))
    port_y = int(numpy.ceil(port_y_m / d))

    # Положение пробника
    probe_x = int(numpy.ceil(probe_x_m / d))
    probe_y = int(numpy.ceil(probe_y_m / d))

    # Компоненты поля
    Hx = numpy.zeros((sizeX, sizeY))
    Hy = numpy.zeros((sizeX, sizeY))
    Ez = numpy.zeros((sizeX, sizeY))

    # Параметры среды
    # Диэлектрическая проницаемость среды
    eps = numpy.ones((sizeX, sizeY))

    # Магнитная проницаемость среды
    mu = numpy.ones((sizeX, sizeY))

    # Проводимость среды
    sigma = numpy.zeros((sizeX, sizeY))

    # "Магнитная проводимость" среды
    sigma_m = numpy.zeros((sizeX, sizeY))

    # Коэффициенты для конечно-разностной схемы
    loss_m = sigma_m * dt / (2 * mu * mu0)
    loss = sigma * dt / (2 * eps * eps0)

    Chxh = ((1 - sigma_m * dt / (2 * mu * mu0)) /
            (1 + sigma_m * dt / (2 * mu * mu0)))

    Chxe = 1 / (1 + (sigma_m * dt / (2 * mu * mu0))) * dt / (mu * mu0 * d)

    Chyh = Chxh
    Chye = Chxe

    Ceze = ((1 - sigma * dt / (2 * eps * eps0)) /
            (1 + sigma * dt / (2 * eps * eps0)))

    Cezh = (1 / (1 + (sigma * dt / (2 * eps * eps0))) *
            dt / (eps * eps0 * d))

    # Какую компоненту поля будем отображать
    visualize_field = Ez

    # Поле, зарегистрированное в датчике в зависимости от времени
    probeTimeHx = numpy.zeros(maxTime)
    probeTimeHy = numpy.zeros(maxTime)
    probeTimeEz = numpy.zeros(maxTime)

    pylab.ion()
    fig = pylab.figure()

    for t in range(maxTime):
        Hx[:, :-1] = (Chxh[:, :-1] * Hx[:, :-1] -
                      Chxe[:, :-1] * (Ez[:, 1:] - Ez[:, :-1]))

        Hy[:-1, :] = (Chyh[:-1, :] * Hy[:-1, :] +
                      Chye[:-1, :] * (Ez[1:, :] - Ez[:-1, :]))

        Ez[1:-1, 1:-1] = (Ceze[1:-1, 1:-1] * Ez[1:-1, 1:-1] +
                          Cezh[1:-1, 1:-1] * ((Hy[1:-1, 1:-1] - Hy[:-2, 1:-1]) -
                                              (Hx[1:-1, 1:-1] - Hx[1:-1, :-2])))

        Ez[port_x, port_y] = (Ez[port_x, port_y] +
                              numpy.sin(2 * numpy.pi * t / wavelength + phi_0))

        probeTimeHx[t] = Hx[probe_x, probe_y]
        probeTimeHy[t] = Hy[probe_x, probe_y]
        probeTimeEz[t] = Ez[probe_x, probe_y]

        if t % 2 == 0:
            pylab.clf()
            pylab.imshow(visualize_field.transpose(),
                         vmin=-0.1,
                         vmax=0.1,
                         cmap='jet')
            pylab.draw()
            pylab.pause(0.01)

    pylab.ioff()
    pylab.figure()
    pylab.subplot(3, 1, 1)
    pylab.plot(probeTimeHx, 'b')
    pylab.xlabel('t, отсчет')
    pylab.ylabel('Hx, А/м')
    pylab.grid()

    pylab.subplot(3, 1, 2)
    pylab.plot(probeTimeHy, 'b')
    pylab.xlabel('t, отсчет')
    pylab.ylabel('Hy, А/м')
    pylab.grid()

    pylab.subplot(3, 1, 3)
    pylab.plot(probeTimeEz, 'r')
    pylab.xlabel('t, отсчет')
    pylab.ylabel('Ez, В/м')
    pylab.grid()

    pylab.show()