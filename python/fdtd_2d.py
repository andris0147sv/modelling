# -*- coding: utf-8 -*-
'''
Двумерный FDTD. Версия 1.0
Поляризация TEz. Граничные условия - PEC
'''

import pylab
import numpy

if __name__ == '__main__':
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

    # Параметры гауссова сигнала
    gauss_width_sec = 2e-11
    gauss_delay_sec = 2.5 * gauss_width_sec

    # Физические константы
    # Магнитная постоянная
    mu0 = numpy.pi * 4e-7

    # Электрическая постоянная
    eps0 = 8.854187817e-12

    # Скорость света в вакууме
    c = 1.0 / numpy.sqrt(mu0 * eps0)

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

    gauss_width = gauss_width_sec / dt
    gauss_delay = gauss_delay_sec / dt

    # Компоненты поля
    Ex = numpy.zeros((sizeX, sizeY))
    Ey = numpy.zeros((sizeX, sizeY))
    Hz = numpy.zeros((sizeX, sizeY))

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

    Chzh = (1 - loss_m) / (1 + loss_m)

    Chze = 1 / (1 + loss_m) * (dt / (mu * mu0 * d))

    Cexe = (1 - loss) / (1 + loss)

    Cexh = 1 / (1 + loss) * (dt / (eps * eps0 * d))

    Ceye = Cexe
    Ceyh = Cexh

    # Какую компоненту поля будем отображать
    visualize_field = Ey

    # Поле, зарегистрированное в датчике в зависимости от времени
    probeTimeEx = numpy.zeros(maxTime)
    probeTimeEy = numpy.zeros(maxTime)
    probeTimeHz = numpy.zeros(maxTime)

    pylab.ion()
    fig = pylab.figure()

    for t in range(maxTime):
        Hz[:-1, :-1] = (Chzh[:-1, :-1] * Hz[:-1, :-1] +
                        Chze[:-1, :-1] * (Ex[:-1, 1:] - Ex[:-1, :-1] -
                                          (Ey[1:, :-1] - Ey[:-1, :-1])))

        Ex[:-1, 1:-1] = (Cexe[:-1, 1:-1] * Ex[:-1, 1:-1] +
                         Cexh[:-1, 1:-1] * (Hz[:-1, 1:-1] - Hz[:-1, :-2]))

        Ey[1:-1, :-1] = (Ceye[1:-1, :-1] * Ey[1:-1, :-1] -
                         Ceyh[1:-1, :-1] * (Hz[1:-1, :-1] - Hz[:-2, :-1]))

        Hz[port_x, port_y] = Hz[port_x, port_y] + numpy.exp(-(t - gauss_delay) ** 2 / (gauss_width ** 2))

        probeTimeEx[t] = Ex[probe_x, probe_y]
        probeTimeEy[t] = Ey[probe_x, probe_y]
        probeTimeHz[t] = Hz[probe_x, probe_y]

        if t % 2 == 0:
            pylab.clf()
            pylab.imshow(Ex.transpose(), vmin=-10.0, vmax=10.0, cmap='jet')
            pylab.draw()
            pylab.pause(0.01)

    pylab.ioff()
    pylab.figure()
    pylab.subplot(3, 1, 1)
    pylab.plot(probeTimeEx, 'b')
    pylab.xlabel('t, отсчет')
    pylab.ylabel('Ex, В/м')
    pylab.grid()

    pylab.subplot(3, 1, 2)
    pylab.plot(probeTimeEy, 'b')
    pylab.xlabel('t, отсчет')
    pylab.ylabel('Ey, В/м')
    pylab.grid()

    pylab.subplot(3, 1, 3)
    pylab.plot(probeTimeHz, 'r')
    pylab.xlabel('t, отсчет')
    pylab.ylabel('Hz, А/м')
    pylab.grid()

    pylab.show()
