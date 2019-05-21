# -*- coding: utf-8 -*-

import pylab
import numpy

if __name__ == '__main__':
    # Волновое сопротивление свободного пространства
    W0 = 120.0 * numpy.pi

    # Число Куранта
    Sc = 1.0

    # Время расчета в отсчетах
    maxTime = 1000

    # Размер области моделирования в отсчетах
    maxSize = 200

    # Положение датчика, регистрирующего поля
    probePos = 50

    Ez = numpy.zeros(maxSize)
    Hy = numpy.zeros(maxSize)

    # Поле, зарегистрированное в датчике в зависимости от времени
    probeTimeEz = numpy.zeros(maxTime)

    # Подготовка к отображению поля в пространстве
    xlist = numpy.arange(maxSize)

    # Включить интерактивный режим для анимации
    pylab.ion()

    # Создание окна для графика
    fig, ax = pylab.subplots()

    # Установка отображаемых интервалов по осям
    ax.set_xlim(0, maxSize)
    ax.set_ylim(-1.1, 1.1)

    # Установка меток по осям
    ax.set_xlabel('x, отсчет')
    ax.set_ylabel('Ez, В/м')

    # Включить сетку на графике
    ax.grid()

    # Отобразить поле в начальный момент времени
    line, = ax.plot(xlist, Ez)

    # Отобразить положение пробника
    ax.plot(probePos, 0, 'xr')

    for t in range(maxTime):
        # Расчет компоненты поля H
        for x in range(0, maxSize - 1):
            Hy[x] = Hy[x] + (Ez[x + 1] - Ez[x]) * Sc / W0

        # Расчет компоненты поля E
        for x in range(1, maxSize):
            Ez[x] = Ez[x] + (Hy[x] - Hy[x - 1]) * Sc * W0

        # Источник возбуждения
        Ez[0] = numpy.exp(-(t - 30.0) ** 2 / 100.0)

        # Регистрация поля в точке
        probeTimeEz[t] = Ez[probePos]

        if t % 2 == 0:
            # Обновить данные на графике
            line.set_ydata(Ez)
            fig.canvas.draw()
            fig.canvas.flush_events()

    # Отключить интерактивный режим по завершению анимации
    pylab.ioff()

    # Отображение сигнала, сохраненного в пробнике
    tlist = numpy.arange(maxTime)
    fig, ax = pylab.subplots()
    ax.set_xlim(0, maxTime)
    ax.set_ylim(-1.1, 1.1)
    ax.set_xlabel('t, отсчет')
    ax.set_ylabel('Ez, В/м')
    ax.plot(tlist, probeTimeEz)
    ax.grid()
    pylab.show()
