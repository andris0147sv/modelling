# -*- coding: utf-8 -*-
'''
История изменений:
    * Первая версия программы с использованием метода FDTD.
'''

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

    # Положение источника
    sourcePos = 75

    Ez = numpy.zeros(maxSize)
    Hy = numpy.zeros(maxSize)

    # Поле, зарегистрированное в датчике в зависимости от времени
    probeTimeEz = numpy.zeros(maxTime)
    probeTimeEz[0] = Ez[probePos]

    for t in range(1, maxTime):
        # Расчет компоненты поля H
        for x in range(0, maxSize - 1):
            Hy[x] = Hy[x] + (Ez[x + 1] - Ez[x]) * Sc / W0

        # Расчет компоненты поля E
        for x in range(1, maxSize):
            Ez[x] = Ez[x] + (Hy[x] - Hy[x - 1]) * Sc * W0

        # Источник возбуждения
        Ez[sourcePos] += numpy.exp(-(t - 0.5 - 30.0) ** 2 / 100.0)

        # Регистрация поля в точке
        probeTimeEz[t] = Ez[probePos]

    # Отображение сигнала, сохраненного в датчике
    tlist = numpy.arange(maxTime)
    fig, ax = pylab.subplots()
    ax.set_xlim(0, maxTime)
    ax.set_ylim(-1.1, 1.1)
    ax.set_xlabel('t, отсчет')
    ax.set_ylabel('Ez, В/м')
    ax.plot(tlist, probeTimeEz)
    ax.grid()
    pylab.show()
