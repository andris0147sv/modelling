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

    xlist = numpy.arange(maxSize)
    pylab.ion()

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
            pylab.clf()
            pylab.plot(xlist, Ez)
            pylab.xlim(0, maxSize)
            pylab.ylim(-1.1, 1.1)
            pylab.xlabel('t, отсчет')
            pylab.ylabel('Ez, В/м')
            pylab.draw()
            pylab.pause(0.01)

    pylab.close()
