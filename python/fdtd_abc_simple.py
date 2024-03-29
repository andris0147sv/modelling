# -*- coding: utf-8 -*-
'''
Простейшие граничные условия
'''

import numpy

import tools


if __name__ == '__main__':
    # Волновое сопротивление свободного пространства
    W0 = 120.0 * numpy.pi

    # Число Куранта
    Sc = 1.0

    # Время расчета в отсчетах
    maxTime = 400

    # Размер области моделирования в отсчетах
    maxSize = 200

    # Положение источника в отсчетах
    sourcePos = 75

    # Датчики для регистрации поля
    probesPos = [50, 100]
    probes = [tools.Probe(pos, maxTime) for pos in probesPos]

    Ez = numpy.zeros(maxSize)
    Hy = numpy.zeros(maxSize)

    for probe in probes:
        probe.addData(Ez, Hy)

    # Параметры отображения поля
    # Для поля E
    display_field = Ez
    display_ylabel = 'Ez, В/м'
    display_ymin = -1.1
    display_ymax = 1.1

    # Создание экземпляра класса для отображения
    # распределения поля в пространстве
    display = tools.AnimateFieldDisplay(maxSize,
                                        display_ymin, display_ymax,
                                        display_ylabel)

    display.activate()
    display.drawProbes(probesPos)
    display.drawSources([sourcePos])

    for q in range(1, maxTime):
        # Граничные условия для поля H
        # H(q + 1/2)[-1] = H(q - 1/2)[-2]
        Hy[-1] = Hy[-2]

        # Расчет компоненты поля H в момент времени q + 1/2
        Ez_shift = Ez[1:]
        Hy[:-1] = Hy[:-1] + (Ez_shift - Ez[:-1]) * Sc / W0

        # Граничные условия для поля E
        # E(q+1)[0] = E(q)[1]
        Ez[0] = Ez[1]

        # Расчет компоненты поля E в момент времени q + 1
        Hy_shift = Hy[:-1]
        Ez[1:] = Ez[1:] + (Hy[1:] - Hy_shift) * Sc * W0

        # Источник возбуждения
        Ez[sourcePos] += numpy.exp(-(q - 0.5 - 30.0) ** 2 / 100.0)

        # Регистрация поля в датчиках
        for probe in probes:
            probe.addData(Ez, Hy)

        if q % 2 == 0:
            display.updateData(display_field, q)

    display.stop()

    # Отображение сигнала, сохраненного в датчиках
    tools.showProbeSignals(probes, -1.1, 1.1)
