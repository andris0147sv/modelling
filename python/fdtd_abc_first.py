# -*- coding: utf-8 -*-
'''
Моделирование распространения ЭМ волны, падающей на границу
вакуум - идеальный диэлектрик.
Используются граничные условия ABC первой степени.
'''

import numpy

import tools


if __name__ == '__main__':
    # Волновое сопротивление свободного пространства
    W0 = 120.0 * numpy.pi

    # Число Куранта
    Sc = 1.0

    # Время расчета в отсчетах
    maxTime = 750

    # Размер области моделирования в отсчетах
    maxSize = 200

    # Положение источника в отсчетах
    sourcePos = 50

    # Датчики для регистрации поля
    probesPos = [160]
    probes = [tools.Probe(pos, maxTime) for pos in probesPos]

    # Положение начала диэлектрика
    layer_x = 100

    # Параметры среды
    # Диэлектрическая проницаемость
    eps = numpy.ones(maxSize)
    eps[layer_x:] = 9.0

    # Магнитная проницаемость
    mu = numpy.ones(maxSize - 1)

    Ez = numpy.zeros(maxSize)
    Hy = numpy.zeros(maxSize - 1)

    # Ez[1] в предыдущий момент времени
    oldEzLeft = Ez[1]

    # Ez[-2] в предыдущий момент времени
    oldEzRight = Ez[-2]

    # Расчет коэффициентов для граничных условий
    tempLeft = Sc / numpy.sqrt(mu[0] * eps[0])
    koeffABCLeft = (tempLeft - 1) / (tempLeft + 1)

    tempRight = Sc / numpy.sqrt(mu[-1] * eps[-1])
    koeffABCRight = (tempRight - 1) / (tempRight + 1)

    # Параметры отображения поля E
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
    display.drawBoundary(layer_x)

    for q in range(maxTime):
        # Расчет компоненты поля H
        Hy = Hy + (Ez[1:] - Ez[:-1]) * Sc / (W0 * mu)

        # Источник возбуждения с использованием метода
        # Total Field / Scattered Field
        Hy[sourcePos - 1] -= (Sc / (W0 * mu[sourcePos - 1]) *
                              numpy.exp(-(q - 30.0) ** 2 / 100.0))

        # Расчет компоненты поля E
        Hy_shift = Hy[: -1]
        Ez[1:-1] = Ez[1: -1] + (Hy[1:] - Hy_shift) * Sc * W0 / eps[1: -1]

        # Источник возбуждения с использованием метода
        # Total Field / Scattered Field
        Ez[sourcePos] += (Sc / (numpy.sqrt(eps[sourcePos] * mu[sourcePos])) *
                          numpy.exp(-((q + 0.5) - (-0.5 * numpy.sqrt(eps[sourcePos] * mu[sourcePos]) / Sc) - 30.0) ** 2 / 100.0))

        # Граничные условия ABC первой степени
        Ez[0] = oldEzLeft + koeffABCLeft * (Ez[1] - Ez[0])
        oldEzLeft = Ez[1]

        Ez[-1] = oldEzRight + koeffABCRight * (Ez[-2] - Ez[-1])
        oldEzRight = Ez[-2]

        # Регистрация поля в датчиках
        for probe in probes:
            probe.addData(Ez, Hy)

        if q % 2 == 0:
            display.updateData(display_field, q)

    display.stop()

    # Отображение сигнала, сохраненного в датчиках
    tools.showProbeSignals(probes, -1.1, 1.1)
