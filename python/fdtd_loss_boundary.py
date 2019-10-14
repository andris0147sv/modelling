# -*- coding: utf-8 -*-
'''
Моделирование распространения ЭМ волны, падающей на границу
вакуум - диэлектрик с без потерь.
Справа - поглощающий полностью согласованный слой.
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
    probesPos = [60]
    probes = [tools.Probe(pos, maxTime) for pos in probesPos]

    # Положение начала диэлектрика без потерь
    layer_x = 100

    # Где начинается поглощающий диэлектрик
    layer_loss_x = 160

    # Параметры среды
    # Диэлектрическая проницаемость
    eps = numpy.ones(maxSize)
    eps[layer_x:] = 9.0

    # Магнитная проницаемость
    mu = numpy.ones(maxSize - 1)

    # Потери в среде. loss = sigma * dt / (2 * eps * eps0)
    loss = numpy.zeros(maxSize)
    loss[layer_loss_x:] = 0.02

    # Коэффициенты для расчета поля E
    ceze = (1 - loss) / (1 + loss)
    cezh = W0 / (eps * (1 + loss))

    # Коэффициенты для расчета поля H
    chyh = (1 - loss) / (1 + loss)
    chye = 1 / (W0 * (1 + loss))

    Ez = numpy.zeros(maxSize)
    Hy = numpy.zeros(maxSize - 1)

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
    display.drawBoundary(layer_loss_x)

    for t in range(maxTime):
        # Расчет компоненты поля H
        Hy = chyh[:-1] * Hy + chye[:-1] * (Ez[1:] - Ez[:-1])

        # Источник возбуждения с использованием метода
        # Total Field / Scattered Field
        Hy[sourcePos - 1] -= numpy.exp(-(t - 30.0) ** 2 / 100.0) / W0

        # Граничные условия для поля E
        Ez[0] = Ez[1]
        Ez[-1] = Ez[-2]

        # Расчет компоненты поля E
        Hy_shift = Hy[:-1]
        Ez[1:-1] = ceze[1: -1] * Ez[1:-1] + cezh[1: -1] * (Hy[1:] - Hy_shift)

        # Источник возбуждения с использованием метода
        # Total Field / Scattered Field
        Ez[sourcePos] += numpy.exp(-((t + 0.5) - (-0.5) - 30.0) ** 2 / 100.0)

        # Регистрация поля в датчиках
        for probe in probes:
            probe.addData(Ez, Hy)

        if t % 2 == 0:
            display.updateData(display_field, t)

    display.stop()

    # Отображение сигнала, сохраненного в пробнике
    tools.showProbeSignals(probes, -1.1, 1.1)
