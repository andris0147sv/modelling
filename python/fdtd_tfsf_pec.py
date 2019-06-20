# -*- coding: utf-8 -*-
'''
Метод полного поля / рассеянного поля. Две границы.
Моделирование распространения ЭМ волны, пащающей на диэлектрический слой.
'''

import numpy

import tools


if __name__ == '__main__':
    # Волновое сопротивление свободного пространства
    W0 = 120.0 * numpy.pi

    # Число Куранта
    Sc = 1.0

    # Время расчета в отсчетах
    maxTime = 450

    # Размер области моделирования в отсчетах
    maxSize = 250

    # Датчики для регистрации поля
    probesPos = [60]
    probes = [tools.Probe(pos, maxTime) for pos in probesPos]

    # Левая граница TFSF
    tfsf_left = 50

    # Правая граница TFSF
    tfsf_right = 150

    # Положение металлической плоскости (PEC)
    layer_x = 100

    # Параметры среды
    # Диэлектрическая проницаемость
    eps = numpy.ones(maxSize)

    # Магнитная проницаемость
    mu = numpy.ones(maxSize - 1)

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
    display.drawSources([tfsf_left, tfsf_right])
    display.drawBoundary(layer_x)

    for t in range(maxTime):
        # Расчет компоненты поля H
        Hy = Hy + (Ez[1:] - Ez[:-1]) * Sc / (W0 * mu)

        # Источник возбуждения с использованием метода
        # Total Field / Scattered Field
        Hy[tfsf_left - 1] -= numpy.exp(-(t - 30.0 - (tfsf_left - tfsf_left)) ** 2 / 100.0) / W0
        Hy[tfsf_right - 1] += numpy.exp(-(t - 30.0 - (tfsf_right - tfsf_left)) ** 2 / 100.0) / W0

        # Граничные условия для поля E
        Ez[0] = Ez[1]
        Ez[-1] = Ez[-2]

        # Расчет компоненты поля E
        Hy_shift = Hy[:-1]
        Ez[1:-1] = Ez[1:-1] + (Hy[1:] - Hy_shift) * Sc * W0 / eps[1:-1]

        # Источник возбуждения с использованием метода
        # Total Field / Scattered Field
        Ez[tfsf_left] += numpy.exp(-(t + 0.5 - (tfsf_left - tfsf_left - 0.5) - 30.0) ** 2 / 100.0)
        Ez[tfsf_right] -= numpy.exp(-(t + 0.5 - (tfsf_right - tfsf_left - 0.5) - 30.0) ** 2 / 100.0)

        Ez[layer_x] = 0

        # Регистрация поля в датчиках
        for probe in probes:
            probe.addData(Ez, Hy)

        if t % 2 == 0:
            display.updateData(display_field, t)

    display.stop()

    # Отображение сигнала, сохраненного в пробнике
    tools.showProbeSignals(probes, -1.1, 1.1)
