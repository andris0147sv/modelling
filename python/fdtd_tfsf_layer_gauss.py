# -*- coding: utf-8 -*-
'''
Метод полного поля / рассеянного поля. Две границы.
Моделирование распространения ЭМ волны, падающей на диэлектрический слой.
'''

import numpy

import tools


if __name__ == '__main__':
    # Волновое сопротивление свободного пространства
    W0 = 120.0 * numpy.pi

    # Число Куранта
    Sc = 1.0

    # Время расчета в отсчетах
    maxTime = 2000

    # Размер области моделирования в отсчетах
    maxSize = 250

    # Датчики для регистрации поля
    probesPos = [25]
    probes = [tools.Probe(pos, maxTime) for pos in probesPos]

    # Левая граница TFSF
    tfsf_left = 50

    # Правая граница TFSF
    tfsf_right = 200

    # Положение начала диэлектрика
    layer_start = 100
    layer_end = 150

    # Параметры среды
    # Диэлектрическая проницаемость
    eps = numpy.ones(maxSize)
    eps[layer_start: layer_end] = 9.0

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
    display.drawBoundary(layer_start)
    display.drawBoundary(layer_end)

    for t in range(maxTime):
        # Расчет компоненты поля H
        Hy = Hy + (Ez[1:] - Ez[:-1]) * Sc / (W0 * mu)

        # Источник возбуждения с использованием метода
        # Total Field / Scattered Field
        Hy[tfsf_left - 1] -= Sc / (W0 * mu[tfsf_left - 1]) * numpy.exp(-(t - 30.0 - (tfsf_left - tfsf_left)) ** 2 / 100.0)
        Hy[tfsf_right - 1] += Sc / (W0 * mu[tfsf_right - 1]) * numpy.exp(-(t - 30.0 - (tfsf_right - tfsf_left)) ** 2 / 100.0)

        # Граничные условия для поля E
        Ez[0] = Ez[1]
        Ez[-1] = Ez[-2]

        # Расчет компоненты поля E
        Hy_shift = Hy[:-1]
        Ez[1:-1] = Ez[1:-1] + (Hy[1:] - Hy_shift) * Sc * W0 / eps[1:-1]

        # Источник возбуждения с использованием метода
        # Total Field / Scattered Field
        Ez[tfsf_left] += Sc / (numpy.sqrt(eps[tfsf_left] * mu[tfsf_left])) * numpy.exp(-(t + 0.5 - (tfsf_left - tfsf_left - 0.5) - 30.0) ** 2 / 100.0)
        Ez[tfsf_right] -= Sc / (numpy.sqrt(eps[tfsf_right] * mu[tfsf_right])) * numpy.exp(-(t + 0.5 - (tfsf_right - tfsf_left - 0.5) - 30.0) ** 2 / 100.0)

        # Регистрация поля в датчиках
        for probe in probes:
            probe.addData(Ez, Hy)

        if t % 2 == 0:
            display.updateData(display_field, t)

    display.stop()

    # Отображение сигнала, сохраненного в пробнике
    tools.showProbeSignals(probes, -1.1, 1.1)
