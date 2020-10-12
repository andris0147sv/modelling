# -*- coding: utf-8 -*-
'''
Метод полного поля / рассеянного поля (TF / SF).
Гауссов импульс.
Две границы TF / SF.
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

    # Левая граница TF/SF
    tfsf_left = 50

    # Правая граница TF/SF
    tfsf_right = 150

    # Датчики для регистрации поля
    probesPos = [75]
    probes = [tools.Probe(pos, maxTime) for pos in probesPos]

    Ez = numpy.zeros(maxSize)
    Hy = numpy.zeros(maxSize)

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

    for q in range(maxTime):
        # Граничные условия для поля H
        Hy[-1] = Hy[-2]

        # Расчет компоненты поля H
        Ez_shift = Ez[1:]
        Hy[:-1] = Hy[:-1] + (Ez_shift - Ez[:-1]) * Sc / W0

        # Источник возбуждения с использованием метода
        # Total Field / Scattered Field
        Hy[tfsf_left - 1] -= (Sc / W0) * numpy.exp(-(q - 30.0 - tfsf_left) ** 2 / 100.0)
        Hy[tfsf_right - 1] += (Sc / W0) * numpy.exp(-(q - 30.0 - tfsf_right) ** 2 / 100.0)

        # Hy[tfsf_left - 1] -= (Sc / W0) * numpy.exp(-(q - 30.0 - (tfsf_left - tfsf_left)) ** 2 / 100.0)
        # Hy[tfsf_right - 1] += (Sc / W0) * numpy.exp(-(q - 30.0 - (tfsf_right - tfsf_left)) ** 2 / 100.0)

        # Граничные условия для поля E
        Ez[0] = Ez[1]

        # Расчет компоненты поля E
        Hy_shift = Hy[:-1]
        Ez[1:] = Ez[1:] + (Hy[1:] - Hy_shift) * Sc * W0

        # Источник возбуждения с использованием метода
        # Total Field / Scattered Field
        Ez[tfsf_left] += Sc * numpy.exp(-(q + 0.5 - (tfsf_left - 0.5) - 30.0) ** 2 / 100.0)
        Ez[tfsf_right] -= Sc * numpy.exp(-(q + 0.5 - (tfsf_right - 0.5) - 30.0) ** 2 / 100.0)

        # Ez[tfsf_left] += Sc * numpy.exp(-(q + 0.5 - (tfsf_left - tfsf_left - 0.5) - 30.0) ** 2 / 100.0)
        # Ez[tfsf_right] -= Sc * numpy.exp(-(q + 0.5 - (tfsf_right - tfsf_left - 0.5) - 30.0) ** 2 / 100.0)

        # Регистрация поля в датчиках
        for probe in probes:
            probe.addData(Ez, Hy)

        if q % 2 == 0:
            display.updateData(display_field, q)

    display.stop()

    # Отображение сигнала, сохраненного в датчиках
    tools.showProbeSignals(probes, -1.1, 1.1)
