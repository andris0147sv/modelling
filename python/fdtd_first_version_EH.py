# -*- coding: utf-8 -*-
'''
История изменений:
    * Добавлена возможность переключения отображения поля E или H в анимации.
    * Функции и классы, не связанные напрямую с методом FDTD, вынесены в модуль tools.
    * Работа с пробниками вынесена в класс Probe.
    * Рисование графиков вынесено в класс AnimateFieldDisplay и функцию showProbeSignals.
    * Циклы по пространству заменены на поэлементные операции с массивами.
    * Добавлена анимация поля E.
    * Первая версия программы с использованием метода FDTD.
'''

import numpy

import tools


if __name__ == '__main__':
    # Волновое сопротивление свободного пространства
    W0 = 120.0 * numpy.pi

    # Число Куранта
    Sc = 1.0

    # Время расчета в отсчетах
    maxTime = 1000

    # Размер области моделирования в отсчетах
    maxSize = 200

    # Датчики для регистрации поля
    probesPos = [50, 100]
    probes = [tools.Probe(pos, maxTime) for pos in probesPos]

    Ez = numpy.zeros(maxSize)
    Hy = numpy.zeros(maxSize)

    # Параметры отображения поля
    # Для поля E
    # display_field = Ez
    # display_ylabel = 'Ez, В/м'
    # display_ymin = -1.1
    # display_ymax = 1.1

    # Для поля H
    display_field = Hy
    display_ylabel = 'Hy, А/м'
    display_ymin = -1.1 / W0
    display_ymax = 1.1 / W0

    # Создание экземпляра класса для отображения
    # распределения поля в пространстве
    display = tools.AnimateFieldDisplay(maxSize,
                                        display_ymin, display_ymax,
                                        display_ylabel)

    display.activate()
    display.drawProbes(probesPos)

    for t in range(maxTime):
        # Расчет компоненты поля H
        Ez_shift = Ez[1:]
        Hy[:-1] = Hy[:-1] + (Ez_shift - Ez[:-1]) * Sc / W0

        # Расчет компоненты поля E
        Hy_shift = Hy[:-1]
        Ez[1:] = Ez[1:] + (Hy[1:] - Hy_shift) * Sc * W0

        # Источник возбуждения
        Ez[0] = numpy.exp(-(t - 30.0) ** 2 / 100.0)

        # Регистрация поля в датчиках
        for probe in probes:
            probe.addData(Ez, Hy)

        if t % 2 == 0:
            display.updateData(display_field, t)

    display.stop()

    # Отображение сигнала, сохраненного в пробнике
    tools.showProbeSignals(probes, -1.1, 1.1)