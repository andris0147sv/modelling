# -*- coding: utf-8 -*-
'''
Гауссов импульс распространяется в свободном пространстве.
Sc < 1
'''

import numpy
from numpy.fft import fft
import matplotlib.pyplot as plt

import tools


if __name__ == '__main__':
    # Волновое сопротивление свободного пространства
    W0 = 120.0 * numpy.pi

    c = 299792458.0

    # Число Куранта
    Sc = 0.9

    # Время расчета в отсчетах
    maxTime = 1024

    # Размер области моделирования в отсчетах
    maxSize = 2000

    # Положение источника в отсчетах
    sourcePos = 1200

    # Положения датчиков
    probe1Pos = 1300
    # probe2Pos = 1350
    probe2Pos = 1500

    # Расстояние между датчиками
    probeDist = probe2Pos - probe1Pos

    # Датчики для регистрации поля
    probesPos = [probe1Pos, probe2Pos]
    probes = [tools.Probe(pos, maxTime) for pos in probesPos]

    Ez = numpy.zeros(maxSize)
    Hy = numpy.zeros(maxSize)

    # Параметры отображения поля
    # Для поля E
    display_field = Ez
    display_ylabel = 'Ez, В/м'
    display_ymin = -0.2
    display_ymax = 0.6

    # Создание экземпляра класса для отображения
    # распределения поля в пространстве
    display = tools.AnimateFieldDisplay(maxSize,
                                        display_ymin, display_ymax,
                                        display_ylabel)

    display.activate()
    display.drawSources([sourcePos])
    display.drawProbes(probesPos)

    # Номер индекса в спектре, в котором считаем, что фаза меняется линейно
    k_index = 10

    for t in range(maxTime):
        # Расчет компоненты поля H
        Ez_shift = Ez[1:]
        Hy[:-1] = Hy[:-1] + (Ez_shift - Ez[:-1]) * Sc / W0

        # Расчет компоненты поля E
        Hy_shift = Hy[:-1]
        Ez[1:] = Ez[1:] + (Hy[1:] - Hy_shift) * Sc * W0

        # Источник возбуждения
        Ez[sourcePos] += numpy.exp(-(t - 30.0) ** 2 / (5.0 ** 2)) * Sc

        # Регистрация поля в датчиках
        for probe in probes:
            probe.addData(Ez, Hy)

        if t % 10 == 0:
            display.updateData(display_field, t)

    display.stop()

    # Отображение сигналов и их спектров, зарегистрированных в пробниках
    plt.figure()
    plt.suptitle('Сигналы в датчиках')

    for n, probe in enumerate(probes):
        EzField = probe.E
        spectrum = fft(EzField)
        spectrum_abs = numpy.abs(spectrum)
        spectrum_phase = numpy.unwrap(numpy.angle(spectrum))

        # Расчет фазы с вычтенным линейным членом
        k = spectrum_phase[k_index] / k_index
        spectrum_phase_line = spectrum_phase - \
            k * numpy.arange(0, len(spectrum_phase))

        # Отображение сигнала в датчике
        plt.subplot(4, 1, 1)
        plt.plot(EzField, label='Пробник N {}'.format(n + 1))
        plt.xlabel('t, отсчет')
        plt.ylabel('Ez, В/м')
        plt.grid()
        plt.legend()

        # Амплитудный спектр сигнала в датчике
        plt.subplot(4, 1, 2)
        plt.plot(spectrum_abs, label='Пробник N {}'.format(n + 1))
        plt.xlabel('f')
        plt.ylabel('|Ez|, В/(м*Гц)')
        plt.xlim(0, 300)
        plt.grid()
        plt.legend()

        # Фазовый спектр сигнала в датчике без вычитания линейного члена
        plt.subplot(4, 1, 3)
        plt.plot(spectrum_phase, label='Пробник N {}'.format(n + 2))
        plt.xlabel('f')
        plt.ylabel('Phase(Ez), рад.')
        plt.xlim(0, 300)
        plt.grid()
        plt.legend()

        # Фазовый спектр сигнала в датчике
        plt.subplot(4, 1, 4)
        plt.plot(spectrum_phase_line, label='Пробник N {}'.format(n + 2))
        plt.xlabel('f')
        plt.ylabel('Phase(Ez), рад.')
        plt.xlim(0, 300)
        plt.ylim(-5, 5)
        plt.grid()
        plt.legend()

    # Сигнал и спектр в первом пробнике
    EzField_0 = probes[0].E
    spectrum_0 = fft(EzField_0)
    spectrum_0_phase = numpy.unwrap(numpy.angle(spectrum_0))

    # Сигнал и спектр во втором пробнике
    EzField_1 = probes[1].E
    spectrum_1 = fft(EzField_1)
    spectrum_1_phase = numpy.unwrap(numpy.angle(spectrum_1))

    phase_delta = numpy.abs(spectrum_1_phase - spectrum_0_phase)

    # Коэффициент прохождения среды между двумя пробниками
    R = spectrum_1 / spectrum_0
    R_abs = numpy.abs(R)
    R_phase = numpy.unwrap(numpy.angle(R))

    # Вычитание линейного члена из фазы коэффициента прохождения
    k = R_phase[k_index] / k_index
    R_phase_line = R_phase - k * numpy.arange(0, len(R))

    # Отображение спектра падающего сигнала
    plt.figure()
    plt.suptitle('Коэффициент прохождения')
    plt.subplot(3, 1, 1)
    plt.plot(numpy.abs(spectrum_0))
    plt.xlabel('f')
    plt.ylabel('|Ez|, В/(м*Гц)')
    plt.xlim(0, 300)
    plt.grid()

    # Отображение амплитуды коэффициента прохождения
    plt.subplot(3, 1, 2)
    plt.plot(R_abs)
    plt.xlabel('f')
    plt.ylabel('|R|')
    plt.xlim(0, 300)
    plt.ylim(0, 1.1)
    plt.grid()

    # Отображение фазы коэффициента прохождения
    plt.subplot(3, 1, 3)
    plt.plot(R_phase_line)
    plt.xlabel('f')
    plt.ylabel('Phase(R), рад.')
    plt.xlim(0, 300)
    plt.ylim(-4, 4)
    plt.grid()

    # Расчет фазовой скорости на каждой частоте
    v = (c * probeDist * 2 * numpy.pi *
         numpy.arange(0, len(spectrum_0)) / (Sc * phase_delta * len(EzField_1)))

    # Отображение спектра падающего сигнала
    plt.figure()
    plt.suptitle('Фазовая скорость')
    plt.subplot(2, 1, 1)
    plt.plot(numpy.abs(spectrum_0))
    plt.xlabel('f')
    plt.ylabel('|Ez|, В/(м*Гц)')
    plt.xlim(0, 300)
    plt.grid()

    plt.subplot(2, 1, 2)
    plt.plot(v)
    plt.xlabel('f')
    plt.ylabel('v, м/c')
    plt.xlim(0, 300)
    plt.ylim(c * 0.9, c * 1.1)
    plt.grid()

    plt.show()
