# -*- coding: utf-8 -*-
'''
Моделирование отражения гармонического сигнала от слоя диэлектрика
'''

import math
from typing import List

import numpy

import tools
import sources
import boundary
from objects import LayerContinuous, LayerDiscrete, Probe


class Sampler:
    def __init__(self, discrete: float):
        self.discrete = discrete

    def sample(self, x: float) -> int:
        return math.floor(x / self.discrete + 0.5)


def sampleLayer(layer_cont: LayerContinuous, sampler: Sampler) -> LayerDiscrete:
    start_discrete = sampler.sample(layer_cont.xmin)
    end_discrete = (sampler.sample(layer_cont.xmax)
                    if layer_cont.xmax is not None
                    else None)
    return LayerDiscrete(start_discrete, end_discrete,
                         layer_cont.eps, layer_cont.mu)


def fillMedium(layer: LayerDiscrete, eps: List[float], mu: List[float]):
    if layer.xmax is not None:
        eps[layer.xmin: layer.xmax] = layer.eps
        mu[layer.xmin: layer.xmax] = layer.mu
    else:
        eps[layer.xmin:] = layer.eps
        mu[layer.xmin:] = layer.mu


if __name__ == '__main__':
    # Используемые константы
    # Волновое сопротивление свободного пространства
    W0 = 120.0 * numpy.pi

    # Скорость света в вакууме
    c = 299792458.0

    # Электрическая постоянная
    eps0 = 8.854187817e-12

    # Параметры моделирования
    # Частота сигнала, Гц
    f_Hz = 1.2e9

    # Дискрет по пространству в м
    dx = 2e-3

    wavelength = c / f_Hz
    Nl = wavelength / dx

    # Число Куранта
    Sc = 1.0

    # Размер области моделирования в м
    maxSize_m = 2.2

    # Время расчета в секундах
    maxTime_s = 20e-9

    # Положение источника в м
    sourcePos_m = 0.8

    # Координаты датчиков для регистрации поля в м
    probesPos_m = [0.3, 0.8, 1.3]

    # Параметры слоев
    layers_cont = [LayerContinuous(1.1, eps=4.4)]

    # Скорость обновления графика поля
    speed_refresh = 15

    # Переход к дискретным отсчетам
    # Дискрет по времени
    dt = dx * Sc / c

    sampler_x = Sampler(dx)
    sampler_t = Sampler(dt)

    # Время расчета в отсчетах
    maxTime = sampler_t.sample(maxTime_s)

    # Размер области моделирования в отсчетах
    maxSize = sampler_x.sample(maxSize_m)

    # Положение источника в отсчетах
    sourcePos = sampler_x.sample(sourcePos_m)

    layers = [sampleLayer(layer, sampler_x) for layer in layers_cont]

    # Датчики для регистрации поля
    probesPos = [sampler_x.sample(pos) for pos in probesPos_m]
    probes = [Probe(pos, maxTime) for pos in probesPos]

    # Вывод параметров моделирования
    print(f'Число Куранта: {Sc}')
    print(f'Размер области моделирования: {maxSize_m} м')
    print(f'Время расчета: {maxTime_s * 1e9} нс')
    print(f'Координата источника: {sourcePos_m} м')
    print(f'Частота сигнала: {f_Hz * 1e-9} ГГц')
    print(f'Длина волны: {wavelength} м')
    print(f'Количество отсчетов на длину волны (Nl): {Nl}')
    probes_m_str = ', '.join(['{:.6f}'.format(pos) for pos in probesPos_m])
    print(f'Дискрет по пространству: {dx} м')
    print(f'Дискрет по времени: {dt * 1e9} нс')
    print(f'Координаты пробников [м]: {probes_m_str}')
    print()
    print(f'Размер области моделирования: {maxSize} отсч.')
    print(f'Время расчета: {maxTime} отсч.')
    print(f'Координата источника: {sourcePos} отсч.')
    probes_str = ', '.join(['{}'.format(pos) for pos in probesPos])
    print(f'Координаты пробников [отсч.]: {probes_str}')

    # Параметры среды
    # Диэлектрическая проницаемость
    eps = numpy.ones(maxSize)

    # Магнитная проницаемость
    mu = numpy.ones(maxSize - 1)

    for layer in layers:
        fillMedium(layer, eps, mu)

    # Источник
    # amp = 2.5e-3 * dt / (eps0 * eps[sourcePos] * dx)
    amp = 1.0
    amp_E = amp * Sc / numpy.sqrt(eps[sourcePos] * mu[sourcePos])
    amp_H = amp * Sc / (W0 * mu[sourcePos])

    source_E = sources.make_harmonic(amp_E, f_Hz, Sc, dt)
    source_H = sources.make_harmonic(amp_H, f_Hz, Sc, dt)

    Ez = numpy.zeros(maxSize)
    Hy = numpy.zeros(maxSize - 1)

    # Создание экземпляров классов граничных условий
    boundary_left = boundary.ABCSecondLeft(eps[0], mu[0], Sc)
    boundary_right = boundary.ABCSecondRight(eps[-1], mu[-1], Sc)

    # Параметры отображения поля E
    display_field = Ez
    display_ylabel = 'Ez, В/м'
    display_ymin = -2.1
    display_ymax = 2.1

    # Создание экземпляра класса для отображения
    # распределения поля в пространстве
    display = tools.AnimateFieldDisplay(dx, dt,
                                        maxSize,
                                        display_ymin, display_ymax,
                                        display_ylabel,
                                        title='fdtd_dielectric')

    display.activate()
    display.drawSources([sourcePos])
    display.drawProbes(probesPos)
    for layer in layers:
        display.drawBoundary(layer.xmin)
        if layer.xmax is not None:
            display.drawBoundary(layer.xmax)

    for t in range(1, maxTime):
        # Расчет компоненты поля H
        Hy = Hy + (Ez[1:] - Ez[:-1]) * Sc / (W0 * mu)

        # Источник возбуждения с использованием метода
        # Total Field / Scattered Field
        Hy[sourcePos - 1] -= source_H.getField(t)

        # Расчет компоненты поля E
        Ez[1:-1] = Ez[1: -1] + (Hy[1:] - Hy[: -1]) * Sc * W0 / eps[1: -1]

        # Источник возбуждения с использованием метода
        # Total Field / Scattered Field
        Ez[sourcePos] += source_E.getField(t + 0.5 - (-0.5))

        boundary_left.updateField(Ez, Hy)
        boundary_right.updateField(Ez, Hy)

        # Регистрация поля в датчиках
        for probe in probes:
            probe.addData(Ez, Hy)

        if t % speed_refresh == 0:
            display.updateData(display_field, t)

    display.stop()

    # Отображение сигнала, сохраненного в пробнике
    tools.showProbeSignals(probes, dx, dt, -2.1, 2.1)
