# -*- coding: utf-8 -*-
'''
ЭМ волна в виде модулированного гассова импульса распространяется
в одну сторону (TFSF boundary).
Источник находится в диэлектрике.
'''

import numpy

import tools


class GaussianModPlaneWave:
    ''' Класс с уравнением плоской волны для модулированного гауссова сигнала в дискретном виде
    d - определяет задержку сигнала.
    w - определяет ширину сигнала.
    Nl - количество ячеек на длину волны.
    Sc - число Куранта.
    eps - относительная диэлектрическая проницаемость среды, в которой расположен источник.
    mu - относительная магнитная проницаемость среды, в которой расположен источник.
    '''

    def __init__(self, d, w, Nl, Sc=1.0, eps=1.0, mu=1.0):
        self.d = d
        self.w = w
        self.Nl = Nl
        self.Sc = Sc
        self.eps = eps
        self.mu = mu

    def getE(self, m, q):
        '''
        Расчет поля E в дискретной точке пространства m
        в дискретный момент времени q
        '''
        return (numpy.sin(2 * numpy.pi / self.Nl * (q * self.Sc - m * numpy.sqrt(self.eps * self.mu))) *
                numpy.exp(-(((q - m * numpy.sqrt(self.eps * self.mu) / self.Sc) - self.d) / self.w) ** 2))


if __name__ == '__main__':
    # Волновое сопротивление свободного пространства
    W0 = 120.0 * numpy.pi

    # Число Куранта
    Sc = 1.0

    # Время расчета в отсчетах
    maxTime = 500

    # Размер области моделирования в отсчетах
    maxSize = 200

    # Положение источника в отсчетах
    sourcePos = 50

    # Параметры сигнала
    dg = 150
    wg = 50
    Nl = 40

    # Датчики для регистрации поля
    probesPos = [60]
    probes = [tools.Probe(pos, maxTime) for pos in probesPos]

    # Параметры среды
    # Диэлектрическая проницаемость
    eps = numpy.ones(maxSize)
    eps[:] = 9.0

    # Магнитная проницаемость
    mu = numpy.ones(maxSize - 1)

    Ez = numpy.zeros(maxSize)
    Hy = numpy.zeros(maxSize - 1)

    source = GaussianModPlaneWave(
        dg, wg, Nl, Sc, eps[sourcePos], mu[sourcePos])

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

    for q in range(maxTime):
        # Расчет компоненты поля H
        Hy = Hy + (Ez[1:] - Ez[:-1]) * Sc / (W0 * mu)

        # Источник возбуждения с использованием метода
        # Total Field / Scattered Field
        Hy[sourcePos - 1] -= Sc / (W0 * mu[sourcePos - 1]) * source.getE(0, q)

        # Граничные условия для поля E - PEC

        # Расчет компоненты поля E
        Hy_shift = Hy[:-1]
        Ez[1:-1] = Ez[1:-1] + (Hy[1:] - Hy_shift) * Sc * W0 / eps[1:-1]

        # Источник возбуждения с использованием метода
        # Total Field / Scattered Field
        Ez[sourcePos] += (Sc / (numpy.sqrt(eps[sourcePos] * mu[sourcePos])) *
                          source.getE(-0.5, q + 0.5))

        # Регистрация поля в датчиках
        for probe in probes:
            probe.addData(Ez, Hy)

        if q % 2 == 0:
            display.updateData(display_field, q)

    display.stop()

    # Отображение сигнала, сохраненного в датчиках
    tools.showProbeSignals(probes, -1.1, 1.1)
