# -*- coding: utf-8 -*-

import pylab
import numpy
from typing import List


class AnimateFieldDisplay:
    '''
    Класс для отображения анимации распространения ЭМ волны в пространстве
    '''
    def __init__(self,
                 maxXSize: int,
                 minYSize: float, maxYSize: float,
                 probePos: List[int]):
        '''
        maxXSize - размер области моделирования в отсчетах.
        minYSize, maxYSize - интервал отображения графика по оси Y.
        probePos - список координат пробников для регистрации временных
            сигналов.
        '''
        self.maxXSize = maxXSize
        self.minYSize = minYSize
        self.maxYSize = maxYSize
        self.probePos = probePos
        self._xList = None
        self._line = None
        self._xlabel = 'x, отсчет'
        self._ylabel = 'Ez, В/м'
        self._probeStyle = 'xr'

    def activate(self):
        '''
        Инициализировать окно с анимацией
        '''
        self._xList = numpy.arange(maxSize)

        # Включить интерактивный режим для анимации
        pylab.ion()

        # Создание окна для графика
        self._fig, self._ax = pylab.subplots()

        # Установка отображаемых интервалов по осям
        self._ax.set_xlim(0, self.maxXSize)
        self._ax.set_ylim(self.minYSize, self.maxYSize)

        # Установка меток по осям
        self._ax.set_xlabel(self._xlabel)
        self._ax.set_ylabel(self._ylabel)

        # Включить сетку на графике
        self._ax.grid()

        # Отобразить поле в начальный момент времени
        self._line, = self._ax.plot(self._xList, numpy.zeros(self.maxXSize))

        # Отобразить положение пробника
        self._ax.plot(self.probePos,
                      [0] * len(self.probePos),
                      self._probeStyle)

    def stop(self):
        '''
        Остановить анимацию
        '''
        pylab.ioff()

    def updateData(self, data):
        '''
        Обновить данные с распределением поля в пространстве
        '''
        self._line.set_ydata(data)
        self._fig.canvas.draw()
        self._fig.canvas.flush_events()


def showProbeSignals(signals: List[List[float]],
                     minYSize: float,
                     maxYSize: float):
    '''
    Показать графики сигналов, зарегистрированых в датчиках.

    signals - список сигналов, зарегистрированных датчиками.
    minYSize, maxYSize - интервал отображения графика по оси Y.
    '''
    fig, ax = pylab.subplots()
    ax.set_xlim(0, len(signals[0]))
    ax.set_ylim(minYSize, maxYSize)
    ax.set_xlabel('t, отсчет')
    ax.set_ylabel('Ez, В/м')
    ax.grid()

    for signal in signals:
        ax.plot(signal)

    pylab.show()


if __name__ == '__main__':
    # Волновое сопротивление свободного пространства
    W0 = 120.0 * numpy.pi

    # Число Куранта
    Sc = 1.0

    # Время расчета в отсчетах
    maxTime = 1000

    # Размер области моделирования в отсчетах
    maxSize = 200

    # Положение датчика, регистрирующего поля
    probePos = 50

    Ez = numpy.zeros(maxSize)
    Hy = numpy.zeros(maxSize)

    # Поле, зарегистрированное в датчике в зависимости от времени
    probeTimeEz = numpy.zeros(maxTime)

    # Создание экземпляра класса для отображения
    # распределения поля в пространстве
    display = AnimateFieldDisplay(maxSize, -1.1, 1.1, [probePos])
    display.activate()

    for t in range(maxTime):
        # Расчет компоненты поля H
        Ez_shift = Ez[1:]
        Hy[:-1] = Hy[:-1] + (Ez_shift - Ez[:-1]) * Sc / W0

        # Расчет компоненты поля E
        Hy_shift = Hy[:-1]
        Ez[1:] = Ez[1:] + (Hy[1:] - Hy_shift) * Sc * W0

        # Источник возбуждения
        Ez[0] = numpy.exp(-(t - 30.0) ** 2 / 100.0)

        # Регистрация поля в точке
        probeTimeEz[t] = Ez[probePos]

        if t % 2 == 0:
            display.updateData(Ez)

    display.stop()

    # Отображение сигнала, сохраненного в пробнике
    showProbeSignals([probeTimeEz], -1.1, 1.1)
