# -*- coding: utf-8 -*-
'''
Модуль с классами источников разного типа
'''

from abc import ABCMeta, abstractmethod

import numpy as np


class Source1D(metaclass=ABCMeta):
    '''
    Базовый класс для всех источников одномерного метода FDTD
    '''
    @abstractmethod
    def getField(self, time):
        '''
        Метод должен возвращать значение поля источника в момент времени time
        '''
        pass


class Gaussian(Source1D):
    '''
    Источник, создающий гауссов импульс
    '''

    def __init__(self, magnitude, dg, wg):
        '''
        magnitude - максимальное значение в источнике;
        dg - коэффициент, задающий начальную задержку гауссова импульса;
        wg - коэффициент, задающий ширину гауссова импульса;
        '''
        self._magnitude = magnitude
        self._dg = dg
        self._wg = wg

    def getField(self, time):
        return self._magnitude * np.exp(-((time - self._dg) / self._wg) ** 2)
