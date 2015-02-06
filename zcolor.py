#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created by Zoltan Bozoky on 2013.02.09.
Under GPL licence.

Purpose:
========
* An easy handle color managament
* Color schemes
* Generate color sets
"""

listnames = ['4seasons', 'OnlineCulture']


paletta8 = ((219,94,86),
            (219,194,86),
            (145,219,86),
            (86,219,127),
            (86,211,219),
            (86,111,219),
            (160,86,219),
            (219,86,178))

paletta10 = ((231, 47, 39),
            (238, 113, 25),
            (255, 200, 8),
            (170, 189, 27),
            (19, 166, 50),
            (4, 148, 87),
            (1, 134, 141),
            (3, 86, 155),
            (46, 20, 141),
            (204, 63, 92))

heatmap1 = ((10, 50, 120),
            (15, 75, 165),
            (30, 110, 200),
            (60, 160, 240),
            (80, 180, 250),
            (130, 210, 255),
            (160, 240, 255),
            (200, 250, 255),
            (230, 255, 255),
            (255, 250, 220),
            (255, 232, 120),
            (255, 192, 60),
            (255, 160, 0),
            (255, 96, 0),
            (255, 50, 0),
            (225, 20, 0),
            (192, 0, 0),
            (165, 0, 0))



def _4seasons(number):
    szinek = ['DE1841', 'F5861A', 'F5DE1A', '5DD517', '1680CA', '6016CA',
              'F6451A', 'F6B41A', 'D6EB19', '15CA7F', '152ACA', 'CA16BD']
    return szinek[number % len(szinek)]

def _OnlineCulture(number):
    szinek = ['8CC63F', '0A803B', '25AAE1', '0F75BC', 'F7941E', 'F1592A',
              'ED217C', 'BF1E2D']
    return szinek[number % len(szinek)]

def Colorlist(listname, number):
    if listname in listnames:
        szinkod = eval('_'+listname+'('+str(number)+')')
        szinek = []
        for i in range(3):
            szinek.append(eval('0x'+szinkod[i*2:i*2+2]))
    return szinek


def colorwheel_constants(colorname):
    red =   [100, 100, 100, 100, 100, 100, 100, 100, 100, 100,  #   0 -   9
             100, 100, 100, 100, 100, 100, 100, 100, 100, 100,  #  10 -  19
             100, 100, 100, 100, 100, 100, 100, 100, 100, 100,  #  20 -  29
             100, 100, 100, 100, 100, 100, 100, 100, 100, 100,  #  30 -  39
             100, 100, 100, 100, 100, 100, 100, 100, 100, 100,  #  40 -  49
             100, 100, 100, 100, 100, 100, 100, 100, 100, 100,  #  50 -  59
             100, 100, 100, 100, 100, 100, 100, 100, 100, 100,  #  60 -  69
             100, 100, 100, 100, 100, 100, 100, 100, 100, 100,  #  70 -  79
             100, 100, 100, 100, 100, 100, 100, 100, 100, 100,  #  80 -  89
             100, 100, 100, 100, 100, 100, 100, 100, 100, 100,  #  90 -  99
             100, 100, 100, 100, 100, 100, 100, 100, 100, 100,  # 100 - 109
             100, 100, 100, 100, 100, 100, 100, 100, 100, 100,  # 110 - 119
              99,  97,  96,  95,  94,  93,  92,  91,  89,  88,  # 120 - 129
              87,  86,  85,  84,  83,  82,  81,  80,  79,  78,  # 130 - 139
              76,  75,  75,  73,  72,  71,  70,  69,  68,  67,  # 140 - 149
              65,  64,  63,  62,  60,  59,  58,  56,  55,  54,  # 150 - 159
              52,  51,  49,  47,  45,  44,  42,  39,  37,  35,  # 160 - 169
              33,  30,  27,  24,  21,  17,  14,   9,   5,   0,  # 170 - 179
               0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  # 180 - 189
               0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  # 190 - 199
               0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  # 200 - 209
               0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  # 210 - 219
               0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  # 220 - 229
               0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  # 230 - 239
               0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  # 240 - 249
               0,   0,   0,   0,   0,   2,   4,   5,   7,   9,  # 250 - 259
              11,  12,  14,  15,  17,  18,  20,  21,  22,  24,  # 260 - 269
              25,  26,  28,  29,  30,  31,  33,  34,  35,  36,  # 270 - 279
              38,  39,  40,  42,  43,  44,  45,  47,  48,  49,  # 280 - 289
              51,  52,  54,  55,  56,  58,  60,  61,  63,  64,  # 290 - 299
              66,  68,  70,  72,  74,  76,  78,  80,  83,  85,  # 300 - 309
              88,  91,  93,  96, 100, 100, 100, 100, 100, 100,  # 310 - 319
             100, 100, 100, 100, 100, 100, 100, 100, 100, 100,  # 320 - 329
             100, 100, 100, 100, 100, 100, 100, 100, 100, 100,  # 330 - 339
             100, 100, 100, 100, 100, 100, 100, 100, 100, 100,  # 340 - 349
             100, 100, 100, 100, 100, 100, 100, 100, 100]       # 350 - 359
    #######################################################################
    green = [  3,   5,   7,   9,  12,  14,  16,  17,  19,  21,  #   0 -   9
              22,  24,  25,  27,  29,  30,  31,  33,  34,  35,  #  10 -  19
              36,  37,  38,  39,  40,  42,  42,  44,  44,  45,  #  20 -  29
              46,  47,  48,  49,  50,  51,  51,  52,  53,  54,  #  30 -  39
              55,  55,  56,  56,  57,  58,  58,  59,  60,  60,  #  40 -  49
              61,  62,  62,  63,  64,  64,  65,  65,  66,  67,  #  50 -  59
              67,  68,  68,  69,  69,  70,  71,  71,  72,  72,  #  60 -  69
              73,  73,  74,  75,  75,  75,  76,  76,  77,  78,  #  70 -  79
              78,  79,  79,  80,  80,  81,  81,  82,  82,  83,  #  80 -  89
              84,  84,  84,  85,  85,  86,  87,  87,  87,  88,  #  90 -  99
              89,  89,  90,  91,  91,  91,  92,  93,  93,  94,  # 100 - 109
              95,  95,  96,  96,  96,  97,  98,  99,  99, 100,  # 110 - 119
             100, 100, 100, 100, 100, 100, 100, 100, 100, 100,  # 120 - 129
             100, 100, 100, 100, 100, 100, 100, 100, 100, 100,  # 130 - 139
             100, 100, 100, 100, 100, 100, 100, 100, 100, 100,  # 140 - 149
             100, 100, 100, 100, 100, 100, 100, 100, 100, 100,  # 150 - 159
             100, 100, 100, 100, 100, 100, 100, 100, 100, 100,  # 160 - 169
             100, 100, 100, 100, 100, 100, 100, 100, 100, 100,  # 170 - 179
             100, 100, 100, 100, 100, 100, 100, 100, 100, 100,  # 180 - 189
             100, 100, 100, 100, 100, 100, 100, 100, 100, 100,  # 190 - 199
             100, 100, 100, 100, 100, 100, 100, 100, 100, 100,  # 200 - 209
              96,  91,  88,  84,  81,  78,  75,  73,  70,  67,  # 210 - 219
              65,  63,  61,  58,  56,  55,  53,  51,  49,  47,  # 220 - 229
              45,  44,  42,  40,  39,  37,  35,  34,  32,  30,  # 230 - 239
              29,  27,  25,  23,  22,  20,  18,  16,  14,  12,  # 240 - 249
               9,   7,   5,   3,   0,   0,   0,   0,   0,   0,  # 250 - 259
               0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  # 260 - 269
               0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  # 270 - 279
               0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  # 280 - 289
               0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  # 290 - 299
               0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  # 300 - 309
               0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  # 310 - 319
               0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  # 320 - 329
               0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  # 330 - 339
               0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  # 340 - 349
               0,   0,   0,   0,   0,   0,   0,   0,   0]       # 350 - 359
    #######################################################################
    blue = [   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  #   0 -   9
               0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  #  10 -  19
               0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  #  20 -  29
               0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  #  30 -  39
               0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  #  40 -  49
               0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  #  50 -  59
               0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  #  60 -  69
               0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  #  70 -  79
               0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  #  80 -  89
               0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  #  90 -  99
               0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  # 100 - 109
               0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  # 110 - 119
               0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  # 120 - 129
               0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  # 130 - 139
               0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  # 140 - 149
               0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  # 150 - 159
               0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  # 160 - 169
               0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  # 170 - 179
               7,  12,  17,  22,  26,  30,  34,  37,  40,  44,  # 180 - 189
              46,  49,  52,  55,  57,  60,  62,  65,  67,  70,  # 190 - 199
              72,  75,  78,  80,  83,  86,  89,  93,  96, 100,  # 200 - 209
             100, 100, 100, 100, 100, 100, 100, 100, 100, 100,  # 210 - 219
             100, 100, 100, 100, 100, 100, 100, 100, 100, 100,  # 220 - 229
             100, 100, 100, 100, 100, 100, 100, 100, 100, 100,  # 230 - 239
             100, 100, 100, 100, 100, 100, 100, 100, 100, 100,  # 240 - 249
             100, 100, 100, 100, 100, 100, 100, 100, 100, 100,  # 250 - 259
             100, 100, 100, 100, 100, 100, 100, 100, 100, 100,  # 260 - 269
             100, 100, 100, 100, 100, 100, 100, 100, 100, 100,  # 270 - 279
             100, 100, 100, 100, 100, 100, 100, 100, 100, 100,  # 280 - 289
             100, 100, 100, 100, 100, 100, 100, 100, 100, 100,  # 290 - 299
             100, 100, 100, 100, 100, 100, 100, 100, 100, 100,  # 300 - 309
             100, 100, 100, 100, 100,  96,  91,  88,  84,  81,  # 310 - 319
              78,  75,  73,  70,  67,  65,  63,  61,  58,  56,  # 320 - 329
              55,  53,  51,  49,  47,  45,  44,  42,  40,  39,  # 330 - 339
              37,  35,  34,  32,  30,  29,  27,  25,  23,  22,  # 340 - 349
              20,  18,  16,  14,  12,   9,   7,   5,   3]       # 350 - 359
    #######################################################################
    return eval(colorname)




class MyColor(object):
    """
    Everything which is related to colors
    """
    __R = colorwheel_constants('red')
    __G = colorwheel_constants('green')
    __B = colorwheel_constants('blue')
    def __init__(self):
        """
        Initializing the color wheel and constants
        """
        self.__lowest = 0.0
        self.__highest = 100.0

        self.__triadangle = 30.0
        self.__tetradangle = 30.0
        self.__analogicangle = 30.0
        return None
    def rgb(self, angle, brightness = None):
        """
        Get one color rgb values
        Parameters:
        -----------
            * angle: 0 - 360 deg (red - green - blue)
            * brightness: 0 - 100 (dark to light)
        Return:
        -------
            * rgb value [red, green, blue] values: 0.0 - 1.0
        """
        angle %= 360.0
        #
        if not brightness:
            brightness = 100.0
        #
        drgb = [(self.lowest + self.__R[int(angle) % len(self.__R)]
                             * (self.highest - self.lowest)) / 100.0,
                (self.lowest + self.__G[int(angle) % len(self.__G)]
                             * (self.highest - self.lowest)) / 100.0,
                (self.lowest + self.__B[int(angle) % len(self.__B)]
                             * (self.highest - self.lowest)) / 100.0]
        for rgb_i in range(len(drgb)):
            drgb[rgb_i] = drgb[rgb_i] * brightness / 10000.0
        return drgb
    def complement(self, angle):
        """
        Get two complementer color
        """
        return (self.rgb(angle),
                self.rgb(angle + 180.0))
    def triad(self, angle):
        """
        Get three colors
        """
        return (self.rgb(angle),
                self.rgb(angle + 180.0 - self.triadangle / 2.0),
                self.rgb(angle + 180.0 + self.triadangle / 2.0))
    def tetrad(self, angle):
        """
        Get four matching colors
        """
        return (self.rgb(angle),
                self.rgb(angle + self.tetradangle),
                self.rgb(angle + 180.0 - self.tetradangle),
                self.rgb(angle + 180.0))
    def analog(self, angle):
        """
        Get three analog colors
        """
        return (self.rgb(angle),
                self.rgb(angle + self.analogicangle),
                self.rgb(angle - self.analogicangle))
    def series(self, number, totalnumber, startangle, endangle, brightness = None):
        return self.rgb(startangle + number*(endangle - startangle)/float(totalnumber-1), brightness)
    def __setlowest(self, value):
        """
        Define the darkest rgb value : 0 - 100
        """
        self.__lowest = value % 100.0
        return None
    def __sethighest(self, value):
        """
        Define the lightest rgb value : 0 - 100
        """
        self.__highest = value % 100.0
        return None
    def __settriadangle(self, value):
        """
        Triadangle goes from 20 to 90 deg - default = 30 deg
        """
        self.__triadangle = value
        if self.__triadangle < 20.0:
            self.__triadangle = 20.0
        if self.__triadangle > 90.0:
            self.__triadangle = 90.0
        return None
    def __settetradangle(self, value):
        """
        Tetradangle goes from -90 to 90 deg, but not between -5 to 5 deg
        - default = 30 deg
        """
        self.__tetradangle = value
        if self.__tetradangle < -90.0:
            self.__tetradangle = -90.0
        if self.__tetradangle > 90.0:
            self.__tetradangle = 90.0
        if abs(self.__tetradangle) < 5.0:
            self.__tetradangle = 5.0 * self.__tetradangle / self.__tetradangle
        return None
    def __setanalogicangle(self, value):
        """
        Tetradangle goes from -90 to 90 deg, but not between -5 to 5 deg
        - default = 30 deg
        """
        self.__analogicangle = value
        if self.__analogicangle < 5.0:
            self.__analogicangle = 5.0
        if self.__analogicangle > 90.0:
            self.__analogicangle = 90.0
        return None
    def __getlowest(self):
        """
        Get the darkest rgb value
        """
        return self.__lowest
    def __gethighest(self):
        """
        Get the ligthest rgb value
        """
        return self.__highest
    def __gettriadangle(self):
        """
        Return the triad angle value
        """
        return self.__triadangle
    def __gettetradangle(self):
        """
        Return the tetrad angle value
        """
        return self.__tetradangle
    def __getanalogicangle(self):
        """
        Return the tetrad angle value
        """
        return self.__analogicangle
    lowest = property(__getlowest, __setlowest)
    highest = property(__gethighest, __sethighest)
    triadangle = property(__gettriadangle, __settriadangle)
    tetradangle = property(__gettetradangle, __settetradangle)
    analogicangle = property(__getanalogicangle)


