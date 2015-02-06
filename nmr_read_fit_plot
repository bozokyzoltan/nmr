#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created by Zoltan Bozoky on 2013.03.13.
Under GPL licence.

Module handling a Sparky NMR file without nmrglue or additional nmr proccessing 
or handling softwares

Purpose:
========
* One program that reads, visualizes multidimensional NMR data, finds peaks, 
  fit peak shapes, calculate volumes and intensities

"""
import struct
import random
import math
from pylab import plt
# Color handling module
import zcolor
# Fitting module
import GeneralFit as GF

################################################################################
# TODO Make the plot axes ratio to the right N-H-C ratios
################################################################################

def ppm2fid(ppm, spectral_data):
    """
    Convert ppm information into fid number
    """
    Frequency, MiddlePPM, SpectralWidth, NumberOfPoints = spectral_data
    return int((NumberOfPoints/2 - ((ppm-MiddlePPM) * Frequency * NumberOfPoints) / SpectralWidth) % NumberOfPoints)
###########################
def fid2ppm(fid, spectral_data):
    """
    Convert fid number into ppm information
    """
    Frequency, MiddlePPM, SpectralWidth, NumberOfPoints = spectral_data
    return MiddlePPM + (NumberOfPoints*SpectralWidth - 2*fid*SpectralWidth) / (2.0*Frequency*NumberOfPoints)
###########################
def distance(pos1, pos2):
    """
    Calculate Euclidean distance between two points
    """
    distance_value = 0.0
    for (p1, p2) in (pos1, pos2):
        distance_value += (p1 - p2)**2
    return math.sqrt(distance_value)
###########################
def ceil(number):
    """
    Return the closest higher integer to the number
    """
    if number - int(number) != 0:
        number = int(number) + 1
    return int(number)
###########################
def Gauss(Peak_coordinate, Coordinate, Szoras):
    """ gaussian peak """
    return (1/(Szoras*math.sqrt(2.0*math.pi)))*math.exp(-1*(Peak_coordinate-Coordinate)**2/(2.0*Szoras**2))
def Lorentz(Peak_coordinate, Coordinate, Szoras):
    """ lorentzian peak """
    return 1/(math.pi*Szoras*(1+((Peak_coordinate-Coordinate)/float(Szoras))**2))
###########################
def parabolic(x, p):
    """
    Fit a parabolic to the tip of the peaks
    """
    c, b, a = p
    y = a*(x-b)**2 + c
    return y
###########################
def linewidth_fit(x, p):
    """
    Linewidth fit error function
    """
    lw, [height, peak_pos] = p

    value = height * math.exp(-1 * (peak_pos - x)**2 / (2.0 * lw**2))
    print 'param',lw, height, peak_pos, value
    return value

def linewidth_fit2(x, p):
    """
    Linewidth fit error function 2
    """
    height, lw, peak_pos = p

    value = height * math.exp(-1 * (peak_pos - x)**2 / (2.0 * lw**2)) #- 38237.4296875
    return value


################################################################################
################################################################################
################################################################################

############################
## Sparky file header class
############################

class SparkyFileHeader(object):
    """
    Everything about the sparky file header

    - the first 180 bytes in the sparky file
    """
    def __init__(self, headerinfo):
        """
        """
        self._extract_file_header_information(headerinfo)
        
        self._header_info = {}
        #
        return None
    ##########
    def _extract_file_header_information(self, header):
        """
        """
        infos = struct.unpack('>10s 4c 9s 26s 80s 3x l 40s 4x', header)
        self._header_info['Sparky ID'           ] = str(infos[0]).strip('\x00')
        self._header_info['Number of Axis'      ] = ord(infos[1]) #
        self._header_info['Number of Components'] = ord(infos[2]) # = 1 for real data
        self._header_info['Encoding'            ] = ord(infos[3])
        self._header_info['Version'             ] = ord(infos[4]) # = 2 for current format
        self._header_info['Owner'               ] = str(infos[5]).strip('\x00')
        self._header_info['Date'                ] = str(infos[6]).strip('\x00')
        self._header_info['Comment'             ] = str(infos[7]).strip('\x00')
        self._header_info['Seek Position'       ] = str(infos[8]).strip('\x00')
        self._header_info['Scratch'             ] = str(infos[9]).strip('\x00')
        return None
    ##########
    def _get_number_of_axis(self):
        return self._header_info['Number of Axis']
    ##########
    number_of_axis = property(_get_number_of_axis)
    ##########

################################################################################
################################################################################
################################################################################

############################
## Sparky axis class
############################

class SparkyFileAxis(object):
    """
    Everything what axis must know
    - 128 bytes for each axis
    """
    def __init__(self, axisinfo):
        """
        """
        self._extract_axis_information(axisinfo)
        return None
    ##########
    def _extract_axis_information(self, axisdata):
        """
        """
        infos = struct.unpack('>6s h 3I 6f 84s', axisdata)
        self._axis_info = {}
        self._axis_info['Nucleus'               ] = str(infos[0]).strip('\x00') # nucleus name (1H, 13C, 15N, 31P, ...
        self._axis_info['Spectral Shift'        ] = infos[1] # to left or right shift
        self._axis_info['Number of Points'      ] = infos[2] # # of active data points - integer number of data points along this axis
        self._axis_info['Size'                  ] = infos[3] # total size of axis
        self._axis_info['BlockSize'             ] = infos[4] # # of points per cache block - integer tile size along this axis
        self._axis_info['Spectrometer frequency'] = infos[5] # MHz - float spectrometer frequency for this nucleus (MHz)
        self._axis_info['Spectral width'        ] = infos[6] # Hz  - float spectral width
        self._axis_info['xmtr frequency'        ] = infos[7] # transmitter offset (ppm) - float center of data (ppm)
        self._axis_info['Zero order'            ] = infos[8] # phase corrections
        self._axis_info['First order'           ] = infos[9] # phase corrections
        self._axis_info['First pt scale'        ] = infos[10] # scaling for first point
        self._axis_info['Extended'              ] = str(infos[11]).strip('\x00') #
        #
        self._axis_info['Scale'] = []
        for fid in range(0, int(self._axis_info['Number of Points']) + 1, 1):
            self._axis_info['Scale'].append(fid2ppm(fid, self.frq_carrier_sw_np))
        return None
    ##########
    def _get_parameter(self, parameter_name):
        return self._axis_info[parameter_name]
    ##########
    def _get_blocksize(self):
        return self._get_parameter('BlockSize')
    def _get_nucleus(self):
        return self.nucleus_info[-1]
    def _get_nucleus_info(self):
        return self._get_parameter('Nucleus')
    def _get_numberofpoints(self):
        return self._get_parameter('Number of Points')
    def _get_scale(self):
        return self._get_parameter('Scale')
    def _get_size(self):
        return self._get_parameter('Size')
    def _get_spectrometer_frequency(self):
        return self._get_parameter('Spectrometer frequency')
    def _get_spectral_width(self):
        return self._get_parameter('Spectral width')
    def _get_xmtr_frequency(self):
        return self._get_parameter('xmtr frequency')
    def _get_infos(self):
        return (self.spectrometer_frequency, self.xmtr_frequency,
                self.spectral_width, self.number_of_points)
    def ppm2index(self, ppm):
        index = 0
        while (index < self.number_of_points) and (self.scale[index] > ppm):
            index += 1
        return index
    def index2ppm(self, index):
        return fid2ppm(index, self.frq_carrier_sw_np)
    ##########
    blocksize = property(_get_blocksize)
    nucleus = property(_get_nucleus)
    nucleus_info = property(_get_nucleus_info)
    number_of_points = property(_get_numberofpoints)
    scale = property(_get_scale)
    size = property(_get_size)
    spectral_width = property(_get_spectral_width)
    spectrometer_frequency = property(_get_spectrometer_frequency)
    xmtr_frequency = property(_get_xmtr_frequency)
    frq_carrier_sw_np = property(_get_infos)
    ##########

################################################################################
################################################################################
################################################################################

############################
## Sparky spectral object
############################

class SparkySpectrum(object):
    """
    """
    def __init__(self, spectralinfo, blocksize_size_for_each_axis, log = True):
        """
        Parameters:
        ===========
            * spectralinfo = sparky file content with the spectral information
            * blocksize_size_for_each_axis = blocksize,size pairs for the axis
            * log = print out file processing information on the fly
        """
        self._log = log
        self._number_of_dimensions = len(blocksize_size_for_each_axis)
        self._d1 = None
        self._d2 = None
        self._d3 = None
        self._Spectrum = []
        self._noise_level = None
        #
        if self._log:
            print 'File read has started:',
        eval('self._extract_' + str(self.number_of_dimensions) +
             'D_data(spectralinfo, blocksize_size_for_each_axis)')
        if self._log:
            print '100% file read is done.'
        return None
    ##########
    def _extract_1D_data(self, Filecontent, Blocksize):
        """
        """
        self._Spectrum = list(struct.unpack
                             ('>'+'f'*(len(Filecontent)/4), Filecontent))
        return None
    ##########
    def _extract_2D_data(self, Filecontent, Blocksize):
        """
        """
        # First dimensional data
        FirstDimensionBlockSize = Blocksize[0]['BlockSize']
        FirstDimensionSpectralSize = Blocksize[0]['Size']
        # Second dimensional data
        SecondDimensionBlockSize = Blocksize[1]['BlockSize']
        SecondDimensionSpectralSize = Blocksize[1]['Size']
        #
        Blocksize = FirstDimensionBlockSize * SecondDimensionBlockSize
        # The number of blocks needed for a spectral size is
        # not necessary an integer number
        NumberOfBlocksInSecondDimension = (
            ceil(SecondDimensionSpectralSize / float(SecondDimensionBlockSize)))
        #---------------------------------
        # Rearrange the data from a list to an array
        for i_FirstDimension in range(FirstDimensionSpectralSize):
            # Print out info to follow the processing
            if self._log and i_FirstDimension % 50 == 0:
                print '{0:3.2f}%'.format(100.0 * i_FirstDimension
                                         / FirstDimensionSpectralSize),
            #---------------------------------
            BlockNumber = (i_FirstDimension / FirstDimensionBlockSize
                           * NumberOfBlocksInSecondDimension)
            PositionWithinBlock = (i_FirstDimension
                                   % FirstDimensionBlockSize
                                   * SecondDimensionBlockSize)
            # Concatenate the block portions in a list
            SpectralInfo1D = []
            #---------------------------------
            # Go through all second dimension protion to get a line
            for i_SecondDimension in range(NumberOfBlocksInSecondDimension):
                # If this is the last Block in line then the dimension is
                # not necessary the blocksize
                if i_SecondDimension < NumberOfBlocksInSecondDimension:
                    SecondDimension = SecondDimensionBlockSize
                else:
                    SecondDimension = (SecondDimensionSpectralSize
                                       % SecondDimensionBlockSize)
                #---------------------------------
                # The actual position within the block; 1 float number = 4 bytes
                pos = (4 * (Blocksize * (BlockNumber + i_SecondDimension)
                       + PositionWithinBlock))
                #---------------------------------
                # Unpack the data. Note that the coding is big endian ">"
                SpectralInfo1D += list(struct.unpack('>'+'f'*SecondDimension,
                                  Filecontent[pos : pos + 4 * SecondDimension]))
            #---------------------------------
            # Add a line into the spectrum
            self._Spectrum.append(SpectralInfo1D)
        return None
    ##########
    def _extract_3D_data(self, Filecontent, Blocksize):
        """
        """
        # Third dimensional data
        ThirdDimensionBlockSize         = Blocksize[0]['BlockSize']
        ThirdDimensionSpectralSize      = Blocksize[0]['Size']
        # Second dimensional data
        SecondDimensionBlockSize        = Blocksize[1]['BlockSize']
        SecondDimensionSpectralSize     = Blocksize[1]['Size']
        # First dimensional data
        FirstDimensionBlockSize         = Blocksize[2]['BlockSize']
        FirstDimensionSpectralSize      = Blocksize[2]['Size']
        #
        Blocksize = (FirstDimensionBlockSize
                    * SecondDimensionBlockSize
                    * ThirdDimensionBlockSize)
        #---------------------------------
        # The number of blocks needed for a spectral size is not necessary
        # an integer number
        NumberOfBlocksInFirstDimension  = ceil(FirstDimensionSpectralSize
                                              / float(FirstDimensionBlockSize ))
        NumberOfBlocksInSecondDimension = ceil(SecondDimensionSpectralSize
                                              / float(SecondDimensionBlockSize))
        #---------------------------------
        # Rearrange the data from a list to an 3D array
        for i_ThirdDimension in range(ThirdDimensionSpectralSize):
            # Print out log information
            if self._log and i_ThirdDimension % 10 == 0:
                print '{0:3.2f}%'.format(100.0 * i_ThirdDimension
                                        / ThirdDimensionSpectralSize),
            #---------------------------------
            BlockNumberDim3 = ((i_ThirdDimension
                                / ThirdDimensionBlockSize)
                                * NumberOfBlocksInSecondDimension
                                * NumberOfBlocksInFirstDimension)
            PositionWithinBlockDim3 = ((i_ThirdDimension
                                        % ThirdDimensionBlockSize)
                                        * SecondDimensionBlockSize
                                        * FirstDimensionBlockSize)
            #---------------------------------
            # Collect data of 2D in a variable
            SpectralInfo2D = []
            # Go through each block in 2D
            #for i_SecondDimension in range(SecondDimensionBlockSize * NumberOfBlocksInSecondDimension):
            for i_SecondDimension in range(SecondDimensionSpectralSize):
                #
                BlockNumberDim2 = (BlockNumberDim3
                                   + (i_SecondDimension / SecondDimensionBlockSize)
                                   * NumberOfBlocksInFirstDimension)
                PositionWithinBlockDim2 = (PositionWithinBlockDim3
                                          + (i_SecondDimension % SecondDimensionBlockSize)
                                          * FirstDimensionBlockSize)
                #---------------------------------
                # Collect data of 1D in a variable
                SpectralInfo1D = []
                # Go through each block in 1D
                for i_FirstDimension in range(NumberOfBlocksInFirstDimension):
                    # The last block size might be smaller than a blocksize
                    if i_FirstDimension < NumberOfBlocksInFirstDimension-1:
                        FirstDimension = FirstDimensionBlockSize
                    else:
                        FirstDimension = FirstDimensionSpectralSize % FirstDimensionBlockSize
                    #---------------------------------
                    # Position within block; 1 float number = 4 bytes
                    pos = 4 * (Blocksize * (BlockNumberDim2 + i_FirstDimension) + PositionWithinBlockDim2)
                    #---------------------------------
                    # Unpack the data. NOTE: big endian data storage ">"
                    SpectralInfo1D += list(struct.unpack('>'+'f'*FirstDimension,Filecontent[pos: pos + 4*FirstDimension]))
                #---------------------------------
                # Put each 1D slice into the 2D
                SpectralInfo2D.append(SpectralInfo1D)
            #---------------------------------
            # Store a 2D slice into the final array
            self._Spectrum.append(SpectralInfo2D)
        return None
    ##########
    def intensity(self, position):
        """
        Return an intensity value corresponds to the position
        """
        data_intensity = 0.0
        if self.number_of_dimensions == 1:
            data_intensity = (self._Spectrum[position[0] % self.dim1])
        if self.number_of_dimensions == 2:
            data_intensity = (self._Spectrum[position[1] % self.dim1]
                                            [position[0] % self.dim2])
        if self.number_of_dimensions == 3:
            data_intensity = (self._Spectrum[position[2] % self.dim1]
                                            [position[1] % self.dim2]
                                            [position[0] % self.dim3])
        return data_intensity
    ##########
    def calculate_noise_level(self, number_of_points = 10000):
        """
        """
        noise = 0.0
        # calculate the average level on a small subset of data
        average = 0.0
        pre_set_size = 100
        for i in range(pre_set_size):
            if self.number_of_dimensions == 1:
                average += self.intensity([random.randint(0, self.dim1 - 1)])
            if self.number_of_dimensions == 2:
                average += self.intensity([random.randint(0, self.dim1 - 1),
                                           random.randint(0, self.dim2 - 1)])
            if self.number_of_dimensions == 3:
                average += self.intensity([random.randint(0, self.dim1 - 1),
                                           random.randint(0, self.dim2 - 1),
                                           random.randint(0, self.dim3 - 1)])
        average /= float(pre_set_size)
        # Calculate the actual noise level
        numberofdata = 0
        sumofdata = 0.0
        highestvalue = 0.0
        i = 0
        while (i <= number_of_points*2) and (numberofdata <= number_of_points):
            if self.number_of_dimensions == 1:
                value = abs(self.intensity([random.randint(0, self.dim1 - 1)]))
            if self.number_of_dimensions == 2:
                value = abs(self.intensity([random.randint(0, self.dim1 - 1),
                                            random.randint(0, self.dim2 - 1)]))
            if self.number_of_dimensions == 3:
                value = abs(self.intensity([random.randint(0, self.dim1 - 1),
                                            random.randint(0, self.dim2 - 1),
                                            random.randint(0, self.dim3 - 1)]))
            # Only count a value if that is not far from the average
            # (= not a peak)
            if value < average * 5:
                numberofdata  += 1
                sumofdata += value
                average = sumofdata / float(numberofdata)
                if value > highestvalue:
                    highestvalue = value
            i += 1
        # Cut back from the highest to have a bit of noise
        noise = highestvalue / 1.2
        # Return the value as well
        return noise
    ##########
    def slice_1d(self, minmax, orderXY):
        """
        Return a 1D sub spectrum
        """
        highestvalue = 0.0
        lowestvalue = 0.0
        spectrum = []
        #---------------------------------
        # 1D
        if self.number_of_dimensions == 1:
            for x in range(min(minmax['X']), max(minmax['X']), 1):
                value = self.intensity([x])
                spectrum.append(value)
                if highestvalue < value:
                    highestvalue = value
                if lowestvalue > value:
                    lowestvalue = value
        #---------------------------------
        # 2D
        if self.number_of_dimensions == 2:
            y = min(minmax['Y'])
            for x in range(min(minmax['X']), max(minmax['X']), 1):
                if orderXY[0] == '0':
                    value = self.intensity([y, x])
                else:
                    value = self.intensity([x, y])
                spectrum.append(value)
                if highestvalue < value:
                    highestvalue = value
                if lowestvalue > value:
                    lowestvalue = value
        #---------------------------------
        # 3D
        if self.number_of_dimensions == 3:
            y = min(minmax['Y'])
            z = min(minmax['Z'])
            for x in range(min(minmax['X']), max(minmax['X']), 1):
                if orderXY[0:2] == '02':
                    value = self.intensity([y, z, x])
                elif orderXY[0:2] == '01':
                    value = self.intensity([z, y, x])
                elif orderXY[0:2] == '20':
                    value = self.intensity([x, z, y])
                elif orderXY[0:2] == '21':
                    value = self.intensity([x, y, z])
                elif orderXY[0:2] == '10':
                    value = self.intensity([z, x, y])
                elif orderXY[0:2] == '12':
                    value = self.intensity([y, x, z])
                else:
                    value = 0.0
                spectrum.append(value)
                if highestvalue < value:
                    highestvalue = value
                if lowestvalue > value:
                    lowestvalue = value

        return highestvalue, lowestvalue, spectrum
    ##########
    def slice_2d(self, minmax, orderXY):
        """
        Return a 2D sub spectrum
        """
        highestvalue = 0.0
        lowestvalue = 0.0
        spectrum = []
        #---------------------------------
        # 2D
        if self.number_of_dimensions == 2:
            for y in range(min(minmax['Y']), max(minmax['Y']), 1):
                fid = []
                for x in range(min(minmax['X']), max(minmax['X']), 1):
                    if orderXY[0] == '0':
                        value = self.intensity([y, x])
                    else:
                        value = self.intensity([x, y])
                    fid.append(value)
                    if highestvalue < value:
                        highestvalue = value
                    if lowestvalue > value:
                        lowestvalue = value
                spectrum.append(fid)
        # 3D
        if self.number_of_dimensions == 3:
            z = min(minmax['Z'])
            for y in range(min(minmax['Y']), max(minmax['Y']), 1):
                fid = []
                for x in range(min(minmax['X']), max(minmax['X']), 1):
                    if orderXY[0:2] == '02':
                        value = self.intensity([y, z, x])
                    elif orderXY[0:2] == '01':
                        value = self.intensity([z, y, x])
                    elif orderXY[0:2] == '20':
                        value = self.intensity([x, z, y])
                    elif orderXY[0:2] == '21':
                        value = self.intensity([x, y, z])
                    elif orderXY[0:2] == '10':
                        value = self.intensity([z, x, y])
                    elif orderXY[0:2] == '12':
                        value = self.intensity([y, x, z])
                    else:
                        value = 0.0
                    fid.append(value)
                    if highestvalue < value:
                        highestvalue = value
                    if lowestvalue > value:
                        lowestvalue = value
                spectrum.append(fid)
        return highestvalue, lowestvalue, spectrum
    ##########
    def slice_3d(self, minmax, orderXY):
        """
        Return a 3D sub spectrum
        """
        highestvalue = 0.0
        lowestvalue = 0.0
        spectrum = []
        #---------------------------------
        # 3D
        if self.number_of_dimensions == 3:
            for z in range(min(minmax['Z']), max(minmax['Z']), 1):
                fid2d = []
                for y in range(min(minmax['Y']), max(minmax['Y']), 1):
                    fid = []
                    for x in range(min(minmax['X']), max(minmax['X']), 1):
                        if orderXY[0:2] == '02':
                            value = self.intensity([y, z, x])
                        elif orderXY[0:2] == '01':
                            value = self.intensity([z, y, x])
                        elif orderXY[0:2] == '20':
                            value = self.intensity([x, z, y])
                        elif orderXY[0:2] == '21':
                            value = self.intensity([x, y, z])
                        elif orderXY[0:2] == '10':
                            value = self.intensity([z, x, y])
                        elif orderXY[0:2] == '12':
                            value = self.intensity([y, x, z])
                        else:
                            value = 0.0
                        fid.append(value)
                        if highestvalue < value:
                            highestvalue = value
                        if lowestvalue > value:
                            lowestvalue = value
                fid2d.append(fid)
            spectrum.append(fid2d)
        return highestvalue, lowestvalue, spectrum
    ##########
    def _get_dimension1(self):
        if not self._d1 and self.number_of_dimensions >= 1:
            self._d1 = len(self._Spectrum)
        return self._d1
    def _get_dimension2(self):
        if not self._d2 and self.number_of_dimensions >= 2:
            self._d2 = len(self._Spectrum[0])
        return self._d2
    def _get_dimension3(self):
        if not self._d3 and self.number_of_dimensions >= 3:
            self._d3 = len(self._Spectrum[0][0])
        return self._d3
    def _get_dimensions(self):
        return self._number_of_dimensions
    def _get_noiselevel(self):
        if not self._noise_level:
            self._noise_level = self.calculate_noise_level()
        return self._noise_level
    ##########
    dim1 = property(_get_dimension1)
    dim2 = property(_get_dimension2)
    dim3 = property(_get_dimension3)
    number_of_dimensions = property(_get_dimensions)
    noise_level = property(_get_noiselevel)
    ##########

################################################################################
################################################################################
################################################################################

class ChemicalShift(object):
    """
    Storing one chemical shift
    """
    def __init__(self):
        self._value = None
        self._original_value = None
        return None
    ##########
    def shift(self, value):
        self.chemical_shift = self.chemical_shift + value
        return None
    ##########
    def _get_cs(self):
        if not self._value:
            value = 0.0
        else:
            value = self._value
        return value
    def _set_cs(self, newvalue):
        self._value = newvalue
        if not self._original_value:
            self._original_value = newvalue
        return None
    def _get_original_cs(self):
        if not self._original_value:
            value = 0.0
        else:
            value = self._original_value
        return value
    ##########
    chemical_shift = property(_get_cs, _set_cs)
    original_chemical_shift = property(_get_original_cs)
    ##########

class Peak(object):
    """
    Storing all chemical shift for a peak:
    Parameters:
    ===========
        * adjusted
        * info
        * name
        * nucleus
        * chemical_shift
        * original_chemical_shift
        * add_chemical_shift
        * shift
    """
    def __init__(self):
        self.CSs = {}
        self._adjusted = False
        self._intensity = None
        return None
    def add_chemical_shift(self, nucleus, chemical_shift):
        if not nucleus in self.CSs:
            self.CSs[nucleus] = ChemicalShift()
        self.CSs[nucleus].chemical_shift = float(chemical_shift)
        return None
    def chemical_shift(self, nucleus):
        if nucleus in self.CSs:
            value = self.CSs[nucleus].chemical_shift
        else:
            value = 0.0
        return value
    def original_chemical_shift(self, nucleus):
        if nucleus in self.CSs:
            value = self.CSs[nucleus].original_chemical_shift
        else:
            value = 0.0
        return value
    def shift(self, nucleus, value):
        if nucleus in self.CSs:
            self.CSs[nucleus].shift(value)
        return None
    def set_peak_info(self, line, peaknameinfo):
        colomns = line.split()
        self._info = colomns[0]
        spins = self.info.split('-')
        self._peakname = eval('spins[0]' + peaknameinfo)
        for i,spin in enumerate(spins):
            self.add_chemical_shift(spin[-1], colomns[i+1])
        return None
    ##########
    def _get_adjusted(self):
        return self._adjusted
    def _set_adjusted( self, value):
        self._adjusted = value
        return None
    def _get_info(self):
        return self._info
    def _set_info(self, value):
        self._info = value
        return None
    def _get_intensity(self):
        if not self._intensity:
            value = 0.0
        else:
            value = self._intensity
        return value
    def _set_intensity(self, value):
        self._intensity = value
        return None
    def _get_name(self):
        return self._peakname
    def _set_name(self, value):
        self._peakname = value
        return None
    def _get_nucleuses(self):
        return self.CSs.keys()
    ##########
    adjusted = property(_get_adjusted, _set_adjusted)
    info = property(_get_info, _set_info)
    intensity = property(_get_intensity, _set_intensity)
    name = property(_get_name, _set_name)
    nucleus = property(_get_nucleuses)
    ##########


class Peaklist(object):
    """
    Everything about peaklists
    """
    def __init__(self):
        self._peaks = {}
        self._number_of_peaks = -1
        return None
    # public procedures
    def read_peaklist(self, filename, info):
        """
        """
        self.filename = filename
        try:
            peaklist_file = open(filename, 'r')
        except IOError:
            print 'Error opening ' + filename + '!!! Please check it!'
            exit()
        lines = peaklist_file.readlines()
        peaklist_file.close()
        num = 0
        for line in lines:
            if (not 'Assignment' in line) and (len(line) > 5):
                line.strip()
                self._peaks[num] = Peak()
                self._peaks[num].set_peak_info(line, info)
                num += 1
        self._number_of_peaks = num - 1
        return None
    def print_peaklist(self, filename = None):
        """
        """
        if filename:
            fil = open(filename,'w')
        for i in range(self.number_of_peaks):
            nucleus = self._peaks[i].nucleus
            nucleus.reverse()
            line = self._peaks[i].name
            for j, nuc in enumerate(nucleus):
                if j == 0:
                    line = ''.join([line, '_', nuc])
                else:
                    line = ''.join([line, '-', nuc])
            for nuc in nucleus:
                line = ' '.join([line, str(self._peaks[i].chemical_shift(nuc))])
            if filename:
                line = ''.join([line, '\n'])
                fil.write(line)
            else:
                print line
        if filename:
            fil.close()
        return None
    def add_peak(self, peak_info):
        """
        Needs a hash line {'N':129.3,'H':8.5,'C':178.2}
        """
        number = self.number_of_peaks
        self._peaks[number] = Peak()
        for info in peak_info:
            self._peaks[number].add_chemical_shift(info, peak_info[info])
        self._peaks[number].info = str(number + 1)
        self._peaks[number].name = str(number + 1)
        self._number_of_peaks += 1
        return None
    def adjust(self, num):
        self._peaks[num].adjusted = True
        return None
    def adjusted(self, num):
        return self._peaks[num].adjusted
    def add_cs(self, num, nucleus, value):
        self._peaks[num].add_chemical_shift(nucleus, value)
        return None
    def cs(self, num, nucleus):
        return self._peaks[num].chemical_shift(nucleus)
    def add_intensity(self, num, value):
        self._peaks[num].intensity = value
        return None
    def intensity(self, num):
        return self._peaks[num].intensity
    def original_cs(self, num, nucleus):
        return self._peaks[num].original_chemical_shift(nucleus)
    def name(self, num):
        return self._peaks[num].name
    def nucleus(self, num):
        return self._peaks[num].nucleus
    def info(self, num):
        return self._peaks[num].info
    def shift(self, num, nucleus, value):
        self._peaks[num].shift(nucleus, value)
        return None
    # private procedures
    ##########
    def _get_number_of_peaks(self):
        return self._number_of_peaks + 1
    ##########
    number_of_peaks = property(_get_number_of_peaks)
    ##########


################################################################################
################################################################################
################################################################################

class Sparky_plot(object):
    """
    Make a plot of 1d or 2d spectrum
    """
    _plot_already = None
    def __init__(self):
        self._noiselevel = 0.0
        self._number_of_contours = 25
        self._factor = 1.1
        self._colors = []
        self._first_plot_on_figure = False
        self._plot_negatives = True
        #
        self.mycolor = zcolor.MyColor()
        self.colors = [self.mycolor.series(i, self.number_of_contours, 0, 330, 100.0) for i in range(self.number_of_contours)]
        #
        if not self._plot_already:
            self._plot_already = 1
            self.newfigure()
        return None
    ##########
    def newfigure(self):
        plt.figure()
        self._plot_already = 1
        self._first_plot_on_figure = True
        return None
    ##########
    def plot_1d(self, x_axis_scale, spectral_data, axes_label, color = None):
        """
        Plot a 1D slice
        """
        if self._first_plot_on_figure:
            # plot zero
            plt.plot([x_axis_scale[0],x_axis_scale[-1]],[0.0,0.0],'k-')
            # plot noise level
            plt.plot([x_axis_scale[0],x_axis_scale[-1]],[self.noise_level,self.noise_level],'k--')
        #----------------
        # color selection
        if not color:
            try:
                plotcolor = self.colors[0]
            except IndexError:
                plotcolor = 'k'
        else:
            plotcolor = color
        # plot the data
        plt.plot(x_axis_scale, spectral_data, color = plotcolor)
        # set the x axis limit
        plt.xlim(x_axis_scale[0],x_axis_scale[-1])
        #
        if self._first_plot_on_figure:
            plt.xlabel(axes_label[0] + ' (ppm)', size = 15)
            plt.ylabel('Intensity', size = 15)
        return None
    ##########
    def plot_2d(self, x_axis_scale, y_axis_scale, spectral_data, axes_label, color = None):
        """
        Plot a 2D spectrum
        """
        # Colors
        if not color:
            if len(self.colors) < self.number_of_contours:
                plotcolors = []
                for i in range(self.number_of_contours):
                    plotcolors.append([0.0, 0.0, 0.0])
            else:
                plotcolors = self.colors
        else:
            plotcolors = color
        # Contour levels
        contourlevels = [self.noise_level * self.factor**i for i in range(self.number_of_contours)]

        # plot positive contours
        plt.contour(x_axis_scale, y_axis_scale, spectral_data, contourlevels, colors = plotcolors)
        if self.plot_negatives:
            # plot negatives if needed!
            plt.contour(x_axis_scale, y_axis_scale, spectral_data, [-1*i for i in contourlevels], colors = [[0.0,0.0,0.0] for i in range(self.number_of_contours)])

        if self._first_plot_on_figure:
            # Invert the axis direction
            plt.gca().invert_xaxis()
            plt.gca().invert_yaxis()
            # Put label on axes
            plt.xlabel(axes_label[0] + ' (ppm)', size = 15)
            plt.ylabel(axes_label[1] + ' (ppm)', size = 15)
        return None
    ##########
    def show(self, filename = None):
        """
        Show or save the figure depending on whether filename is provided
        """
        if not filename:
            plt.show()
        else:
            plt.savefig(filename)
        return None
    ##########
    def plot_peaklist_2d(self, peaklist, orderXY):
        """
        """
        print 'Peaks on the plot:'
        print '  #   name     cs1      cs2        intensity    adjusted'
        print '--------------------------------------------------------'
        for number, peak in enumerate(peaklist):
            #
            info = peak
            loc_x = peaklist[peak][orderXY[0]]
            loc_y = peaklist[peak][orderXY[1]]
            adj = peaklist[peak]['Adjusted']
            intensity = peaklist[peak]['Intensity']
            #
            print '{0:3d}. {1:>5s}   {2:7.3f}  {3:7.3f} {4:14.3f}     '.format(number + 1, peak, loc_y, loc_x, intensity),
            if adj:
                print 'true'
                labelcolor = 'black'
            else:
                print 'false'
                labelcolor = 'red'
            #
            dx = 0.0
            dy = 0.2
            plt.gca().annotate(info,
                               xy = (loc_x, loc_y),
                               color = labelcolor,
                               xytext = (loc_x - dx,loc_y - dy),
                               arrowprops = dict(arrowstyle = "-|>",
                                                 connectionstyle = "arc3",
                                                 facecolor = labelcolor))
        print '--------------------------------------------------------'
        return None
    ##########
    def set_factor(self, highestvalue):
        #
        self.factor = math.exp(math.log(highestvalue /float(self.noise_level)) * 1.0/(float(self.number_of_contours)))
        return self.factor
    ##########
    def _set_colors(self, levels):
        self._colors = levels
        return None
    def _get_colors(self):
        return self._colors
    def _set_factor(self, level):
        self._factor = level
        return None
    def _get_factor(self):
        return self._factor
    def _set_noiselevel(self, level):
        self._noiselevel = level
        return None
    def _get_noiselevel(self):
        return self._noiselevel
    def _set_number_of_contours(self, level):
        self._number_of_contours = level
        return None
    def _get_number_of_contours(self):
        return self._number_of_contours
    def _set_plot_negatives(self, level):
        self._plot_negatives = level
        return None
    def _get_plot_negatives(self):
        return self._plot_negatives
    ##########
    colors = property(_get_colors, _set_colors)
    factor = property(_get_factor, _set_factor)
    noise_level = property(_get_noiselevel, _set_noiselevel)
    number_of_contours = property(_get_number_of_contours, _set_number_of_contours)
    plot_negatives = property(_get_plot_negatives, _set_plot_negatives)
    ##########


################################################################################
################################################################################
################################################################################


class ZB_spectrum(object):
    """
    """
    def __init__(self, filename):
        """
        """
        self._peaklist = Peaklist()
        #
        try:
            filehandler = open(filename, 'rb')
        except IOError:
            print ('ERROR!!!\nPlease check the ' + filename + ' location, '
                   'because an error happened during the file open...\n')
            exit()
        #---------------------------------
        # Read the file header information
        self.header = SparkyFileHeader(filehandler.read(180))
        #---------------------------------
        # Read the axes information
        self.axis = {}
        self.axis_order = ''
        blocksize_info = []
        for i in range(self.header.number_of_axis):
            axis = SparkyFileAxis(filehandler.read(128))
            self.axis_order += axis.nucleus
            self.axis[self.axis_order[-1]] = axis
            blocksize_info.append({'BlockSize':axis.blocksize, 'Size':axis.size})
        #---------------------------------
        # Read the spectral information
        self.spectrum = SparkySpectrum(filehandler.read(), blocksize_info)
        #---------------------------------
        filehandler.close()
        #
        return None
    ##########
    def _get_limits(self, limit, nucleus):
        if limit[nucleus] == []:
            index_min = 0
            index_max = self.axis[nucleus].number_of_points - 1
        else:
            index_min = self.axis[nucleus].ppm2index(max(limit[nucleus]))
            index_max = self.axis[nucleus].ppm2index(min(limit[nucleus]))
        return index_min, index_max + 1
    ##########
    def plot1d(self, limits, orderXY):
        """
        Parameters:
        ===========
            * limits: a hash with the limits in ppm
            * orderxY:

            example: plot1d({'H':[5.5,9.2],'N':[105,122]})
        """
        if not orderXY:
            orderXY = 'HN'
        # Dealing with the limits
        xy_limits = {}
        xy_limits['X'] = self._get_limits(limits, orderXY[0])
        if self.header.number_of_axis > 1:
            xy_limits['Y'] = self._get_limits(limits, orderXY[1])
        if self.header.number_of_axis > 2:
            xy_limits['Z'] = self._get_limits(limits, orderXY[2])
        # Dealing with the order
        axes_order = ''
        for i in range(len(orderXY)):
            axes_order += str(self.axis_order.index(orderXY[i]))
        #
        highest, lowest, spectral_data = self.spectrum.slice_1d(xy_limits, axes_order)

        scale = self.axis[orderXY[0]].scale[xy_limits['X'][0] : xy_limits['X'][1]]

        self.figure = Sparky_plot()

        self.figure.noise_level = self.spectrum.noise_level

        self.figure.plot_1d(scale, spectral_data, self.axis[orderXY[0]].nucleus_info, 'b')

        print '#############################################'
        print '###   P L O T   #   P A R A M E T E R S   ###'
        print '#############################################'
        print 'Noise level   =', self.figure.noise_level
        print 'Highest value =', highest
        print 'Lowest value  =', lowest
        print '#############################################'

        return None
    ##########
    def plot(self, limits, orderXY = None):
        """
        Parameters:
        ===========
            * limits: a hash with the limits in ppm
            * orderxY:

            example: plot1d({'H':[5.5,9.2],'N':[105,122]})
        """
        if not orderXY:
            orderXY = 'HN'
        # Dealing with the limits
        xy_limits = {}
        xy_limits['X'] = self._get_limits(limits, orderXY[0])
        xy_limits['Y'] = self._get_limits(limits, orderXY[1])
        if self.header.number_of_axis > 2:
            xy_limits['Z'] = self._get_limits(limits, orderXY[2])
        # Dealing with the order
        axes_order = ''
        for i in range(len(orderXY)):
            axes_order += str(self.axis_order.index(orderXY[i]))
        # Axis labels
        labels = []
        for o in orderXY:
            labels.append(self.axis[o].nucleus_info)
        #----------------
        highest, lowest, spectral_data = self.spectrum.slice_2d(xy_limits, axes_order)

        x_scale = self.axis[orderXY[0]].scale[xy_limits['X'][0] : xy_limits['X'][1]]
        y_scale = self.axis[orderXY[1]].scale[xy_limits['Y'][0] : xy_limits['Y'][1]]


        self.figure = Sparky_plot()

        self.figure.noise_level = self.spectrum.noise_level
        self.figure.set_factor(highest)

        print '#############################################'
        print '###   P L O T   #   P A R A M E T E R S   ###'
        print '#############################################'
        print 'Noise level   =', self.figure.noise_level
        print 'Factor        =', self.figure.factor
        print 'Highest value =', highest
        print 'Lowest value  =', lowest
        print '#############################################'


        self.figure.plot_2d(x_scale, y_scale, spectral_data, labels)
        # prepare peaklist
        peaklist = {}
        for i in range(self._peaklist.number_of_peaks):
            within = True
            for o in orderXY:
                cs = self._peaklist.cs(i, o)
                if limits[o] != []:
                    if (cs < min(limits[o])) or (max(limits[o]) < cs):
                        within = False
            if within:
                peaklist[self._peaklist.name(i)] = {}
                peaklist[self._peaklist.name(i)][orderXY[0]] = self._peaklist.cs(i, orderXY[0])
                peaklist[self._peaklist.name(i)][orderXY[1]] = self._peaklist.cs(i, orderXY[1])
                peaklist[self._peaklist.name(i)]['Adjusted'] = self._peaklist.adjusted(i)
                peaklist[self._peaklist.name(i)]['Intensity'] = self.spectrum.intensity([
                    self.axis[orderXY[0]].ppm2index(self._peaklist.cs(i, orderXY[0])),
                    self.axis[orderXY[1]].ppm2index(self._peaklist.cs(i, orderXY[1]))])
        # plot peaklist
        self.figure.plot_peaklist_2d(peaklist, orderXY)


        return None
    ###########################
    def show(self, filename = ''):
        """
        """
        self.figure.show(filename)
        return None
    ###########################
    def _extremes_finder(self, position, dimension, axis_order, find_max = True):
        """
        find positive and negative extremes on the spectrum
        Parameters:
        ===========
            * position = spectrum starting position for the peak finding,
                         order must be same as in the spectrum
            * dimension = find local maximum or minimum in 2D or 3D
            * find_max = maximum or minimum finding
        Return:
        =======
            * local extreme
        """
        checklist = [[-1, 0, 0],[+1, 0, 0], # x
                     [ 0,-1, 0],[ 0,+1, 0], # y
                     [-1,-1, 0],[+1,-1, 0], # xy
                     [-1,+1, 0],[+1,+1, 0], # xy
                     [ 0, 0,-1],[ 0, 0,+1], # z
                     [-1, 0,-1],[+1, 0,-1], # xz
                     [-1, 0,+1],[+1, 0,+1], # xz
                     [ 0,-1,-1],[ 0,-1,-1], # yz
                     [ 0,+1,+1],[ 0,+1,+1]] # yz
        #
        spectral_width = []
        for o in axis_order:
            spectral_width.append(eval('self.spectrum.dim' + str(int(o)+1)))
        #spectral_width = [self.spectrum.dim2, self.spectrum.dim1, self.spectrum.dim3]
        # If the dimension 2D, we find check the value in x,y otherwise in x,y,z
        if dimension == 2:
            checklist_size = 4
        else:
            checklist_size = len(checklist)
        # minimum or maximum finder
        finder_type = [['min','<'],['max','>']][find_max]
        # It goes till it finds a local maximum
        not_on_an_extreme_value = True
        while not_on_an_extreme_value:
            # check all values according to the checklist
            checked_values = []
            for check in checklist[0 : checklist_size]:
                checked_values.append(self.spectrum.intensity([pos + ch for (pos, ch) in zip(position[0 : dimension], check[0 : dimension])]))
            # if the position data is the most extreme, than we are done
            most_extreme_in_array = eval(eval('finder_type[0]') + '(checked_values)')
            if eval('self.spectrum.intensity(position)' + eval('finder_type[1]') + 'most_extreme_in_array'):
                not_on_an_extreme_value = False
            else:
                # modifiy the position to the highest
                checked_values_max_index = checked_values.index(most_extreme_in_array)
                for i in range(dimension):
                    position[i] += checklist[checked_values_max_index][i]
                    position[i] %= spectral_width[i]
        return position
    ###########################
    def _find_peak_1d(self, data, noiselevel):
        hits = []
        direction = True
        for i in range(len(data)-1):
            if data[i] > data[i+1] and data[i] > noiselevel and direction:
                hits.append(i)
                direction = False
            if data[i] < data[i+1]:
                direction = True

        if len(hits) > 0 and False:
            plt.figure()

            plt.plot(range(len(data)),data)
            plt.plot(hits,[50000 for j in range(len(hits))], 'k', marker= 'o', linestyle = '')

            plt.show()
        return hits
    ###########################
    def _find_peak_2d(self, data2d, noiselevel, order):
        hits = {}
        for i, data1d in enumerate(data2d):
            hit1d = self._find_peak_1d(data1d, noiselevel)
            for hit in hit1d:
                hits[' '.join(str(d) for d in self._extremes_finder([hit, i], 2, order))] = 0
        peaks = []
        for hit in hits:
            peaks.append(hit.split())
        return peaks
    ###########################
    def peak_finder_2d(self, orderXY = 'HN', times_noiselevel = 1.5):
        # Dealing with the order
        axes_order = ''
        for i in range(len(orderXY)):
            axes_order += str(self.axis_order.index(orderXY[i]))
        #
        xy = {}
        xy['X'] = [0, self.axis[orderXY[0]].number_of_points - 1]
        xy['Y'] = [0, self.axis[orderXY[1]].number_of_points - 1]
        #
        print 'Finding peaks...',
        peaklist = {}
        for i,peak in enumerate(self._find_peak_2d(self.spectrum.slice_2d(xy, axes_order)[-1],self.spectrum.noise_level*times_noiselevel, axes_order)):
            peak_info = {}
            for j, o in enumerate(orderXY):
                peak_info[o] = self.axis[o].index2ppm(float(peak[j]))
            self._peaklist.add_peak(peak_info)
            self._peaklist.adjust(self._peaklist.number_of_peaks - 1)
        print str(i + 1) + ' peaks found!'
        return peaklist
    ###########################
    def _one_peak(self, peaknumber, orderXY):
        if (0 <= peaknumber) and (peaknumber < self._peaklist.number_of_peaks):
            window = {'H':0.08,'N':0.5,'C':0.5}
            limit = {}
            for o in orderXY:
                limit[o] = [self._peaklist.cs(peaknumber, o) - window[o],self._peaklist.cs(peaknumber, o) + window[o]]
            self.plot(limit, orderXY)

            lim1d = {}
            o = orderXY[0]
            lim1d[o] = [self._peaklist.cs(peaknumber, o) - window[o], self._peaklist.cs(peaknumber, o) + window[o]]
            o = orderXY[1]
            lim1d[o] = [self._peaklist.cs(peaknumber, o)]
            self.plot1d(lim1d,orderXY)

            lim1dd = {}
            o = orderXY[1]
            lim1dd[o] = [self._peaklist.cs(peaknumber, o) - window[o], self._peaklist.cs(peaknumber, o) + window[o]]
            o = orderXY[0]
            lim1dd[o] = [self._peaklist.cs(peaknumber, o)]

            self.plot1d(lim1dd,orderXY[1]+orderXY[0])
        return None
    ###########################
    def _get_spectrum_around_peak(self, axis_order, position):
        """
        It returns 1d slices of the spectrum for peak fitting

        Parameters:
        ===========
            * axis_order = nucleus order in XYZ format
            * position = info as in spectrum

        Returns:
        ========
            * One dimensional slices: all, left, right, top
        """
        topwindow = 2
        permutation = {1 : {0: {'left':[        -1], 'right':[        +1]}},
                       2 : {0: {'left':[     0, -1], 'right':[     0, +1]},
                            1: {'left':[    -1,  0], 'right':[    +1,  0]}},
                       3 : {0: {'left':[ 0,  0, -1], 'right':[ 0,  0, +1]},
                            1: {'left':[ 0, -1,  0], 'right':[ 0, +1,  0]},
                            2: {'left':[-1,  0,  0], 'right':[+1,  0,  0]}}}
        slices = {}
        for dimension in axis_order:
            slices[dimension] = {}
            # Get the left and the right side of the peak separately
            for direction in ['left','right']:
                # Start from the original postion
                pos = []
                for p in position:
                    pos.append(p)
                # Collect the data
                tomb = []
                while self.spectrum.intensity(pos) >= self.spectrum.noise_level:
                    tomb.append(self.spectrum.intensity(pos))
                    for j in range(len(pos)):
                        pos[j] += permutation[len(position)][axis_order.index(dimension)][direction][j]
                # Store the data
                slices[dimension][direction] = tomb
            # extract the whole peak and just the top part
            slices[dimension]['all'] = []
            slices[dimension]['top'] = []
            for i in range(len(slices[dimension]['left'])):
                slices[dimension]['all'].append(slices[dimension]['left'][len(slices[dimension]['left']) - i - 1])
                if i <= topwindow:
                    slices[dimension]['top'].append(slices[dimension]['left'][topwindow - i])
            for i in range(1,len(slices[dimension]['right'])):
                slices[dimension]['all'].append(slices[dimension]['right'][i])
                if i <= topwindow:
                    slices[dimension]['top'].append(slices[dimension]['right'][i])
        return slices
    ###########################
    def _fit_one_slice_around_peak(self, spectrum, pos):
        """
        Fit a 1d array with a gaussian or lorentian curve
        """

        fit = GF.Fit_general(range(len(spectrum)),
                             spectrum,
                             [max(spectrum), len(spectrum)*0.7],
                             linewidth_fit2,
                             z = [pos for i in range(len(spectrum))],
                             Log = False,
                             NoErrors = 'NoErrors')
        print fit.Value, fit.Chi2/len(spectrum)
#        a,b = fit.GenerateCurve(0,len(spectrum))

#        plt.figure()
#        plt.plot(range(len(spectrum)), spectrum,'k',linestyle='',marker='o')
#        plt.plot(a,b)
#        plt.show()

        return fit.Value
    ###########################
    def peakfit(self, peaknumber):
        """
        """
        peakposition = []
        cs = []
        axisorder = ''
        for i in range(len(self.axis_order), 0, -1):
            ax = self.axis_order[i - 1]
            axisorder += ax
            cs.append(self._peaklist.cs(peaknumber, ax))
            peakposition.append(self.axis[ax].ppm2index(self._peaklist.cs(peaknumber, ax)))
        #
        slices = self._get_spectrum_around_peak(axisorder, peakposition)
        # fitting the tip of the peak
        intensity = []
        new_index = []
        linewidth = {}
        for i,ax in enumerate(axisorder):
            print ax
            linewidth[ax] = []
            spectrum = slices[ax]['top']
            fit = GF.Fit_general(range(len(spectrum)), spectrum, [max(spectrum), len(spectrum)//2, -1E+5], parabolic, Log = False, NoErrors = 'NoErrors')
            intensity.append(fit.Value[0])
            new_index.append(fit.Value[1] - len(spectrum)//2)

#            a,b = fit.GenerateCurve(0,len(spectrum)-1)
    
#            plt.figure()
#            plt.plot(range(len(spectrum)), spectrum,'k',linestyle='',marker='o')
#            plt.plot(a,b)
#            plt.show()
            # fit the sides of the peak
            for side in ['left','right','all']:
                spectrum = slices[ax][side]
                lw = self._fit_one_slice_around_peak(spectrum, spectrum.index(max(spectrum)) + new_index[-1])
                linewidth[ax].append(lw[1])
                

        print 'intensity:',sum(intensity) / len(intensity), intensity
        for i,ax in enumerate(axisorder):
            print 'position:',ax, self.axis[ax].index2ppm(peakposition[i] + new_index[i])
            print 'lw:',min(linewidth[ax]),self.axis[ax].index2ppm(min(linewidth[ax]))*self.axis[ax].spectrometer_frequency

        print axisorder
        print cs
        print peakposition
        print new_index


        exit()


        window = 3
        max_window_peak = 8
        order = {1:['0'], 2:['10','01'],3:['210','102','021']}
        axis = ['X','Y','Z']
        nucleuses = self._peaklist.nucleus(peaknumber)
        #
        index = {}
        for nuc in nucleuses:
            index[nuc] = self.axis[nuc].ppm2index(self._peaklist.cs(peaknumber, nuc))
        for orderXY in order[len(nucleuses)]:
            xy = {}
            xypeak_left = {}
            xypeak_right = {}
            for i, o in enumerate(orderXY):
                nuc = nucleuses[int(o)]
                ax = axis[i]
                xy[ax] = [index[nuc]]
                xypeak_left[ax] = [index[nuc]]
                xypeak_right[ax] = [index[nuc]]
            xy['X'] = [xy['X'][0] - window, xy['X'][0] + window + 1]
            xypeak_left['X'] = [xypeak_left['X'][0] - max_window_peak, xypeak_left['X'][0]]
            xypeak_right['X'] = [xypeak_right['X'][0], xypeak_right['X'][0] + max_window_peak + 1]

            rev_order = ''
            for o in orderXY:
                rev_order = ''.join([o, rev_order])

            # Fitting the tip of the peak with a parabolic
            spectrum = self.spectrum.slice_1d(xy, rev_order)[2]
            spectrum_peak_left = self.spectrum.slice_1d(xypeak_left, rev_order)[2]
            spectrum_peak_right = self.spectrum.slice_1d(xypeak_right, rev_order)[2]

            fit = GF.Fit_general(range(len(spectrum)), spectrum, [max(spectrum), window, -1E+5], parabolic, Log = False, NoErrors = 'NoErrors')

            xaxisnuc = nucleuses[int(orderXY[0])]

            index_diff = fit.Value[1] - window
            new_index = index[xaxisnuc] + index_diff
            ppm = self.axis[xaxisnuc].index2ppm(new_index)
            intensity = fit.Value[0]
            self._peaklist.add_cs(peaknumber, xaxisnuc, ppm)

            if xaxisnuc == 'H':
                self._peaklist.add_intensity(peaknumber, intensity)
            # Fitting the peak with a gaussian
##            fit_left = GF.Fit_general(range(len(spectrum_peak_left)),
##                                      spectrum_peak_left,
##                                      [intensity, 2],
##                                      linewidth_fit2,
##                                      z = [max_window_peak + index_diff for i in range(len(spectrum_peak_left))],
##                                      #z = [(max_window_peak + index_diff, min(spectrum_peak_left)) for i in range(len(spectrum_peak_left))],
##                                      Log = False,
##                                      NoErrors = 'NoErrors')
###            fit_left = GF.Fit_general(range(len(spectrum_peak_left)), spectrum_peak_left, [window_peak], linewidth_fit, z = [[intensity, window_peak + index_diff] for i in range(len(spectrum_peak_left))], Log = True, NoErrors = 'NoErrors')
###            fit_right = GF.Fit_general(range(len(spectrum_peak_right)), spectrum_peak_right, [window_peak], linewidth_fit, z = [[intensity, index_diff] for i in range(len(spectrum_peak_right))], Log = False, NoErrors = 'NoErrors')
##
##            print fit_left.Value
#            print fit_right.Value


##            a,b = fit_left.GenerateCurve(0,7)


            xy = {}
            for i, o in enumerate(orderXY):
                nuc = nucleuses[int(o)]
                ax = axis[i]
                xy[ax] = [index[nuc]]

            left = []
            y = xy['Y'][0]
            x = xy['X'][0]

            dd = self._get_spectrum_around_peak([y,x], 2)
#            print dd
            exit()

            while self.spectrum.intensity([y,x]) >= self.spectrum.noise_level:
                left.append(self.spectrum.intensity([y,x]))
                x = x - 1
            left.append(self.spectrum.intensity([y,x]))

            left.reverse()
            print len(left) + index_diff
            left_fit = GF.Fit_general(range(len(left)),
                                      left,
                                      [max(left), 1.0],
                                      linewidth_fit2,
                                      z = [len(left) - 1 + index_diff for i in range(len(left))],
                                      Log = True,
                                      NoErrors = 'NoErrors')

            e,f = left_fit.GenerateCurve(0,7)


            plt.figure()
##            plt.plot(range(len(spectrum_peak_left)), spectrum_peak_left,'k',marker = 'o',linestyle= '')
##            plt.plot(a,b)

            plt.plot(range(len(left)), left,'r',marker = 'o',linestyle= '')
            plt.plot(e,f,'r--')

            plt.show()
            exit()
        return None
    ###########################
    def read_peaklist(self, peaklist_filename, info = '[:-2]'):
        self._peaklist.read_peaklist(peaklist_filename, info)
        return None
    ###########################
    def print_peaklist(self):
        self._peaklist.print_peaklist()
        return None
    ###########################
    def save_peaklist(self, peaklist_filename):
        self._peaklist.print_peaklist(peaklist_filename)
        return None
    ###########################
    def _get_noiselevel(self):
        return self.spectrum.noise_level
    ###########################
    noise_level = property(_get_noiselevel)
    ###########################






################################################################################
################################################################################
################################################################################


class SparkyFile(object):
    """
    """
    ###########################
    Plotting_parameters = []
    ###########################
    def __init__(self, filename, log = True):
        """
        Parameters:
            * filename = A sparky file with path information
            * log      = True to print out information during processing
        """
        # Information on dimensionarity of the measured data
        self._FileHeader_ = {}
        # Information on measured axes
        self.myaxis = {}
        #
        self._AxisOrder_ = []
        # Spectrum data
        self._Spectrum_ = []
        # Peaklist information
        self._Peaklist_ = {}
        # Store the peaklist keys in order of the read in
        self._Peaklistkeysorder_ = []
        #
        self._Peaklistchemicalshiftorder_ = []
        #
        self._PeaklistDoNotCare_ = []
        #
        self.Noiselevel = None

        #---------------------------------
        self.log = log
        # Open the sparky file
        try:
            filehandler = open(filename, 'rb')
        except IOError:
            print ('ERROR!!!\nPlease check the ' + filename + ' location, '
                   'because an error happened during the file open...\n')
            exit()
        #---------------------------------
        # Read the file header information
        data = filehandler.read(180)
        head = SparkyFileHeader(data)
        print head.number_of_axis

        self.GetFileHeaderInformation(data)
        #---------------------------------
        # Read all axis information
        for AxisNumber in range(self._FileHeader_['Number of Axis']):
            datax = filehandler.read(128)
            self.GetAxisInformation(datax, AxisNumber)
            self._AxisOrder_.append(self.myaxis[AxisNumber]['Nucleus'][-1])
#        exit()
        #---------------------------------
        # Only 2D and 3D are ok
        if not self.NumberOfAxis in [2, 3]:
            print ('Sorry! The dimension of your spectrum (' +
                   str(self.NumberOfAxis) +
                   'D) is not handled by this program...\n')
            exit()
        #---------------------------------
        # Calculate the block size information
        Blocksize = 1
        for AxisNumber in range(self.NumberOfAxis):
            Blocksize *= self.myaxis[AxisNumber]['BlockSize']
        #---------------------------------
        # Read the data from the file
        Filecontent = filehandler.read()
        #---------------------------------
        # Close the file
        filehandler.close()
        #---------------------------------
        # Get the actual specral information
        if self.log:
            print 'File read has started',
        self._Spectrum_ = []
        # It can read 2D and 3D datafile
        if self.NumberOfAxis in [2,3]:
            eval('self._extract_'+str(self.NumberOfAxis)+'D_data(Filecontent, Blocksize)')
        #---------------------------------
        # Calculate a noise level for the spectrum
        self.CalculateNoiseLevel()
        #---------------------------------
        if self.log:
            print '100% file read is done.'
        return None
    ###########################
    def GetFileHeaderInformation(self, data):
        infos = struct.unpack('>10s 4c 9s 26s 80s 3x l 40s 4x',data)
        self._FileHeader_['Sparky ID'           ] = str(infos[0]).strip('\x00')
        self._FileHeader_['Number of Axis'      ] = ord(infos[1]) #
        self._FileHeader_['Number of Components'] = ord(infos[2]) # = 1 for real data
        self._FileHeader_['Encoding'            ] = ord(infos[3])
        self._FileHeader_['Version'             ] = ord(infos[4]) # = 2 for current format
        self._FileHeader_['Owner'               ] = str(infos[5]).strip('\x00')
        self._FileHeader_['Date'                ] = str(infos[6]).strip('\x00')
        self._FileHeader_['Comment'             ] = str(infos[7]).strip('\x00')
        self._FileHeader_['Seek Position'       ] = str(infos[8]).strip('\x00')
        self._FileHeader_['Scratch'             ] = str(infos[9]).strip('\x00')
        return None
    ###########################
    def GetAxisInformation(self, data, Number):
        infos = struct.unpack('>6s h 3I 6f 84s',data)
        self.myaxis[Number] = {}
        self.myaxis[Number]['Nucleus'               ] = str(infos[0]).strip('\x00') # nucleus name (1H, 13C, 15N, 31P, ...
        self.myaxis[Number]['Spectral Shift'        ] = infos[1] # to left or right shift
        self.myaxis[Number]['Number of Points'      ] = infos[2] # # of active data points - integer number of data points along this axis
        self.myaxis[Number]['Size'                  ] = infos[3] # total size of axis
        self.myaxis[Number]['BlockSize'             ] = infos[4] # # of points per cache block - integer tile size along this axis
        self.myaxis[Number]['Spectrometer frequency'] = infos[5] # MHz - float spectrometer frequency for this nucleus (MHz)
        self.myaxis[Number]['Spectral width'        ] = infos[6] # Hz  - float spectral width
        self.myaxis[Number]['xmtr frequency'        ] = infos[7] # transmitter offset (ppm) - float center of data (ppm)
        self.myaxis[Number]['Zero order'            ] = infos[8] # phase corrections
        self.myaxis[Number]['First order'           ] = infos[9] # phase corrections
        self.myaxis[Number]['First pt scale'        ] = infos[10] # scaling for first point
        self.myaxis[Number]['Extended'              ] = str(infos[11]).strip('\x00') #

        self.myaxis[Number]['Scale'] = []
        for i in range(0, int(self.myaxis[Number]['Number of Points']) + 1, 1):
            self.myaxis[Number]['Scale'].append(self._fid2ppm(i, infos[5], infos[7], infos[6], infos[2]))
        return None
    ###########################
    def _extract_2D_data(self, Filecontent, Blocksize):
        """
        """
        # First dimensional data
        FirstDimensionBlockSize         = self.myaxis[0]['BlockSize']
        FirstDimensionSpectralSize      = self.myaxis[0]['Size']
        # Second dimensional data
        SecondDimensionBlockSize        = self.myaxis[1]['BlockSize']
        SecondDimensionSpectralSize     = self.myaxis[1]['Size']
        # The number of blocks needed for a spectral size is
        # not necessary an integer number
        NumberOfBlocksInSecondDimension = (
            self._ceil(SecondDimensionSpectralSize /
                       float(SecondDimensionBlockSize)))
        print FirstDimensionBlockSize, SecondDimensionBlockSize, Blocksize
        exit()
        #---------------------------------
        # Rearrange the data from a list to an array
        for i_FirstDimension in range(FirstDimensionSpectralSize):
            # Print out info to follow the processing
            if self.log and i_FirstDimension % 50 == 0:
                print '{0:3.2f}%'.format(100.0 * i_FirstDimension
                                         / FirstDimensionSpectralSize),
            #---------------------------------
            BlockNumber = (i_FirstDimension / FirstDimensionBlockSize
                           * NumberOfBlocksInSecondDimension)
            PositionWithinBlock = (i_FirstDimension
                                   % FirstDimensionBlockSize
                                   * SecondDimensionBlockSize)
            # Concatenate the block portions in a list
            SpectralInfo1D = []
            #---------------------------------
            # Go through all second dimension protion to get a line
            for i_SecondDimension in range(NumberOfBlocksInSecondDimension):
                # If this is the last Block in line then the dimension is
                # not necessary the blocksize
                if i_SecondDimension < NumberOfBlocksInSecondDimension:
                    SecondDimension = SecondDimensionBlockSize
                else:
                    SecondDimension = (SecondDimensionSpectralSize
                                       % SecondDimensionBlockSize)
                #---------------------------------
                # The actual position within the block; 1 float number = 4 bytes
                pos = (4 * (Blocksize * (BlockNumber + i_SecondDimension)
                       + PositionWithinBlock))
                #---------------------------------
                # Unpack the data. Note that the coding is big endian ">"
                SpectralInfo1D += list(struct.unpack('>'+'f'*SecondDimension,
                                  Filecontent[pos : pos + 4 * SecondDimension]))
            #---------------------------------
            # Add a line into the spectrum
            self._Spectrum_.append(SpectralInfo1D)

        self.myaxis[0]['Actual size'] = len(self._Spectrum_)
        self.myaxis[1]['Actual size'] = len(self._Spectrum_[0])
        return None
    ###########################
    def _extract_3D_data(self, Filecontent, Blocksize):
        """
        """
        # Third dimensional data
        ThirdDimensionBlockSize         = self.myaxis[0]['BlockSize']
        ThirdDimensionSpectralSize      = self.myaxis[0]['Size']
        # Second dimensional data
        SecondDimensionBlockSize        = self.myaxis[1]['BlockSize']
        SecondDimensionSpectralSize     = self.myaxis[1]['Size']
        # First dimensional data
        FirstDimensionBlockSize         = self.myaxis[2]['BlockSize']
        FirstDimensionSpectralSize      = self.myaxis[2]['Size']
        #---------------------------------
        # The number of blocks needed for a spectral size is not necessary an integer number
        NumberOfBlocksInFirstDimension  = self._ceil(FirstDimensionSpectralSize /float(FirstDimensionBlockSize ))
        NumberOfBlocksInSecondDimension = self._ceil(SecondDimensionSpectralSize/float(SecondDimensionBlockSize))
        #---------------------------------
        # Rearrange the data from a list to an 3D array
        for i_ThirdDimension in range(ThirdDimensionSpectralSize):
            # Print out log information
            if i_ThirdDimension % 10 == 0:
                print '{0:3.2f}%'.format(100.0*i_ThirdDimension/ThirdDimensionSpectralSize),
            #---------------------------------
            BlockNumberDim3         = (i_ThirdDimension / ThirdDimensionBlockSize) * NumberOfBlocksInSecondDimension * NumberOfBlocksInFirstDimension
            PositionWithinBlockDim3 = (i_ThirdDimension % ThirdDimensionBlockSize) * SecondDimensionBlockSize        * FirstDimensionBlockSize
            #---------------------------------
            # Collect data of 2D in a variable
            SpectralInfo2D = []
            # Go through each block in 2D
            #for i_SecondDimension in range(SecondDimensionBlockSize * NumberOfBlocksInSecondDimension):
            for i_SecondDimension in range(SecondDimensionSpectralSize):
                #
                BlockNumberDim2         = BlockNumberDim3         + (i_SecondDimension / SecondDimensionBlockSize) * NumberOfBlocksInFirstDimension
                PositionWithinBlockDim2 = PositionWithinBlockDim3 + (i_SecondDimension % SecondDimensionBlockSize) * FirstDimensionBlockSize
                #---------------------------------
                # Collect data of 1D in a variable
                SpectralInfo1D = []
                # Go through each block in 1D
                for i_FirstDimension in range(NumberOfBlocksInFirstDimension):
                    # The last block size might be smaller than a blocksize
                    if i_FirstDimension < NumberOfBlocksInFirstDimension-1:
                        FirstDimension = FirstDimensionBlockSize
                    else:
                        FirstDimension = FirstDimensionSpectralSize % FirstDimensionBlockSize
                    #---------------------------------
                    # Position within block; 1 float number = 4 bytes
                    pos = 4 * (Blocksize * (BlockNumberDim2 + i_FirstDimension) + PositionWithinBlockDim2)
                    #---------------------------------
                    # Unpack the data. NOTE: big endian data storage ">"
                    SpectralInfo1D += list(struct.unpack('>'+'f'*FirstDimension,Filecontent[pos: pos + 4*FirstDimension]))
                #---------------------------------
                # Put each 1D slice into the 2D
                SpectralInfo2D.append(SpectralInfo1D)
            #---------------------------------
            # Store a 2D slice into the final array
            self._Spectrum_.append(SpectralInfo2D)
        self.myaxis[0]['Actual size'] = len(self._Spectrum_)
        self.myaxis[1]['Actual size'] = len(self._Spectrum_[0])
        self.myaxis[2]['Actual size'] = len(self._Spectrum_[0][0])
        return None
    ###########################
    def DataIntensity(self, position):
        if len(position) == 3:
            intensity =  (self._Spectrum_[position[0] % self.myaxis[0]['Actual size']]
                                         [position[1] % self.myaxis[1]['Actual size']]
                                         [position[2] % self.myaxis[2]['Actual size']])
        else:
            intensity =  (self._Spectrum_[position[0] % self.myaxis[0]['Actual size']]
                                         [position[1] % self.myaxis[1]['Actual size']])
        return intensity
    ###########################
    def distance(self, pos1, pos2):
        distance_value = 0.0
        for (p1, p2) in (pos1, pos2):
            distance_value += (p1 - p2)**2
        return math.sqrt(distance_value)
    ###########################
    def read_peaklist(self, PeaklistFilename, Info='[0:-1]', shifts = [0.0, 0.0, 0.0]):
        """
        Reads a sparky peaklist file
        """
        try:
            pfile = open(PeaklistFilename, 'r')
        except IOError:
            print 'Error opening ' + PeaklistFilename + '!!! Please check it!'
            exit()
        lines = pfile.readlines()
        pfile.close()
        for line in lines:
            if (len(line) > 12) and (not 'Assig' in line):
                data = line.split()
                key = data[0]
                self._Peaklistkeysorder_.append(key)
                self._Peaklist_[key] = {}
                order = key.split('-')
                if self.NumberOfAxis == 2:
                    self._Peaklist_[key]['Info'] = eval('order[0]' + Info)
                    self._Peaklist_[key][order[-1][-1]] = float(data[2]) + shifts[1]
                    self._Peaklist_[key][order[-2][-1]] = float(data[1]) + shifts[0]
                    self._Peaklist_[key]['Adjusted'] = 'red'
                    #
                    if not self._Peaklistchemicalshiftorder_:
                        self._Peaklistchemicalshiftorder_.append(order[-2][-1])
                        self._Peaklistchemicalshiftorder_.append(order[-1][-1])
                else:
                    self._Peaklist_[key]['Info'] = eval('order[0]'+Info)
                    self._Peaklist_[key][order[-1][-1]] = float(data[3]) + shifts[2]
                    self._Peaklist_[key][order[-2][-1]] = float(data[2]) + shifts[1]
                    self._Peaklist_[key][order[-3][-1]] = float(data[1]) + shifts[0]
                    self._Peaklist_[key]['Adjusted'] = 'red'
                    #
                    if not self._Peaklistchemicalshiftorder_:
                        self._Peaklistchemicalshiftorder_.append(order[-3][-1])
                        self._Peaklistchemicalshiftorder_.append(order[-2][-1])
                        self._Peaklistchemicalshiftorder_.append(order[-1][-1])
        return None
    ###########################
    def save_peaklist(self, filename):
        pfile = open(filename, 'w')
        for peak in self._Peaklistkeysorder_:
            line = peak
            for axis in self._Peaklistchemicalshiftorder_:
                line = ' '.join([line, str(self._Peaklist_[peak][axis])])
            line = ' '.join([line, str(self._Peaklist_[peak]['Intensity'])])
            line = ' '.join([line, '{0:5.2f}'.format(self._Peaklist_[peak]['Intensity']/self.Noiselevel)])
            line = ' '.join([line, '\n'])
            pfile.write(line)
        pfile.close()
        return None
    ###########################
    def extremes_finder(self, position, dimension, find_max = True):
        """
        find positive and negative extremes on the spectrum
        Parameters:
        ===========
            * position = spectrum staring position for the peak finding,
                         order must be same as in the spectrum
            * dimension = find local maximum or minimum in 2D or 3D
            * find_max = maximum or minimum finding
        Return:
        =======
            * local extreme
        """
        checklist = [[-1, 0, 0],[+1, 0, 0], # x
                     [ 0,-1, 0],[ 0,+1, 0], # y
                     [-1,-1, 0],[+1,-1, 0], # xy
                     [-1,+1, 0],[+1,+1, 0], # xy
                     [ 0, 0,-1],[ 0, 0,+1], # z
                     [-1, 0,-1],[+1, 0,-1], # xz
                     [-1, 0,+1],[+1, 0,+1], # xz
                     [ 0,-1,-1],[ 0,-1,-1], # yz
                     [ 0,+1,+1],[ 0,+1,+1]] # yz
        # If the dimension 2D, we find check the value in x,y otherwise in x,y,z
        if dimension == 2:
            checklist_size = 4
        else:
            checklist_size = len(checklist)
        # minimum or maximum finder
        finder_type = [['min','<'],['max','>']][find_max]
        # It goes till it finds a local maximum
        not_on_an_extreme_value = True
        while not_on_an_extreme_value:
            # check all values according to the checklist
            checked_values = []
            for check in checklist[0 : checklist_size]:
                checked_values.append(self.DataIntensity([pos + ch for (pos, ch) in zip(position[0 : dimension], check[0 : dimension])]))
            # if the position data is the most extreme, than we are done
            most_extreme_in_array = eval(eval('finder_type[0]') + '(checked_values)')
            if eval('self.DataIntensity(position)' + eval('finder_type[1]') + 'most_extreme_in_array'):
                not_on_an_extreme_value = False
            else:
                # modifiy the position to the highest
                checked_values_max_index = checked_values.index(most_extreme_in_array)
                for i in range(dimension):
                    position[i] += checklist[checked_values_max_index][i]
                    position[i] %= self.myaxis[i]['Actual size']
        return position
    ###########################
    def ClimbUpTheHill3D(self,ResidueKey, Modify = False, delta = [0.0,0.0,0.0]):
        if ResidueKey in self._Peaklistkeysorder_:
            #
            p = []
            original = []
            for i in range(3):
                p.append(int(round(delta[i])) + self._FidNumberbyAxis(self._Peaklist_[ResidueKey][self._Peaklistchemicalshiftorder_[i]],self._AxisOrder_.index(self._Peaklistchemicalshiftorder_[i])))
                original.append(int(round(delta[i])) + self._FidNumberbyAxis(self._Peaklist_[ResidueKey][self._Peaklistchemicalshiftorder_[i]],self._AxisOrder_.index(self._Peaklistchemicalshiftorder_[i])))

            checklist = [[-1, 0, 0],[+1, 0, 0], # x
                         [ 0,-1, 0],[ 0,+1, 0], # y
                         [ 0, 0,-1],[ 0, 0,+1], # z
                         [-1,-1, 0],[+1,-1, 0], # xy
                         [-1,+1, 0],[+1,+1, 0], # xy
                         [-1, 0,-1],[+1, 0,-1], # xz
                         [-1, 0,+1],[+1, 0,+1], # xz
                         [ 0,-1,-1],[ 0,-1,-1], # yz
                         [ 0,+1,+1],[ 0,+1,+1]] # yz
            Iteration = True
            while Iteration:
                tomb = []
                for ch in checklist:
                    tomb.append(self.DataIntensity([p[0] + ch[0],p[1] + ch[1],p[2] + ch[2]]))
                if self.DataIntensity(p) >= max(tomb):
                    Iteration = False
                else:
                    ti = tomb.index(max(tomb))
                    for i in range(3):
                        p[i] = (p[i] + checklist[ti][i]) % self.myaxis[i]['Size']

                    if ResidueKey == 'T680_N-C-H':
                        print 'PPM:',self._PPMNumberbyAxis(p[2],2)
            if Modify:
                for i in range(3):
                    self._Peaklist_[ResidueKey][self._Peaklistchemicalshiftorder_[i]] = self._PPMNumberbyAxis(p[i],self._AxisOrder_.index(self._Peaklistchemicalshiftorder_[i]))
        return p,original
    ###########################
    def AdjustAllPeakPositions3D(self):
        numberofpeaks = 0
        diff = [0.0, 0.0, 0.0]
        for key in self._Peaklistkeysorder_:
            if not (key in self._PeaklistDoNotCare_):
                a,b = self.ClimbUpTheHill3D(key)
                numberofpeaks += 1
                for i in range(3):
                    diff[i] += (a[i]-b[i])
        for i in range(3):
            diff[i] /= float(numberofpeaks)
        print diff
        for key in self._Peaklistkeysorder_:
            if not (key in self._PeaklistDoNotCare_):
                a,b = self.ClimbUpTheHill3D(key, Modify=True, delta= diff)



        return None
    ###########################
    def adjust_peaklist_2d(self):
        numberofpeaks = 0
        diff = [0.0, 0.0, 0.0]
        peaks = {}
        for key in self._Peaklistkeysorder_:
            if not (key in self._PeaklistDoNotCare_):
                position  = [self._FidNumberbyAxis(self._Peaklist_[key]['N'],'N'),
                             self._FidNumberbyAxis(self._Peaklist_[key]['H'],'H')]
                peaks[key] = {}
                peaks[key]['original'] = []
                peaks[key]['firsthit'] = []
                peaks[key]['secondhit'] = []
                #
                for pos in position:
                    peaks[key]['original'].append(pos)
                #
                peaks[key]['firsthit'] = self.extremes_finder(position, 2)
                numberofpeaks += 1
                for i in range(len(position)):
                    diff[i] += (peaks[key]['firsthit'][i] - peaks[key]['original'][i])
        for i in range(len(diff)):
            diff[i] /= numberofpeaks
            diff[i] = round(diff[i])
        #
        for key in self._Peaklistkeysorder_:
            if not (key in self._PeaklistDoNotCare_):
                position = []
                for i,pos in enumerate(peaks[key]['original']):
                    position.append(int(pos + diff[i]))
                peaks[key]['secondhit'] = self.extremes_finder(position, 2)
        #
        for i in range(len(self._Peaklistkeysorder_)):
            key = self._Peaklistkeysorder_[i]
            if not (key in self._PeaklistDoNotCare_):
                multiple = []
                j = 0
                while j < len(self._Peaklistkeysorder_):
                    key2 = self._Peaklistkeysorder_[j]
                    if (peaks[key]['secondhit'] == peaks[key2]['secondhit']) and (i != j):
                        multiple.append(j)
                    j += 1
                if not multiple:
                    # Unique peak found
                    peaks[key]['final'] = peaks[key]['secondhit']
                    peaks[key]['fit'] = 'black'
                else:
                    # Move the peak which is the closest
                    closest = True
                    for j in multiple:
                        key2 = self._Peaklistkeysorder_[j]
                        if (self.distance(peaks[key]['original'], peaks[key]['secondhit']) >=
                            self.distance(peaks[key2]['original'], peaks[key2]['secondhit'])):
                            closest = False
                    # if this peak is the most likely
                    if closest:
                        peaks[key]['final'] = peaks[key]['secondhit']
                        peaks[key]['fit'] = 'black'
                    else:
                        # If other peaks are closer, than just move to the average
                        peaks[key]['final'] = []
                        for (i, o) in enumerate(peaks[key]['original']):
                            peaks[key]['final'].append(int(o + diff[i]))
                        peaks[key]['fit'] = 'red'

#                print key, peaks[key]['original'], peaks[key]['firsthit'], peaks[key]['secondhit'],multiple, peaks[key]['final']
        for key in self._Peaklistkeysorder_:
            if not (key in self._PeaklistDoNotCare_):
                self._Peaklist_[key]['N'] =  self._PPMNumberbyAxis(peaks[key]['final'][0],'N')
                self._Peaklist_[key]['H'] =  self._PPMNumberbyAxis(peaks[key]['final'][1],'H')
                self._Peaklist_[key]['Adjusted'] = peaks[key]['fit']
                self._Peaklist_[key]['Intensity'] = self.DataIntensity(peaks[key]['final'])

        # TODO Fit the tip?
        return None
    ###########################
    def find_peak_1d(self, data, noiselevel):
        hits = []
        direction = True
        for i in range(len(data)-1):
            if data[i] > data[i+1] and data[i] > noiselevel and direction:
                hits.append(i)
                direction = False
            if data[i] < data[i+1]:
                direction = True
        return hits
    ###########################
    def find_peak_2d(self, data2d, noiselevel):
        hits = {}
        for i, data1d in enumerate(data2d):
            hit1d = self.find_peak_1d(data1d, noiselevel)
            for hit in hit1d:
                hits[' '.join(str(d) for d in self.extremes_finder([i, hit], 2))] = 0
        peaks = []
        for hit in hits:
            peaks.append(hit.split())
        return peaks
    ###########################
    def peak_finder(self, times_noiselevel):
        print 'Finding peaks...',
        peaklist = {}
        for i,peak in enumerate(self.find_peak_2d(self._Spectrum_,self.Noiselevel*times_noiselevel)):
            peaklist[i] = {}
            peaklist[i]['Info'] = str(i+1)
            peaklist[i]['N'] = self._PPMNumberbyAxis(float(peak[0]),'N')
            peaklist[i]['H'] = self._PPMNumberbyAxis(float(peak[1]),'H')
            peaklist[i]['Adjusted'] = 'black'
        print str(i + 1) + ' peaks found!'
        return peaklist
    ###########################
    def Plot1D(self, chemicalshift):
        dim = self._ppm2fid(chemicalshift,self.myaxis[0]['Spectrometer frequency'],self.myaxis[0]['xmtr frequency'],self.myaxis[0]['Spectral width'],self.myaxis[0]['Number of Points'])
        data = self._Spectrum_[dim]
        plt.figure()
        plt.plot(data)
        plt.show()
        return None
    ###########################
    def Plot1Dfid(self, fid):
        data = self._Spectrum_[fid]
        plt.figure()
        plt.plot(data)
        plt.show()
        return None
    ###########################
    def PPM_to_index(self,ppm,axisnumber):
        index = 0
        while (index < self.myaxis[axisnumber]['Number of Points']) and (self.myaxis[axisnumber]['Scale'][index] > ppm):
            index += 1
        return index
    ###########################
    def Limits_to_index(self, limits, axisnumber):
        if not limits:
            index_min = 0
            index_max = self.myaxis[axisnumber]['Number of Points']-1
        else:
            index_min = self.PPM_to_index(max(limits), axisnumber)
            index_max = self.PPM_to_index(min(limits), axisnumber)

        if index_max > self.myaxis[axisnumber]['Actual size']:
            index_max = self.myaxis[axisnumber]['Actual size']
        return index_min, index_max
    ###########################
    def spectrum_2d_slice(self, x_axis_min_index, x_axis_max_index,y_axis_min_index, y_axis_max_index, orderXY):
        highestvalue = 0.0
        lowestvalue = 0.0
        spectrum = []
        #---------------------------------
        # 2D
        if self.NumberOfAxis == 2:
            for y in range(y_axis_min_index, y_axis_max_index, 1):
                fid = []
                for x in range(x_axis_min_index, x_axis_max_index, 1):
                    if orderXY[0] == 'H':
                        value = self._Spectrum_[y][x]
                    else:
                        value = self._Spectrum_[x][y]
                    fid.append(value)
                    if highestvalue < value:
                        highestvalue = value
                    if lowestvalue > value:
                        lowestvalue = value
                spectrum.append(fid)
        return highestvalue, lowestvalue, spectrum
    ###########################
    def Plot_peaklist(self, Peaklist, x_min, x_max, y_min, y_max, orderXY):
        print 'Peaks on the plot:'
        number = 0
        for k in Peaklist:
            loc_x = Peaklist[k][orderXY[-2]]
            loc_y = Peaklist[k][orderXY[-1]]

            if ((x_min < loc_x) and (loc_x < x_max) and
                (y_min < loc_y) and (loc_y < y_max)):
                # TODO make is adjustable
                peak_info_pos_x = 0.0
                peak_info_pos_y = 0.0

#                plt.text(loc_x + peak_info_pos_x, loc_y + peak_info_pos_y, Peaklist[k]['Info'])
                number += 1
                print '{0:3d}.'.format(number),Peaklist[k]['Info'], loc_y, loc_x,
                if Peaklist[k]['Adjusted'] == 'black':
                    print 'ok'
                else:
                    print ''
                # TODO Make the dx,dy to be adjustable

                dx = 0.05
                dy = 0.2
                plt.gca().annotate(Peaklist[k]['Info'],
                                   xy=(loc_x,loc_y),
                                   color = Peaklist[k]['Adjusted'],
                                   xytext=(loc_x,loc_y - dy),
                                   arrowprops=dict(arrowstyle="-|>",
                                                   connectionstyle="arc3",
                                                   facecolor = Peaklist[k]['Adjusted']))

#
#                plt.plot([loc_x , loc_x + dx],[loc_y , loc_y + dy], 'k-')
#                plt.plot([loc_x , loc_x + dx],[loc_y , loc_y - dy], 'k-')

        return None
    ###########################
    def Plot(self, limits, orderXY='HN', color = [0, 0, 0], nf = True, peaklist = None):
        #
        axis_x = self._nucleustype2axisindex(orderXY[0])
        axis_y = self._nucleustype2axisindex(orderXY[1])
        # Figure out the limits
        x_axis_min_index, x_axis_max_index = self.Limits_to_index(limits[0],axis_x)
        y_axis_min_index, y_axis_max_index = self.Limits_to_index(limits[1],axis_y)

        x_scale = self.myaxis[axis_x]['Scale'][x_axis_min_index : x_axis_max_index]
        y_scale = self.myaxis[axis_y]['Scale'][y_axis_min_index : y_axis_max_index]

        # 2D
        if self.NumberOfAxis == 2:
            highestvalue, lowestvalue, spectrum = self.spectrum_2d_slice(x_axis_min_index, x_axis_max_index, y_axis_min_index, y_axis_max_index, orderXY)
        #---------------------------------
        mc = zcolor.MyColor()
        contour_start = self.Noiselevel
        contour_number = 25
        contour_factor = math.exp(math.log((highestvalue) /float(contour_start)) * 1.0/(float(contour_number)))
        contourlevels = [contour_start*contour_factor**i for i in range(contour_number)]
        contourcolors = [mc.series(i,contour_number,0,300) for i in range(contour_number)]

        print '#############################################'
        print '###   P L O T   #   P A R A M E T E R S   ###'
        print '#############################################'
        print 'Noise level   =', contour_start
        print 'Factor        =', contour_factor
        print 'Highest value =', highestvalue
        print 'Lowest value  =', lowestvalue
        print '#############################################'

        if nf:
            plt.figure()

        plt.contour(x_scale, y_scale, spectrum, contourlevels, colors = contourcolors)
        # plot negatives if needed!
        plt.contour(x_scale, y_scale, spectrum, [-1*i for i in contourlevels], colors = [[0.0,0.0,0.0] for i in range(contour_number)])

        if nf:
            plt.xlabel(self.myaxis[axis_x]['Nucleus']+' (ppm)',size=15)
            plt.ylabel(self.myaxis[axis_y]['Nucleus']+' (ppm)',size=15)
            plt.gca().invert_xaxis()
            plt.gca().invert_yaxis()

        # If peak labels are needed
        if self._Peaklist_ or peaklist:
            if not peaklist:
                self.Plot_peaklist(self._Peaklist_, x_scale[-1], x_scale[0], y_scale[-1], y_scale[0], orderXY)
            else:
                self.Plot_peaklist(peaklist, x_scale[-1], x_scale[0], y_scale[-1], y_scale[0], orderXY)

#        plt.show()
        return None
    ###########################
    def Plot_ori(self, limits, orderXY='HN', color = [0, 0, 0], Negatives = False, Peaklist=True, negcolors = 'o', ContourNumber = 15, Factor = 0.0, Noiselevel = 0, linewidth = 1.0, newfigure = True, figuresize=(8,5), figdpi=72, textsize=15):
        """
        Parameters:
            * limits           = an array of arrays with the PPM value limits, empty array means the whole spectral width
            * color            = one color value in [r,g,b] format eg. [1.0,0.0,0.0]
                               = array of color values (number must be the same as ContourNumber) eg. [[0.1,0.0,0.0],[0.2,0.0,0.0],...]
                               = built-in color eg. 'blue-cyan'
                               = built-in color + lighting info eg. ['g',0.5]
            * ContourNumber    = Number of contours on the figure
            * Factor           = factor between each contour level, provide 0.0 to calculate the value
            * Noiselevel       = If 0 is provided noise level is calculated from the sepctrum
            * linewidth        = contour line width, increase it when the zoom is high eg. 1.5
            * newfigure        = Boolean depending on the overlay plot option
            * figuresize       = figuresize in inch
            * figdpi           = dpi value, use 72 for screen, 300 for prints
            * textsize         = label size in pt eg. 12

        Examples:
            * Plot2D([[],[]],color = 'rainbow1')
            * Plot2D([[110,125],[7.2,9.5]],color = ['green',0.5], ContourNumber = 20, Factor = 1.2, Noiselevel = 100000, linewidth = 1.5, NumberOfThicksXY=[3,8], newfigure=False, figuresize=(5,5), figdpi=300, textsize=18)
        """

        ShowPeakLabelWithinPPM = [0.15,0.15,0.05] #NCH
        ShiftLabel = [0.0,0.0,0.0]
        #ShiftLabel = [0.05,0.05,0.02]
        CrossSize  = [0.05,0.05,0.01]
        Nucleuses  = ['N','C','H']

        #---------------------------------
        axisorder = []
        for ch in orderXY.upper():
            o = 0
            while (o < self.NumberOfAxis) and self.myaxis[o]['Nucleus'][-1] != ch:
                o += 1
            if o < self.NumberOfAxis:
                axisorder.append(o)
            else:
                print 'Please check the axes: ',orderXY
                exit()
        #---------------------------------
        # Check the limits to be within the spectrum range
        originallimits = limits
        lim = []
        for i in range(2):
            lim.append(self._AxisLimitCheck(axisorder[i],limits[i]))
        limits = lim
        if len(originallimits) == 3:
            limits.append(originallimits[2])
        #---------------------------------
        areamin = []
        areamax = []
        for i in range(2):
            areamax.append(self._ppm2fid(min(limits[i]),self.myaxis[axisorder[i]]['Spectrometer frequency'],self.myaxis[axisorder[i]]['xmtr frequency'],self.myaxis[axisorder[i]]['Spectral width'],self.myaxis[axisorder[i]]['Number of Points']))
            areamin.append(self._ppm2fid(max(limits[i]),self.myaxis[axisorder[i]]['Spectrometer frequency'],self.myaxis[axisorder[i]]['xmtr frequency'],self.myaxis[axisorder[i]]['Spectral width'],self.myaxis[axisorder[i]]['Number of Points']))
        #exit()
        # Get axis chemical shifts
        xscale = []
        for i in range(areamin[0],areamax[0]+1,1):
            xscale.append(self.myaxis[axisorder[0]]['Scale'][len(self.myaxis[axisorder[0]]['Scale'])-i-1])
#        print xscale[0],xscale[-1]
#        exit()
        yscale = []
        for i in range(areamin[1],areamax[1]+1,1):
            yscale.append(self.myaxis[axisorder[1]]['Scale'][len(self.myaxis[axisorder[1]]['Scale'])-i-1])
        print 'limits = ',areamin[0],areamax[0]
        #---------------------------------
        # Get the spectral information to plot
        highestvalue = 0.0
        area = []
        #---------------------------------
        # 2D
        if self.NumberOfAxis == 2:
            # Proton is on x
            if orderXY[0] == 'H':
                #
                for y in range(areamin[1],areamax[1]+1,1):
                    area.append(self._Spectrum_[y][areamin[0]:areamax[0]+1])
                    #
                    if max(self._Spectrum_[y][areamin[0]:areamax[0]+1]) > highestvalue:
                        highestvalue = max(self._Spectrum_[y][areamin[0]:areamax[0]+1])
            # Proton is on y
            else:
                for y in range(areamin[1],areamax[1]+1,1):
                    data = []
                    for x in range(areamin[0],areamax[0]+1,1):
                        value = self._Spectrum_[x][y]
                        data.append(value)
                        if value > highestvalue:
                            highestvalue = value
                    area.append(data)
        #---------------------------------
        # 3D
        if self.NumberOfAxis == 3:
            # Calculate the third dimension fid number
            zfid = self._ppm2fid(limits[2][0],self.myaxis[axisorder[2]]['Spectrometer frequency'],self.myaxis[axisorder[2]]['xmtr frequency'],self.myaxis[axisorder[2]]['Spectral width'],self.myaxis[axisorder[2]]['Number of Points'])
            # Extract the 2D from the 3D
            for y in range(areamin[1],areamax[1]+1,1):
                data = []
                for x in range(areamin[0],areamax[0]+1,1):
                    if orderXY[0:2] == 'HN':
                        value = self._Spectrum_[y][zfid][x]
                    elif orderXY[0:2] == 'HC':
                        value = self._Spectrum_[zfid][y][x]
                    elif orderXY[0:2] == 'NH':
                        value = self._Spectrum_[x][zfid][y]
                    elif orderXY[0:2] == 'NC':
                        value = self._Spectrum_[x][y][zfid]
                    elif orderXY[0:2] == 'CH':
                        value = self._Spectrum_[zfid][x][y]
                    elif orderXY[0:2] == 'CN':
                        value = self._Spectrum_[y][x][zfid]
                    else:
                        value = 0.0
                    # Store the value
                    data.append(value)
                    # Check whether it is the highest
                    if value > highestvalue:
                        highestvalue = value
                area.append(data)
        #---------------------------------
        # If the user did not set up a noise level, use the calculated one
        if Noiselevel == 0:
            contour_start   = self.Noiselevel
        else:
            contour_start   = Noiselevel
        contour_number  = ContourNumber
        #---------------------------------
        # If the user do not provide factor information
        if Factor == 0.0:
            # Calculcate based on the noise level and the highest peak height
            try:
                contour_factor = math.exp(math.log((highestvalue) /float(contour_start)) * 1.0/(float(contour_number)))
            except ValueError:
                contour_factor = 0.0
        # if the user provided the factor information
        else:
            contour_factor  = Factor
        #---------------------------------
        # Set the contour levels
        contourlevels = [contour_start*contour_factor**i for i in range(contour_number)]
        #---------------------------------
        # If the user provided a color
        contourcolors = self._ColorChoise(color,contour_number)
        if Negatives:
            # Colors
            negcontourcolors = self._ColorChoise(negcolors,contour_number)
            # Levels
            negcontourlevels = []
            for level in contourlevels:
                negcontourlevels.append(-1.0*level)
        #---------------------------------
        print '---------------'
        print self.myaxis[axisorder[0]]['Nucleus']+':',min(limits[0]),'-',max(limits[0])
        print self.myaxis[axisorder[1]]['Nucleus']+':',min(limits[1]),'-',max(limits[1])
        if self.NumberOfAxis == 3:
            print self.myaxis[axisorder[2]]['Nucleus']+':',limits[2][0]
        print 'Noise level   =', contour_start
        print 'Factor        =', contour_factor
        print 'Highest value =', highestvalue
        print '---------------'
        #---------------------------------
        # To be able to plot several figure on each other, the new figure is an option
        if newfigure:
            plt.figure(figsize=figuresize,dpi=figdpi)
        #---------------------------------
        # Generate the plot
        plt.contour(xscale,yscale,area,contourlevels,colors = contourcolors,linewidths = linewidth)
        if Negatives:
            plt.contour(xscale,yscale,area,negcontourlevels,colors = negcontourcolors,linewidths = linewidth)
        #---------------------------------
        # Invert the axes direction
        if newfigure:
            plt.gca().invert_xaxis()
            plt.gca().invert_yaxis()
        #---------------------------------
        # Put on axis label
        plt.xlabel(self.myaxis[axisorder[0]]['Nucleus']+' (ppm)',size=textsize)
        plt.ylabel(self.myaxis[axisorder[1]]['Nucleus']+' (ppm)',size=textsize)
        if self.NumberOfAxis == 3:
            plt.title(self.myaxis[axisorder[2]]['Nucleus']+': {0:6.3f} ppm'.format(limits[2][0]),size=textsize)
        #---------------------------------
        # If peak labels are needed
        if Peaklist and (self._Peaklist_ != {}):
            print 'Peaks on the plot:'
            for k in self._Peaklistkeysorder_:
                ItIsOn = True
                p = []
                for i in range(self.NumberOfAxis):
                    p.append(self._Peaklist_[k][self.myaxis[axisorder[i]]['Nucleus'][-1]])
                i = 0
                while (i < 2 ) and ItIsOn:
                    if (areamin[i] > p[i]) or (p[i] > areamax[i]):
                        ItIsOn = False
                    i += 1
                if self.NumberOfAxis == 3:
                    if abs(p[2] - limits[2][0]) > ShowPeakLabelWithinPPM[axisorder[2]]:
                        ItIsOn = False
                if ItIsOn:
                    print self._Peaklist_[k]['Info'],p[0],p[1],self._Peaklist_[k][Nucleuses[axisorder[2]]]
                    plt.text(p[0]-ShiftLabel[axisorder[0]],p[1]-ShiftLabel[axisorder[1]],self._Peaklist_[k]['Info'],size=textsize)
                    # Put on the crosspeak
                    dx = CrossSize[axisorder[0]]
                    dy = CrossSize[axisorder[1]]
                    #
                    plt.plot([p[0]-dx,p[0]+dx],[p[1]-dy,p[1]+dy],'k-')
                    plt.plot([p[0]-dx,p[0]+dx],[p[1]+dy,p[1]-dy],'k-')
        #
        return None
    ###########################
    def Show(self,FileName = ''):
        if FileName == '':
            plt.show()
        else:
            plt.savefig(FileName)
        return None
    ###########################
    def _AxisTicks(self,limits,number,PPMscale = True):
        # Calculate the step size
        step = abs(limits[0]-limits[1])/float(number-1)
        # Store the scales in data
        data = []
        for i in range(number):
            # if it is a ppm scale, then the values go down
            if PPMscale:
                value = max(limits)-i*step
            # if it is point scale then it goes up
            else:
                value = i*step
            #---------------------------------
            # if the value is extreme, then let 3 digits
            if int(value*1000) != value*1000:
                value = '{0:6.3f}'.format(value)
            data.append(value)
        return data
    ###########################
    def _AxisLimitCheck(self,Axisnumber,limits):
        # If there is no data provided, use the full spectrum
        if limits == []:
            limits = [-9.99E-99,+9.99E+99]
        # Store the data
        newlimits = []
        # Spectrum information
        ppmlimit = self.PPM_limit[Axisnumber]
        # Lower limit
        if min(ppmlimit) > min(limits):
            newlimits.append(self.myaxis[Axisnumber]['Scale'][1])
        else:
            newlimits.append(min(limits))
        # Upper limit
        if max(ppmlimit) < max(limits):
            newlimits.append(max(ppmlimit))
        else:
            newlimits.append(max(limits))
        return newlimits
    ###########################
    def _ppm2fid(self, ppm, Frequency, MiddlePPM, SpectralWidth, NumberOfPoints):
        return int((NumberOfPoints/2 - ((ppm-MiddlePPM) * Frequency * NumberOfPoints) / SpectralWidth) % NumberOfPoints)
    ###########################
    def _fid2ppm(self, fid, Frequency, MiddlePPM, SpectralWidth, NumberOfPoints):
        return MiddlePPM + (NumberOfPoints*SpectralWidth - 2*fid*SpectralWidth) / (2.0*Frequency*NumberOfPoints)
    ###########################
    def _nucleustype2axisindex(self, nucleus):
        axis = 0
        while (axis < self.NumberOfAxis) and (self.myaxis[axis]['Nucleus'][-1] != nucleus):
            axis += 1
        return axis
    ###########################
    def _axisindex2nucleustype(self, axisindex):
        return self.myaxis[axisindex]['Nucleus'][-1]
    ###########################
    def _FidNumberbyAxis(self, ppm, Axis):
        if type(Axis) == type(''):
            Axis = self._nucleustype2axisindex(Axis)
        return self._ppm2fid(ppm,
                             self.myaxis[Axis]['Spectrometer frequency'],
                             self.myaxis[Axis]['xmtr frequency'],
                             self.myaxis[Axis]['Spectral width'],
                             self.myaxis[Axis]['Number of Points'])
    ###########################
    def _PPMNumberbyAxis(self, fid, Axis):
        if type(Axis) == type(''):
            Axis = self._nucleustype2axisindex(Axis)
        return self._fid2ppm(fid,
                             self.myaxis[Axis]['Spectrometer frequency'],
                             self.myaxis[Axis]['xmtr frequency'],
                             self.myaxis[Axis]['Spectral width'],
                             self.myaxis[Axis]['Number of Points'])
    ###########################
    def _ceil(self, number):
        if number - int(number) != 0:
            number = int(number) + 1
        return int(number)
    ###########################
    def CalculateNoiseLevel(self,NumberOfDataPoints = 10000):
        Noise = 0.0
        # calculate the average level on a small subset of data
        average = 0.0
        for i in range(100):
            # 2D
            if self.NumberOfAxis == 2:
                average += abs(self._Spectrum_[random.randint(0,self.myaxis[0]['Number of Points']-1)][random.randint(0,self.myaxis[1]['Number of Points']-150)])
            # 3D
            if self.NumberOfAxis == 3:
                average += abs(self._Spectrum_[random.randint(0,self.myaxis[0]['Number of Points']-1)][random.randint(0,self.myaxis[1]['Number of Points']-1)][random.randint(0,self.myaxis[2]['Number of Points']-1)])
        average /= 100.0
        # Calculate the actual noise level
        numberofdata = 0
        sumofdata = 0.0
        highestvalue = 0.0
        i = 0
        while (i <= NumberOfDataPoints*2) and (numberofdata <= NumberOfDataPoints):
            # 2D
            if self.NumberOfAxis == 2:
                value = abs(self._Spectrum_[random.randint(0,self.myaxis[0]['Number of Points']-1)][random.randint(0,self.myaxis[1]['Number of Points']-150)])
            # 3D
            if self.NumberOfAxis == 3:
                value = abs(self._Spectrum_[random.randint(0,self.myaxis[0]['Number of Points']-1)][random.randint(0,self.myaxis[1]['Number of Points']-1)][random.randint(0,self.myaxis[2]['Number of Points']-1)])
            # Only count a value if that is not far from the average (= not a peak)
            if value < average * 5:
                numberofdata  += 1
                sumofdata += value
                average = sumofdata / float(numberofdata)
                if value > highestvalue:
                    highestvalue = value
            i += 1
        # Cut back from the highest to have a bit of noise
        Noise = highestvalue / 1.2
        # Assign the self.Noise to this value
        self.Noiselevel = Noise
        # Return the value as well
        return Noise
    ###########################
    def _ColorChoise(self,color,contour_number):
        if (type(color) == type([])) and (len(color) == 3):
            contourcolors = [color for _ in range(contour_number)]
        # if the user provided all the colors
        elif (type(color) == type([])) and (len(color) == contour_number):
            contourcolors = color
        # if the color is selected and light information is provided as well
        elif (type(color) == type([])) and (len(color) == 2):
            light = color[1]
            if (0.0 < light) or (light < 1.0):
                light = 1.0
            contourcolors = self.ColorSchemer(contour_number,color[0],light)
        # if there is no color information or built in colors has been selected
        else:
            contourcolors = self.ColorSchemer(contour_number,color)
        return contourcolors
    ###########################
    def ColorSchemer(self, Number, color, light = 1.0):
        data = []
        step = 1 / float(Number-1)
        for i in range(Number):
            element = [0.0,0.0,0.0]
            if (color == 'r') or (color == 'red'):
                element = [1.0,0.0,0.0]
            if (color == 'g') or (color == 'green'):
                element = [0.0,1.0,0.0]
            if (color == 'b') or (color == 'blue'):
                element = [0.0,0.0,1.0]
            #---------------------------------
            if (color == 'c') or (color == 'cyan'):
                element = [0.0,1.0,1.0]
            if (color == 'y') or (color == 'yellow'):
                element = [1.0,1.0,0.0]
            if (color == 'p') or (color == 'purple'):
                element = [1.0,0.0,1.0]
            #---------------------------------
            if (color == 'm') or (color == 'magenta'):
                element = [1.0,0.0,0.5]
            if (color == 'pi') or (color == 'pink'):
                element = [1.0,0.5,0.5]
            if (color == 'o') or (color == 'orange'):
                element = [1.0,0.5,0.0]
            #---------------------------------
            if (color == 'g1') or (color == 'grey1'):
                element = [0.1 for _ in range(3)]
            if (color == 'g2') or (color == 'grey2'):
                element = [0.2 for _ in range(3)]
            if (color == 'g3') or (color == 'grey3'):
                element = [0.3 for _ in range(3)]
            if (color == 'g4') or (color == 'grey4'):
                element = [0.4 for _ in range(3)]
            if (color == 'g5') or (color == 'grey5'):
                element = [0.5 for _ in range(3)]
            if (color == 'g6') or (color == 'grey6'):
                element = [0.6 for _ in range(3)]
            if (color == 'g7') or (color == 'grey7'):
                element = [0.7 for _ in range(3)]
            if (color == 'g8') or (color == 'grey8'):
                element = [0.8 for _ in range(3)]
            if (color == 'g9') or (color == 'grey9'):
                element = [0.9 for _ in range(3)]
            #---------------------------------
            if (color == 'w') or (color == 'white'):
                element = [1.0, 1.0, 1.0]
            #---------------------------------
            if (color == 'kr') or (color == 'black-red'):
                element = [0.0 + i * step, 0.0, 0.0]
            if (color == 'kg') or (color == 'black-green'):
                element = [0.0, 0.0 + i * step, 0.0]
            if (color == 'kb') or (color == 'black-blue'):
                element = [0.0, 0.0, 0.0 + i * step]
            #---------------------------------
            if (color == 'kc') or (color == 'black-cyan'):
                element = [0.0, 0.0 + i * step, 0.0 + i * step]
            if (color == 'ky') or (color == 'black-yellow'):
                element = [0.0 + i * step, 0.0 + i * step, 0.0]
            if (color == 'kp') or (color == 'black-purple'):
                element = [0.0 + i * step, 0.0, 0.0 + i * step]
            #---------------------------------
            if (color == 'km') or (color == 'black-magenta'):
                element = [0.0 + i * step, 0.0, 0.0 + (i / 2.0) * step]
            if (color == 'kpi') or (color == 'black-pink'):
                element = [0.0 + i * step, 0.0 + (i / 2.0) * step, 0.0 + (i / 2.0) * step]
            if (color == 'ko') or (color == 'black-orange'):
                element = [0.0 + i * step, 0.0 +(i / 2.0) * step, 0.0]
            #---------------------------------
            if (color == 'kw') or (color == 'black-white'):
                element = [0.0 + i * step, 0.0 + i * step, 0.0 + i * step]
            #---------------------------------
            if (color == 'rr') or (color == 'red-ring'):
                if i % 5 != 0:
                    element = [1.0, 0.0, 0.0]
                else:
                    element = [0.0, 0.0, 0.0]
            if (color == 'gr') or (color == 'green-ring'):
                if i % 5 != 0:
                    element = [0.0, 1.0, 0.0]
                else:
                    element = [0.0, 0.0, 0.0]
            if (color == 'br') or (color == 'blue-ring'):
                if i % 5 != 0:
                    element = [0.0, 0.0, 1.0]
                else:
                    element = [0.0, 0.0, 0.0]
            #---------------------------------
            if (color == 'red-yellow') or (color == 'rainbow1'):
                element = [1.0, 0.0 + i * step, 0.0]
            #---------------------------------
            if (color == 'blue-cyan') or (color == 'rainbow2'):
                element = [0.0, 0.0 + i * step, 1.0]
            #---------------------------------
            if (color == 'green-red') or (color == 'rainbow3'):
                element = [0.0 + i * step, 0.5 - (i / 2.0) * step, 0.0]
            #---------------------------------
            if type(light) != type(1.0):
                light = 1.0
            element = [element[c] * light for c in range(3)]
            #---------------------------------
            data.append(element)
        return data
    ###########################
    def _getNumberOfAxis(self):
        return len(self.myaxis.keys())
    ###########################
    def _getAxisInfo(self, field):
        info = []
        for axisnumber in range(self.NumberOfAxis):
            info.append(self.myaxis[axisnumber][field])
        return info
    ###########################
    def _getNucleus(self):
        return self._getAxisInfo('Nucleus')
    ###########################
    def _getFrequency(self):
        return self._getAxisInfo('Spectrometer frequency')
    ###########################
    def _getSpectralwidth(self):
        return self._getAxisInfo('Spectral width')
    ###########################
    def _getxmtrfreqency(self):
        return self._getAxisInfo('xmtr frequency')
    ###########################
    def _getscales(self):
        return self._getAxisInfo('Scale')
    ###########################
    def _getnumberofpoints(self):
        return self._getAxisInfo('Number of Points')
    ###########################
    def _getlimit(self):
        info = []
        for axisnumber in range(self.NumberOfAxis):
            info.append([self.myaxis[axisnumber]['Scale'][0],self.myaxis[axisnumber]['Scale'][-1]])
        return info
    ###########################
    NumberOfAxis   = property(_getNumberOfAxis)
    Nucleus        = property(_getNucleus)
    Frequency      = property(_getFrequency)
    SpectralWidth  = property(_getSpectralwidth)
    MiddlePPM      = property(_getxmtrfreqency)
    Scale          = property(_getscales)
    NumberOfPoints = property(_getnumberofpoints)
    PPM_limit      = property(_getlimit)
    #########################################

myspectrum= ZB_spectrum('13030_tcs_e.fid_1.ucsf')
print myspectrum.noise_level
# 1D proton plot
myspectrum.plot1d({'H':[] },'H')
myspectrum.show()

# Find peaks and plot a region with peaks and labeles 
peaks = myspecrum.peak_finder(1.5)
print peaks
myspectrum.plot([[6.8,10.2],[]],orderXY = 'HN', color = [],peaklist = peaks)
myspectrum.Show()

