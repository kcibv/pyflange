# pyFlange - python library for large flanges design
# Copyright (C) 2024  KCI The Engineers B.V.,
#                     Siemens Gamesa Renewable Energy B.V.,
#                     Nederlandse Organisatie voor toegepast-natuurwetenschappelijk onderzoek TNO.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License, as published by
# the Free Software Foundation, either version 3 of the License, or any
# later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License version 3 for more details.
#
# You should have received a copy of the GNU General Public License
# version 3 along with this program.  If not, see <https://www.gnu.org/licenses/>.

''' Fatigue calculation tools.

This module defines functions and classes to support structural fatigue
calculations.

In particular, the module contains the following functions ...

- `markov_matrix_from_SGRE_format(pathFile , unitFactor [optional])` which reads
  a .mkv file from SGRE as markov matrix and converts in into a pandas dataframe

... and the following `FatigueCurve` classes:

- `SingleSlopeFatigueCurve`
- `DoubleSlopeFatigueCurve`
- `SNCurve`

Each fatigue curve class exxposes the following methods:

- `fatigue_curve.N(DS)` returns the number of cycles corresponding to the
  given stress range DS
- `fatigue_curve.DS(N)` returns the stress range corresponding to the
  given number of cycles N
- `fatigue_curve.damage(n, DS)` returns the fatigue damage cumulated by
  a stress range DS repeated n times
'''

import numpy as np
import pandas as pd
from math import sqrt, pi, tan, exp, log10

from dataclasses import dataclass
import functools


def markov_matrix_from_SGRE_format (pathFile, unitFactor=1e3):
    ''' Reads a .mkv file into a pandas.DataFrame object.

    Reads a Markov matrix from a SGRE .mkv file and converts in into
    a padas dataframe having the collowing columns:

    - 'Cycles' : Number of cylces
    - 'Mean'   : mean bending moment
    - 'Range'  : range of the bending moment

    Args:
        pathFile (str): The path of the .mkv file to be read

        unitFactor (float): A scalind factor to be applied to the moment values.
            Useful for unit conversion.

    Returns:
        markov_matrix (pandas.DataFrame): the pandas DataFrame representation
            of the Markov natrix contained in the .mkv file.
    '''

    with open(pathFile) as mkv_file:
        MM = mkv_file.readlines()

    mm_dict={'Cycles':[],
             'Mean':[],
             'Range':[]
             }
    
    rowMeans = False
    rowRanges = False
    countStartLines = 0
    
    for row in MM:
    
        if row == '---------------------------\n':
            countStartLines += 1
            if countStartLines == 2:
                rowMeans=True
                continue
                
        if rowMeans:
            
            rowValues=row.replace('\n','').split(' ')
            meanValues = [e for e in rowValues if e not in (' ')]
            #meanValues.pop(0)
            rowMeans=False
            rowRanges=True
            continue
        
        if rowRanges:
            
            rowValues=row.replace('\n','').split(' ')
            rowValues = [e for e in rowValues if e not in (' ')]
            rangeValue=rowValues[0]
            
            for i in range(1,len(rowValues)):
                if float(rowValues[i]) == 0.0: continue
                #moment
                mm_dict['Cycles'].append(float(rowValues[i]))
                mm_dict['Mean'].append(float(meanValues[i])*unitFactor)
                mm_dict['Range'].append(float(rangeValue)*unitFactor)
        
    return pd.DataFrame(mm_dict) 



class FatigueCurve:
    ''' A Wöhler curve.

    This is a base class for creating Wohler curves. It is not supposed to be
    instantiated directly.
    '''

    def N (self, DS):
        ''' Number of cycles.

        Given a stress range DS, this function return the corresponding
        number of cycles that produce a fatigue failure.
        '''
        pass

    def DS (self, N):
        ''' Stress range.

        Given a number of cycles, this function return the corresponding
        stress range that produce a fatigue failure.
        '''
        pass

    def damage (self, n, DS):
        ''' Fatigue damage.

        Given a number of cycles n and a stress range DS, this function returns
        the dorresponding fatigue damage (D = n / N(DS)).
        '''
        return n / self.N(DS)

    def cumulated_damage (self, markov_matrix):
        ''' Cumulated damage according to the Miner's rule

        Args:
            markov_matrix (pandas.DataFrame): This is the load history expressed as a
                Markov matrix, encoded in a pandas DataFrame having three columns:
                - `Cycles`: containing the number of cycles;
                - `Mean`: containing the mean stress in Pascal;
                - `Range`: containing the stress range in Pascal.

        Returns:
            damage (float): the cumulated fatigue damage produced in the
                detail represented by this fatigue curev, by the load history
                represented by the passed Markov matrix.
        '''

        n = markov_matrix['Cycles'].to_numpy()  # array of number of cycles
        DS = markov_matrix['Range'].to_numpy()     # array of stress ranges
        D = self.damage(n, DS)                  # array of damages
        return np.nansum(D)                     # total damage



@dataclass
class SingleSlopeFatigueCurve (FatigueCurve):
    ''' Wöhler curve with single logarithmic slope.

    This class implements the FatigueCurve interface for a curve with single
    slope m.

    Args:
        m (float): The logarithmic slope of the fatigue curve.

        DS_ref (float): Arbitrary reference stress range.

        N_ref (float): The number of cycles that produce failure under the
            stress range D_ref.

    All the constructor parameters are also available as instance attributes
    (i.e. `scn.m`, `snc.DS_ref`, `snc.N_ref`).

    This class implements all the methods of FatigueCurve.
    '''

    m: float
    DS_ref: float
    N_ref: float


    @functools.cached_property
    def a (self):
        return self.DS_ref ** self.m * self.N_ref

    def N (self, DS):
        return self.a / DS**self.m

    def DS (self, N):
        return (self.a / N)**(1/self.m)



class MultiSlopeFatigueCurve(FatigueCurve):
    '''Multi-Slope Fatigue Curve.

    This class is a FatigueCurve with multiple slopes.
    It takes any number of SingleSlopeFatigueCurve objects as arguments.
    '''

    def __init__(self, *fatigue_curves):
        self.curves = fatigue_curves

    def N(self, DS):
        return np.maximum.reduce([curve.N(DS) for curve in self.curves])

    def DS(self, N):
        return np.maximum.reduce([curve.DS(N) for curve in self.curves])



class DoubleSlopeFatigueCurve (MultiSlopeFatigueCurve):
    ''' Wöhler curve with double logarithmic slope.

    This class implements the FatigueCurve interface for a curve with two
    slopes m1 and m2.

    Args:
        m1 (float): The logarithmic slope of the lower cycle values.

        m2 (float): The logarithmic slope of the higher cycle values.

        DS12 (float): The stress range where the two branches of the curve meet.

        N12 (float): The number of cycles to failure corresponding to DS12.

    All the constructor parameters are also available as instance attributes
    (i.e. `scn.m1`, `snc.m2`, `snc.DS12`, `snc.N12`).

    This class implements all the methods of FatigueCurve.
    '''

    def __init__ (self, m1, m2, DS12, N12):
        curve1 = SingleSlopeFatigueCurve(m1, DS12, N12)
        curve2 = SingleSlopeFatigueCurve(m2, DS12, N12)
        super().__init__(curve1, curve2)



class BoltFatigueCurve (DoubleSlopeFatigueCurve):
    ''' Bolt Fatigue Curve according to IEC 61400-6 AMD1.

    Given a bolt diameter, creates a DoubleSlopeFatigueCurve having
    logaritmic slopes m1=3 and m2=5 and change-of-slope at 2 milion cycles
    and stress range depending on the bolt diameter as specified by
    IEC 61400-6 AMD1.

    Args:
        diameter (float): The bolt diameter in meters.

        gamma_M (float): The material factor.

    Thuis class inherits all the properties and methods of the
    `DoubleSlopeFatigueCurve` class.
    '''

    def __init__ (self, diameter, gamma_M=1.1):
        N12 = 2.0e6    # knee point
        m1 = 3
        m2 = 5
        if diameter <= 0.030:
            DSc = 50e6    # reference stress range, in Pa
        elif diameter <= 0.072:
            DSc = 50e6 * (0.030/diameter)**0.1
        else:
            DSc = 50e6 * (0.030/diameter)**0.1 * (0.072/diameter)**0.25

        # Delegate the rest of the initialization to the parent class
        super().__init__(m1, m2, DSc/gamma_M, N12)


