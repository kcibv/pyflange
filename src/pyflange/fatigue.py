''' Fatigue calculation tools

This module defines functions and classes to support structural fatigue
calculations.

In particular, the module contains the following functions ...

- ``markov_matrix_from_SGRE_format(pathFile , unitFactor [optional])`` which reads 
a .mkv file from SGRE as markov matrix and converts in into a padas dataframe

... and the following ``FatigueCurve`` classes:

- ``SingleSlopeFatigueCurve``
- ``DoubleSlopeFatigueCurve``
-``SNCurve``

Each fatigue curve class exxposes the following methods:

- ``fatigue_curve.N(DS)`` returns the number of cycles corresponding to the
  given stress range DS
- ``fatigue_curve.DS(N)`` returns the stress range corresponding to the
  given number of cycles N
- ``fatigue_curve.damage(n, DS)`` returns the fatigue damage cumulated by
  a stress range DS repeated n times
'''

import numpy as np
import pandas as pd
from math import sqrt, pi, tan, exp, log10

from dataclasses import dataclass
import functools

def markov_matrix_from_SGRE_format(pathFile,unitFactor=1e3):
    ''' Reads .mkv file from SGRE as markov matrix and converts in into
    a padas dataframe:
        'Cycles' : Number of cylces
        'Mean' : mean bending moment
        'Range' : range of the bending moment
    
    '''
    MM=open(pathFile)
    
    MM=MM.readlines()
    
    
    mm_dict={'Cycles':[],
             'Mean':[],
             'Range':[]
             }
    
    rowMeans=False
    rowRanges=False
    countStartLines=0
    
    for row in MM:
    
        if row == '---------------------------\n':
            countStartLines+=1
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
    ''' A Wohler curve

    This is a base class for creating Wohler curves. It is not supposed to be
    instantiated directly.
    '''

    def N (self, DS):
        ''' Number of cycles

        Given a stress range DS, this function return the corresponding
        number of cycles that produce a fatigue failure.
        '''
        pass

    def DS (self, N):
        ''' Stress range

        Given a number of cycles, this function return the corresponding
        stress range that produce a fatigue failure.
        '''
        pass

    def damage (self, n, DS):
        ''' Fatigue damage

        Given a number of cycles n and a stress range DS, this function returns
        the dorresponding fatigue damage (D = n / N(DS)).
        '''
        return n / self.N(DS)


@dataclass
class SingleSlopeFatigueCurve (FatigueCurve):
    ''' Wohler curve with single logarithmic slope

    This class implements the FatigueCurve interface for a curve with single
    slope m.

    **Contructor parameters:**

    - `m` : float
        The logarithmic slope of the fatigue curve.

    - `DS_ref` : float
        Arbitrary reference stress range.

    - `N_ref` : float
        The number of cycles that produce failure under the stress range D_ref.

    **Attributes:**

    - `m` : float
        The logarithmic slope of the fatigue curve.

    - `DS_ref` : float
        Reference stress range.

    - `N_ref` : float
        Number of cycles at failure corresponding to the stress range D_ref.

    - `a` : float
        The Wohler curve constant (a = DS_ref**m * N_ref = DS**m * N)

    **Methods:**

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
    '''Multi-Slope Fatigue Curve

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
    ''' Wohler curve with double logarithmic slope

    This class implements the FatigueCurve interface for a curve with two
    slopes m1 and m2.

    **Contructor parameters:**

    - `m1` : float
        The logarithmic slope of the lower cycle values.

    - `m2` : float
        The logarithmic slope of the higher cycle values.

    - `DS12` : float
        The stress range where the two branches of the curve meet.

    - `N12` : float
        The number of cycles to failure corresponding to DS12.

    **Attributes:**

    - `m1` : float
        The logarithmic slope of the lower cycle values.

    - `m2` : float
        The logarithmic slope of the higher cycle values.

    - `DS12` : float
        The stress range where the two branches of the curve meet.

    - `N12` : float
        The number of cycles to failure corresponding to DS12.

    **Methods:**

    This class implements all the methods of FatigueCurve.
    '''

    def __init__ (self, m1, m2, DS12, N12):
        curve1 = SingleSlopeFatigueCurve(m1, DS12, N12)
        curve2 = SingleSlopeFatigueCurve(m2, DS12, N12)
        super().__init__(curve1, curve2)

@dataclass
class SNCurve(MultiSlopeFatigueCurve):
    ''' Wohler curve with double logarithmic slope for bolts

    This class implements the FatigueCurve interface for a curve with two
    slopes m1 and m2.

    - `m1` : float
        The logarithmic slope of the lower cycle values.

    - `m2` : float
        The logarithmic slope of the higher cycle values.
        
    - `N12` : float
        The knee point, where the slope changes to m2.

    - `DC` : float
        The damage category of the bolt.
    
    - `bolt_diameter` : float
        The nominal diameter of the bolt.
    
    - `size_factor_exponent` : float [optional]
        The exponent for the calculation of the size factor ks.
    
    - `gamma_M` : float [optional]
        The safty factor for the material.

    This class implements all the methods of FatigueCurve.
    '''
    
    m1 : float
    m2 : float
    N12 : float
    
    DC : float
    bolt_diameter : float
    size_factor_exponent : float = 0.0
    
    gamma_M : float = 1.0

    def cumulated_damage(self, df_markov):
        ks=min((30/(self.bolt_diameter*1000))**self.size_factor_exponent,1)
        DC_red=self.DC*ks/self.gamma_M
        
        DC_d = DC_red / ((self.N12/float(2e6))**(1/self.m1))
        curve1 = SingleSlopeFatigueCurve(self.m1, DC_d, self.N12)
        curve2 = SingleSlopeFatigueCurve(self.m2, DC_d, self.N12)
        super().__init__(curve1, curve2)
        
        df_markov['Damage']=self.damage(df_markov['Cycles'], df_markov['DS'])
        df_markov['Damage'].fillna(0)
        
        damage_sum=df_markov['Damage'].sum()

        return damage_sum
