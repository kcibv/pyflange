import pytest
from pyflange.fatigue import markov_matrix_from_SGRE_format, BoltFatigueCurve
from pyflange.flangesegments import PolynomialTFlangeSegment, bolt_markov_matrix
from pyflange.bolts import MetricBolt, HexNut
from pyflange.gap import gap_height_distribution
from math import pi
import os

# Units of measurement
m = 1
mm = 0.001*m

kg = 1
t = 1000*kg

s = 1

N = kg*m/s**2
kN = 1000*N

Pa = 1
MPa = 1e6*Pa
GPa = 1e9*Pa

rad = 1
deg = (pi/180)*rad

class TestFatigue:
    
    def calculate_damage(self,gap_angle, gap_shape_factor=1.0, tilt_angle=0):
        # Bolt
        M48 = MetricBolt(
            nominal_diameter = 48*mm,
            thread_pitch = 5*mm,
            shank_diameter_ratio = 44.752/48,
            shank_length = 150*mm,
            yield_stress = 900*MPa,
            ultimate_tensile_stress = 1000*MPa,
            stud = True)

        M48_hex_nut = HexNut(
            nominal_diameter = 48*mm,
            thickness = 64*mm,
            inscribed_diameter = 75*mm,
            circumscribed_diameter = 82.6*mm,
            bearing_diameter = 92*mm
        )

        # Polinomial Segment Model

        D = 7500*mm
        t_sh = 90*mm
        n = 200 # number of bolts
        gap_length = gap_angle * D/2
        gap = gap_height_distribution(D, 0.0014, gap_length)

        fseg = PolynomialTFlangeSegment(

            a = 62.5*mm,           # distance between inner face of the flange and center of the bolt hole
            b = 111.0*mm,           # distance between center of the bolt hole and center-line of the shell
            s = t_sh,               # shell thickness
            t = 120.0*mm,           # flange thickness
            R = D/2,                # shell outer curvature radius
            central_angle = 2*pi/n, # angle subtented by the flange segment arc

            Zg = -81.4*kN,     # load applied to the flange segment shell at rest
                                    # (normally dead weight of tower + RNA, divided by the number of bolts)

            bolt = M48,
            Fv = 928*kN,       # applied bolt preload

            Do = 52*mm,         # bolt hole diameter
            washer = None,      # no washer
            nut = M48_hex_nut,  # bolt nut

            gap_height = gap.ppf(0.95),             # maximum longitudinal gap height
            gap_angle = gap_angle,                  # longitudinal gap length
            gap_shape_factor = gap_shape_factor,    # scaling factor accounting for the gap shape

            tilt_angle = tilt_angle,    # flange tilt angle

            s_ratio = 1.0    # ratio of bottom shell thickness over tower shell thickness
        )

        markov_path=os.path.join(os.path.dirname(__file__), "..\\validations\\flangesegments.PolynomialTFlangeSegment\\tflange-example-markov.mkv")
        df_markov_shell=markov_matrix_from_SGRE_format(markov_path)
        df_markov_bolt = bolt_markov_matrix(fseg, df_markov_shell,
                                               bending_factor=0.601,
                                               macro_geometric_factor=1.0,
                                               mean_factor=1.3,
                                               range_factor=1.5)
        
        bfc=BoltFatigueCurve(M48.nominal_diameter)
        dmg=bfc.cumulated_damage(df_markov_bolt)

        return dmg   

    def test_fatigue(self):
        delta_max=0.03
        
        assert abs(1-self.calculate_damage( 30*deg, 1.0, 0*deg)/0.997) <= delta_max
        assert abs(1-self.calculate_damage( 60*deg, 1.0, 0*deg)/0.743) <= delta_max
        assert abs(1-self.calculate_damage( 90*deg, 1.0, 0*deg)/0.676) <= delta_max
        assert abs(1-self.calculate_damage( 120*deg, 1.0, 0*deg)/0.658) <= delta_max
        
        assert abs(1-self.calculate_damage( 30*deg, 1.5, 0*deg)/3.728) <= delta_max
        assert abs(1-self.calculate_damage( 60*deg, 1.5, 0*deg)/2.814) <= delta_max
        assert abs(1-self.calculate_damage( 90*deg, 1.5, 0*deg)/2.573) <= delta_max
        assert abs(1-self.calculate_damage( 120*deg, 1.5, 0*deg)/2.506) <= delta_max
        
        assert abs(1-self.calculate_damage( 30*deg, 1.0, 1.0*deg)/4.016) <= delta_max
        assert abs(1-self.calculate_damage( 60*deg, 1.0, 1.0*deg)/2.949) <= delta_max
        assert abs(1-self.calculate_damage( 90*deg, 1.0, 1.0*deg)/2.686) <= delta_max
        assert abs(1-self.calculate_damage( 120*deg, 1.0, 1.0*deg)/2.616) <= delta_max
