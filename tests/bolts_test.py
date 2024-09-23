
import pytest
from pyflange.bolts import *



class TestBoltCrossSection:

    def test_diameter (self):
        csec = BoltCrossSection(10)
        assert csec.diameter == 10

    def test_area (self):
        csec = BoltCrossSection(10)
        assert round(csec.area, 6) == 78.539816

    def test_second_moment_of_area (self):
        csec = BoltCrossSection(10)
        assert round(csec.second_moment_of_area, 6) == 490.873852

    def test_elastic_section_modulus (self):
        csec = BoltCrossSection(10)
        assert round(csec.elastic_section_modulus, 6) == 98.174770



class TestMetricBolt:

    def test_designation (self):
        bolt = MetricBolt(
            nominal_diameter = 0.016,
            thread_pitch = 0.002,
            yield_stress = 640e6,
            ultimate_tensile_stress = 800e6)

        assert bolt.designation == "M16"


    def test_shank_diameter (self):

        # Ensure that shank diameter equals nominal diameter by defauls
        bolt = MetricBolt(
            nominal_diameter = 0.016,
            thread_pitch = 0.002,
            yield_stress = 640e6,
            ultimate_tensile_stress = 800e6)

        assert bolt.shank_diameter == 0.016

        # Ensure that shank diameter is a given ratio of the nominal diameter
        bolt = MetricBolt(
            nominal_diameter = 0.016,
            thread_pitch = 0.002,
            yield_stress = 640e6,
            ultimate_tensile_stress = 800e6,
            shank_diameter_ratio = 0.5)

        assert bolt.shank_diameter == 0.008


    def test_thread_height (self):
        bolt = MetricBolt(
            nominal_diameter = 0.016,
            thread_pitch = 0.002,
            yield_stress = 640e6,
            ultimate_tensile_stress = 800e6)

        assert round(bolt.thread_height, 5) == 0.00173


    def test_thread_basic_minor_diameter (self):
        bolt = MetricBolt(
            nominal_diameter = 0.016,
            thread_pitch = 0.002,
            yield_stress = 640e6,
            ultimate_tensile_stress = 800e6)

        assert round(bolt.thread_basic_minor_diameter, 6) == 0.013835


    def test_thread_basic_pitch_diameter (self):
        bolt = MetricBolt(
            nominal_diameter = 0.016,
            thread_pitch = 0.002,
            yield_stress = 640e6,
            ultimate_tensile_stress = 800e6)

        assert round(bolt.thread_basic_pitch_diameter, 6) == 0.014701


    def test_thread_basic_diameter (self):
        bolt = MetricBolt(
            nominal_diameter = 0.016,
            thread_pitch = 0.002,
            yield_stress = 640e6,
            ultimate_tensile_stress = 800e6)

        assert round(bolt.thread_minor_diameter, 5) == 0.01355


    def test_nominal_csec (self):
        bolt = MetricBolt(
            nominal_diameter = 0.016,
            thread_pitch = 0.002,
            yield_stress = 640e6,
            ultimate_tensile_stress = 800e6,
            shank_diameter_ratio = 2.0)

        csec = bolt.nominal_csec
        assert isinstance(csec, BoltCrossSection)
        assert round(csec.diameter, 6) == 0.016000


    def test_shank_csec (self):
        bolt = MetricBolt(
            nominal_diameter = 0.016,
            thread_pitch = 0.002,
            yield_stress = 640e6,
            ultimate_tensile_stress = 800e6,
            shank_diameter_ratio = 2.0)

        csec = bolt.shank_csec
        assert isinstance(csec, BoltCrossSection)
        assert round(csec.diameter, 6) == 0.032000


    def test_thread_csec (self):
        bolt = MetricBolt(
            nominal_diameter = 0.016,
            thread_pitch = 0.002,
            yield_stress = 640e6,
            ultimate_tensile_stress = 800e6)

        csec = bolt.thread_csec
        assert isinstance(csec, BoltCrossSection)
        assert round(csec.diameter, 6) == 0.014124


    def test_shear_modulus (self):
        bolt = MetricBolt(
            nominal_diameter = 0.016,
            thread_pitch = 0.002,
            yield_stress = 640e6,
            ultimate_tensile_stress = 800e6,
            elastic_modulus = 208e9,
            poissons_ratio = 0.3)

        assert bolt.shear_modulus == 80e9


    def test_ultimate_tensile_capacity (self):
        bolt = MetricBolt(
            nominal_diameter = 0.016,
            thread_pitch = 0.002,
            yield_stress = 640e6,
            ultimate_tensile_stress = 800e6)

        assert round(bolt.ultimate_tensile_capacity()/1000) == 90


    def test_axial_stiffness (self):

        # Hex head bolt
        bolt = MetricBolt(
            nominal_diameter = 0.080,
            thread_pitch = 0.006,
            yield_stress = 900e6,
            ultimate_tensile_stress = 1000e6,
            shank_length = 0.270,
            shank_diameter_ratio = 76.1/80)

        assert round(bolt.axial_stiffness(0.400)/1e6) == 1831

        # Stud bolt
        bolt = MetricBolt(
            nominal_diameter = 0.080,
            thread_pitch = 0.006,
            yield_stress = 900e6,
            ultimate_tensile_stress = 1000e6,
            shank_length = 0.270,
            shank_diameter_ratio = 76.1/80,
            stud = True)

        assert round(bolt.axial_stiffness(0.400)/1e6) == 1711


    def test_bending_stiffness (self):

        # Hex head bolt
        bolt = MetricBolt(
            nominal_diameter = 0.080,
            thread_pitch = 0.006,
            yield_stress = 900e6,
            ultimate_tensile_stress = 1000e6,
            shank_length = 0.270,
            shank_diameter_ratio = 76.1/80)

        assert round(bolt.bending_stiffness(0.400)/1e3) == 648

        # Stud bolt
        bolt = MetricBolt(
            nominal_diameter = 0.080,
            thread_pitch = 0.006,
            yield_stress = 900e6,
            ultimate_tensile_stress = 1000e6,
            shank_length = 0.270,
            shank_diameter_ratio = 76.1/80,
            stud = True)

        assert round(bolt.bending_stiffness(0.400)/1e3) == 601


    # DEPRECATED

    def test_tensile_cross_section_area (self):
        bolt = MetricBolt(
            nominal_diameter = 0.016,
            thread_pitch = 0.002,
            yield_stress = 640e6,
            ultimate_tensile_stress = 800e6)

        assert round(bolt.tensile_cross_section_area, 6) == 0.000157

    def test_shank_cross_section_area (self):
        bolt = MetricBolt(
            nominal_diameter = 0.016,
            thread_pitch = 0.002,
            yield_stress = 640e6,
            ultimate_tensile_stress = 800e6,
            shank_diameter_ratio = 2.0)

        assert round(bolt.shank_cross_section_area, 6) == 0.000804





class TestStandardMetricBolt:

    def test_designation (self):
        bolt = StandardMetricBolt("M16", "8.8")
        assert bolt.designation == "M16"


    def test_shank_diameter (self):

        # Ensure that shank diameter equals nominal diameter by defauls
        bolt = StandardMetricBolt("M16", "8.8")
        assert bolt.shank_diameter == 0.016

        # Ensure that shank diameter is a given ratio of the nominal diameter
        bolt = StandardMetricBolt("M16", "8.8", shank_diameter_ratio = 0.5)
        assert bolt.shank_diameter == 0.008


    def test_shank_cross_section_area (self):
        bolt = StandardMetricBolt("M16", "8.8", shank_diameter_ratio = 2.0)
        assert round(bolt.shank_cross_section_area, 6) == 0.000804


    def test_tensile_cross_section_area (self):
        bolt = StandardMetricBolt("M16", "8.8")
        assert round(bolt.tensile_cross_section_area, 6) == 0.000157


    def test_shear_modulus (self):
        bolt = StandardMetricBolt("M16", "8.8")
        assert round(bolt.shear_modulus/1e9, 2) == 80.77



class TestFlatWasher:

    def test_outer_diameter (self):
        washer = FlatWasher(outer_diameter=10, inner_diameter=3, thickness=1)
        assert washer.outer_diameter == 10

    def test_inner_diameter (self):
        washer = FlatWasher(outer_diameter=10, inner_diameter=3, thickness=1)
        assert washer.inner_diameter == 3

    def test_thickness (self):
        washer = FlatWasher(outer_diameter=10, inner_diameter=3, thickness=1)
        assert washer.thickness == 1

    def test_elastic_modulus (self):
        washer = FlatWasher(outer_diameter=10, inner_diameter=3, thickness=1)
        assert washer.elastic_modulus == 210e9

        washer = FlatWasher(outer_diameter=10, inner_diameter=3, thickness=1, elastic_modulus=100)
        assert washer.elastic_modulus == 100

    def test_poissons_ratio (self):
        washer = FlatWasher(outer_diameter=10, inner_diameter=3, thickness=1)
        assert washer.poissons_ratio == 0.3

        washer = FlatWasher(outer_diameter=10, inner_diameter=3, thickness=1, poissons_ratio=0.5)
        assert washer.poissons_ratio == 0.5

    def test_area (self):
        washer = FlatWasher(outer_diameter=10, inner_diameter=3, thickness=1)
        assert round(washer.area,3) == 71.471

    def test_axial_stiffness (self):
        washer = FlatWasher(outer_diameter=10, inner_diameter=3, thickness=2)
        assert round(washer.axial_stiffness, 3) == 7504479451262.618



def test_ISOFlatWasher ():
    washer = ISOFlatWasher("M16")
    assert isinstance(washer, FlatWasher)
    assert washer.outer_diameter == 0.030
    assert washer.inner_diameter == 0.017
    assert washer.thickness == 0.003
    assert washer.elastic_modulus == 210e9
    assert washer.poissons_ratio == 0.3



class TestHexNut:

    def test_nominal_diameter (self):
        nut = HexNut(nominal_diameter=16, thickness=10, inscribed_diameter=20, circumscribed_diameter=30, bearing_diameter=40)
        assert nut.nominal_diameter == 16

    def test_thickness (self):
        nut = HexNut(nominal_diameter=16, thickness=10, inscribed_diameter=20, circumscribed_diameter=30, bearing_diameter=40)
        assert nut.thickness == 10

    def test_inscribed_diameter (self):
        nut = HexNut(nominal_diameter=16, thickness=10, inscribed_diameter=20, circumscribed_diameter=30, bearing_diameter=40)
        assert nut.inscribed_diameter == 20

    def test_circumscribed_diameter (self):
        nut = HexNut(nominal_diameter=16, thickness=10, inscribed_diameter=20, circumscribed_diameter=30, bearing_diameter=40)
        assert nut.circumscribed_diameter == 30

    def test_bearing_diameter (self):
        nut = HexNut(nominal_diameter=16, thickness=10, inscribed_diameter=20, circumscribed_diameter=30, bearing_diameter=40)
        assert nut.bearing_diameter == 40

    def test_elastic_modulus (self):
        nut = HexNut(nominal_diameter=16, thickness=10, inscribed_diameter=20, circumscribed_diameter=30, bearing_diameter=40)
        assert nut.elastic_modulus == 210e9

        nut = HexNut(nominal_diameter=16, thickness=10, inscribed_diameter=20, circumscribed_diameter=30, bearing_diameter=40, elastic_modulus=100)
        assert nut.elastic_modulus == 100

    def test_poissons_ratio (self):
        nut = HexNut(nominal_diameter=16, thickness=10, inscribed_diameter=20, circumscribed_diameter=30, bearing_diameter=40)
        assert nut.poissons_ratio == 0.3

        nut = HexNut(nominal_diameter=16, thickness=10, inscribed_diameter=20, circumscribed_diameter=30, bearing_diameter=40, poissons_ratio=0.5)
        assert nut.poissons_ratio == 0.5



def test_ISOHexNut ():
    nut = ISOHexNut("M16")
    assert isinstance(nut, HexNut)
    assert nut.nominal_diameter == 0.016
    assert nut.thickness == 0.0148
    assert nut.inscribed_diameter == 0.024
    assert nut.circumscribed_diameter == 0.02675
    assert nut.bearing_diameter == 0.0225



def test_RoundNut ():
    nut = RoundNut("M42")
    assert isinstance(nut, HexNut)
    assert nut.nominal_diameter == 0.042
    assert nut.thickness == 0.042
    assert nut.inscribed_diameter == 0.065
    assert nut.circumscribed_diameter == 0.07130
    assert nut.bearing_diameter == 0.078
