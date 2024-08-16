
import pytest
from pyflange.bolts import MetricBolt, StandardMetricBolt, FlatWasher, ISOFlatWasher


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


    def test_shank_cross_section_area (self):
        bolt = MetricBolt(
            nominal_diameter = 0.016,
            thread_pitch = 0.002,
            yield_stress = 640e6,
            ultimate_tensile_stress = 800e6,
            shank_diameter_ratio = 2.0)

        assert round(bolt.shank_cross_section_area, 6) == 0.000804


    def test_tensile_cross_section_area (self):
        bolt = MetricBolt(
            nominal_diameter = 0.016,
            thread_pitch = 0.002,
            yield_stress = 640e6,
            ultimate_tensile_stress = 800e6)

        assert round(bolt.tensile_cross_section_area, 6) == 0.000157


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

    def outer_diameter (self):
        washer = FlatWasher(outer_diameter=10, inner_diameter=3, thickness=1)
        assert washer.outer_diameter == 10

    def inner_diameter (self):
        washer = FlatWasher(outer_diameter=10, inner_diameter=3, thickness=1)
        assert washer.inner_diameter == 3

    def thickness (self):
        washer = FlatWasher(outer_diameter=10, inner_diameter=3, thickness=1)
        assert washer.outer_diameter == 1

    def elastic_modulus (self):
        washer = FlatWasher(outer_diameter=10, inner_diameter=3, thickness=1)
        assert washer.elastic_modulus == 210e9

        washer = FlatWasher(outer_diameter=10, inner_diameter=3, thickness=1, elastic_modulus=100)
        assert washer.elastic_modulus == 100

    def poissons_ratio (self):
        washer = FlatWasher(outer_diameter=10, inner_diameter=3, thickness=1)
        assert washer.poissons_ratio == 0.3

        washer = FlatWasher(outer_diameter=10, inner_diameter=3, thickness=1, poissons_ratio=0.5)
        assert washer.poissons_ratio == 0.5

    def area (self):
        washer = FlatWasher(outer_diameter=10, inner_diameter=3, thickness=1)
        assert round(washer.area,3) == 71.471

    def axial_stiffness (self):
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
