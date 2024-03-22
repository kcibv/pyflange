
import pytest
from pyflange.bolts import MetricBolt, StandardMetricBolt


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

        bolt = MetricBolt(
            nominal_diameter = 0.080,
            thread_pitch = 0.006,
            yield_stress = 900e6,
            ultimate_tensile_stress = 1000e6,
            shank_length = 0.270,
            shank_diameter_ratio = 76.1/80)

        assert round(bolt.axial_stiffness(0.400)/1e6) == 1831


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

        bolt = MetricBolt(
            nominal_diameter = 0.080,
            thread_pitch = 0.006,
            yield_stress = 900e6,
            ultimate_tensile_stress = 1000e6,
            shank_length = 0.270,
            shank_diameter_ratio = 76.1/80)

        assert round(bolt.bending_stiffness(0.400)/1e3) == 648


        bolt = MetricBolt(
            nominal_diameter = 0.080,
            thread_pitch = 0.006,
            yield_stress = 900e6,
            ultimate_tensile_stress = 1000e6,
            shank_length = 0.270,
            shank_diameter_ratio = 76.1/80,
            stud = True)

        assert round(bolt.bending_stiffness(0.400)/1e3) == 601
