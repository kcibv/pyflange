.. _example_montecarlo_simulation:

Example: Montecarlo Simulation
============================== 

This example shows how the Polynomial L-Flange Segment model implementation
can be used to predict a random serie of actualizations of bolt force and
bolt moment, based on the random imput of the following parameters.

- Gap Length, assumed normally distributed with mean value `30°` and COV `10%`
- Gap Height, assumed log-normally distributed according to IEC 61400-6 AMD1, section 6.7.5.2
- Bolt preload, assumed normally distributed with mean value `2876 kN` and COV `10%`

All the other parameters are assumed to be deterministically known, and defined as 
follows:

- Bolt `M80` with meterial grade `10.9`
- Distance between inner face of the flange and center of the bolt hole: `a = 232.5 mm`
- Distance between center of the bolt hole and center-line of the shell: `b = 166.5 mm`
- Shell thickness: `s = 72 mm`
- Flange thickness: `t = 200 mm`
- Flange outer diameter: `D = 7.5 m`
- Bolt hole diameter: `Do = 86 mm`
- Washer outer diameter: `Dw = 140 mm`

Let's first define a function that, given an actualization of gap length, gap height and
bolt preload generates the corresponding FlangeSegment object.

.. code-block:: python

   # Imports
   from math import pi
   from pyflange.bolts import StandardMetricBolt, ISOFlatWasher, ISOHexNut
   from pyflange.flangesegments import PolynomialLFlangeSegment

   # Define some units for the sake of readability
   m = 1
   mm = 0.001*m
   N = 1
   kN = 1000*N

   # Create the fastener parts
   M80_bolt = StandardMetricBolt("M80", "10.9", shank_length=270*mm stud=True)
   M80_washer = ISOFlatWasher("M80")
   M80_nut    = ISOHexNut("M80")
   Nb = 120   # Number of bolts

   # Flange Segment Constructor
   def create_flange_segment (gap_angle, gap_height, bolt_preload):
      return PolynomialLFlangeSegment(
         a = 232.5*mm,              # distance between inner face of the flange and center of the bolt hole
         b = 166.5*mm,              # distance between center of the bolt hole and center-line of the shell
         s = 72*mm,                 # shell thickness
         t = 200.0*mm,              # flange thickness
         R = 7.5*m / 2,             # shell outer curvature radius
         central_angle = 2*pi/Nb,   # angle subtented by the flange segment arc

         Zg = -15000*kN / Nb,    # load applied to the flange segment shell at rest

         bolt = M80_bolt,        # bolt object created above
         Fv = bolt_preload,      # applied bolt preload

         Do = 86*mm,             # bolt hole diameter
         washer = M80_washer,    # washer object created above
         nut = M80_nut,          # nut object created above

         gap_height = gap_height,   # maximum longitudinal gap height
         gap_angle = gap_angle      # longitudinal gap length
      )

Next, we define the stochastic variables `gap_angle`, `gap_height` and `bolt_pretension`:

.. code-block:: python

   # Gap Height Log-Normal Distribution
   from pyflange.gap import gap_height_distribution
   D = 7.5*m
   gap_length = pi/6 * D/2
   gap_height_dist = gap_height_distribution(7.5*m, 1.4*mm/m, gap_length)

   # Gap angle distribution
   from scipy.stats import normal
   mean = pi/6
   std = 0.10 * mean
   gap_angle_dist = norm(loc=mean, scale=std)

   # Bolt pretension distribution
   mean = 2876*kN
   std = 0.10 * mean
   bolt_preload_dist = norm(loc=mean, scale=std)


Next we generate random actualizations of the stochastic parameters and evaluate the
corresponding values of Fs(Z) and Ms(Z), in discrete form.

.. code-block:: python

   # Let's define the discrete domain Z of the Fs(Z) and Ms(Z) functions we
   # want to determine. We define Z as an array of 1000 items, linearly
   # spaced between -1500 kN and 2100 kN.
   import numpy as np 
   Z = np.linspace(-1500*kN, 2100*kN, 1000)

   # Let's generating 25000 actualizations of Fs(Z) and Ms(Z) and store them in
   # two 1000x25000 matrices, where each row is an actualization of the discrete
   # image of Z through Fs and Ms.
   
   Fs = np.array([])    # Initialize Fs with an empty matrix
   Ms = np.arrat([])    # Initialize Ms with an empty matrix

   for i in range(25000):

      # Generate a random gap height
      gap_height = gap_height_dist.rvs()

      # Generate a random gap angle
      gap_angle = gap_angle_dist.rvs()

      # Generate a random bolt pretension
      bolt_preload = bolt_preload_dist.rvs()

      # Generate the corresponding random FlangeSegment actualization, using the
      # factory function defined above
      fseg = create_flange_segment(gap_angle, gap_height, bolt_preload)

      # Generate the Fs image of Z and store it in the Fs matrix
      Fs.append( fseg.bolt_axial_force(Z) ) 

      # Generate the Ms image of Z and store it in the Ms matrix
      Ms.append( fseg.bolt_bending_moment(Z) ) 

The generated data in `Fs` and `Ms` can be then used to fit a distribution to
the data.

