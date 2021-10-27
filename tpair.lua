-- Auto-generated by prep-gas on: 25-Oct-2021 13:53:09

model = 'ThermallyPerfectGas'
species = {'N2', 'O2', }

physical_model = 'thermally-perfect-gas'
energyModes = {'equilibrium'}
db = {}
db['N2'] = {}
db['N2'].atomicConstituents = { N=2, }
db['N2'].charge = 0
db['N2'].M = 2.80134000e-02
db['N2'].sigma = 3.62100000
db['N2'].epsilon = 97.53000000
db['N2'].Lewis = 1.15200000
db['N2'].thermoCoeffs = {
  origin = 'CEA',
  nsegments = 3, 
  T_break_points = { 200.00, 1000.00, 6000.00, 20000.00, },
  T_blend_ranges = { 400.0, 1000.0, },
  segment0 = {
    2.210371497e+04,
   -3.818461820e+02,
    6.082738360e+00,
   -8.530914410e-03,
    1.384646189e-05,
   -9.625793620e-09,
    2.519705809e-12,
    7.108460860e+02,
   -1.076003744e+01,
  },
  segment1 = {
    5.877124060e+05,
   -2.239249073e+03,
    6.066949220e+00,
   -6.139685500e-04,
    1.491806679e-07,
   -1.923105485e-11,
    1.061954386e-15,
    1.283210415e+04,
   -1.586640027e+01,
  },
  segment2 = {
    8.310139160e+08,
   -6.420733540e+05,
    2.020264635e+02,
   -3.065092046e-02,
    2.486903333e-06,
   -9.705954110e-11,
    1.437538881e-15,
    4.938707040e+06,
   -1.672099740e+03,
  },
}
db['N2'].viscosity = {
   model = 'CEA',
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      A =  6.2526577e-01,
      B = -3.1779652e+01,
      C = -1.6407983e+03,
      D =  1.7454992e+00,
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  8.7395209e-01,
      B =  5.6152222e+02,
      C = -1.7394809e+05,
      D = -3.9335958e-01,
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  8.8503551e-01,
      B =  9.0902171e+02,
      C = -7.3129061e+05,
      D = -5.3503838e-01,
   },
}
db['N2'].thermal_conductivity = {
   model = 'CEA',
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      A =  8.5439436e-01,
      B =  1.0573224e+02,
      C = -1.2347848e+04,
      D =  4.7793128e-01,
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  8.8407146e-01,
      B =  1.3357293e+02,
      C = -1.1429640e+04,
      D =  2.4417019e-01,
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  2.4176185e+00,
      B =  8.0477749e+03,
      C =  3.1055802e+06,
      D = -1.4517761e+01,
   },
}
db['O2'] = {}
db['O2'].atomicConstituents = { O=2, }
db['O2'].charge = 0
db['O2'].M = 3.19988000e-02
db['O2'].sigma = 3.45800000
db['O2'].epsilon = 107.40000000
db['O2'].Lewis = 1.08600000
db['O2'].thermoCoeffs = {
  origin = 'CEA',
  nsegments = 3, 
  T_break_points = { 200.00, 1000.00, 6000.00, 20000.00, },
  T_blend_ranges = { 400.0, 1000.0, },
  segment0 = {
   -3.425563420e+04,
    4.847000970e+02,
    1.119010961e+00,
    4.293889240e-03,
   -6.836300520e-07,
   -2.023372700e-09,
    1.039040018e-12,
   -3.391454870e+03,
    1.849699470e+01,
  },
  segment1 = {
   -1.037939022e+06,
    2.344830282e+03,
    1.819732036e+00,
    1.267847582e-03,
   -2.188067988e-07,
    2.053719572e-11,
   -8.193467050e-16,
   -1.689010929e+04,
    1.738716506e+01,
  },
  segment2 = {
    4.975294300e+08,
   -2.866106874e+05,
    6.690352250e+01,
   -6.169959020e-03,
    3.016396027e-07,
   -7.421416600e-12,
    7.278175770e-17,
    2.293554027e+06,
   -5.530621610e+02,
  },
}
db['O2'].viscosity = {
   model = 'CEA',
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      A =  6.0916180e-01,
      B = -5.2244847e+01,
      C = -5.9974009e+02,
      D =  2.0410801e+00,
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  7.2216486e-01,
      B =  1.7550839e+02,
      C = -5.7974816e+04,
      D =  1.0901044e+00,
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A =  7.3981127e-01,
      B =  3.9194906e+02,
      C = -3.7833168e+05,
      D =  9.0931780e-01,
   },
}
db['O2'].thermal_conductivity = {
   model = 'CEA',
   nsegments = 3,
   segment0 = {
      T_lower = 200.0,
      T_upper = 1000.0,
      A =  7.7229167e-01,
      B =  6.8463210e+00,
      C = -5.8933377e+03,
      D =  1.2210365e+00,
   },
   segment1 = {
      T_lower = 1000.0,
      T_upper = 5000.0,
      A =  9.0917351e-01,
      B =  2.9124182e+02,
      C = -7.9650171e+04,
      D =  6.4851631e-02,
   },
   segment2 = {
      T_lower = 5000.0,
      T_upper = 15000.0,
      A = -1.1218262e+00,
      B = -1.9286378e+04,
      C =  2.3295011e+07,
      D =  2.0342043e+01,
   },
}