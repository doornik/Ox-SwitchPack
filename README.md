# Ox-SwitchingAcceleration
Ox code to illustrate acceleration of EM and other switching algorithms.

## Summary

This code supplements 	J.A. Doornik (2017 or 2018), `Accelerated Estimation of Switching
Algorithms: The Cointegrated VAR Model and Other Applications', forthcoming,
*Scandinavian Journal of Statistics*.

Please reference the above paper when using this code.

## Contents

1. Data
  - `Danish_data_Juselius(2006).xlsx`: the Danish data of Juselius (2006).
2. Supporting Ox code
  - `switching.ox, switching.oxh`: Switching class for maximization with a
  choice of line searches.
  - `maxscalar.ox, maxscalar.oxh`: univariate maximization algorithms by Brent and Powell.
  - `coint1.ox, coint1.oxh`: Coint class for simple I(1) CVAR estimation.
3. Ox programs
  - `switching_coint.ox`: Ox program to evaluate the switching algorithms for simplified version of
  the I(1) model using the Danish data.
  - `switching_gxm.ox`: Ox program to evaluate the EM algorithms for the Gaussian mixture model.
  - `switching_tests.ox`: Ox program to evaluate the EM algorithms for the Poisson mixture
  model and the multivariate-t, as well as alternating algorithms for PARAFAC and low-rank matrix approximation.
  - `test_maxscalar.ox`: Ox program to run some examples using maxscalar.  
  - `trace_plots.ox`: Ox program to make evaluation plots of the more extensive output
  that is generated in the trace folder when switching_coint.ox is run with trace switched on.
4. Output
- the `output' folder has the results of some of my runs of the code.
