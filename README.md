# Project-Instantaneous-Phase-Estimator

The goal of this project is to estimate the phase quickly and causally. See <a href="https://github.com/UriEdenLab/Project-Instantaneous-Phase-Estimator/wiki/Overview-and-Goals">wiki</a> for more detailed descriptions.

The file *particle_filter_phase_tracker.m* is the "Eden method". It uses the example data *sample_data.mat*.

Added code to compute CFC.  Here's an example:

```
Vlo =                 %data filtered to low frequency band.
Vhi =                 %data filtered to high frequency band.
nCtlPts = 8;          %A reasonable default choice.
[r, r_CI,r2,r2_CI, nCtlPts] = GLM_PAC(Vlo, Vhi, nCtlPts);
                      %Examine the resulting figure, and inspect the values of r, r_CI
```
