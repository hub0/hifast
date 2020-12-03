# Change Log

## 20201104

- Correct the doc string of coord.load_ky

## 20201127

- Rename chart.pols_trim to chart.xyyx_trim

## 20201201

- Change mean() in calibration.smooth to more explicit np.mean()

## 20201202

Update:

- Chart.freq changes from 2D array to 1D array
  
  The frequency setup in Chart is now uniform for all spectra.

- Add a frequency smooth function in hifast.calibration

TODO:

- Move frequency smooth function to Chart method
