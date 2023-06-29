# brca1-bilayer-duct-model
Math model of breast duct tissue with chemomechanical feedback

This is an ODE model of bilayered epithelial duct tissue (basal and luminal layers) in mammary glands which incorporates basal-luminal intermediate (BLI) cells.
See a project summary as of June 2023 at: https://www.dropbox.com/s/ik1mirf5nimtm8i/brca1-breast-cancerSP2023.pdf?dl=0

### Key files
  * Most recent ODEs are implemented in `bilayerDuct_05.m` (function file)
  * Run `figuresForSlides.m` to see parameter values, example numerical solutions, plots (corresponding to those in the linked Dropbox slides)
  * Named terms in the ODEs (e.g., proliferation rates, activin level) can be computed by passing the numerical solution and parameters to `usefulQuantities_05.m` (function file, see `premaligDetails.m` for example usage)
