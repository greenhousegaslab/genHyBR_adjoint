# genHyBR_adjoint
The repository includes code to run the genHyBR inverse modeling algorithm with an adjoint model like the GEOS-Chem adjoint model.

**Overview**
Generalized hybrid methods (genHyBR) are computationally efficient algorithms for solving very large inverse problems. Specifically, these algorithms do not require inverting any large, dense matrices and avoid intensive matrix-matrix multiplication. These algorithms can be used to solve inverse problems with millions of observations and millions of unknown emissions or fluxes. In addition, geHyBR can be used to estimate fluxes or emissions that vary both spatially and temporally.

Generalized hybrid methods are iterative. The algorithm here will first transform the inverse problem to get rid of large matrix inversions. Subsequently, the algorithm will compute basis vectors that define a Krylov subspace. These basis vectors are calculated using an iterative algorithm known as a Golub–Kahan bidiagonalization process. These vectors can then be used to estimate the fluxes/emissions and the posterior uncertainties in those fluxes/emissions. The greater the number of iterations, the more vectors that define the Krylov subspace, and the better the flux estimate and posterior uncertainty estimate.

Julianne Chung created a Github repository with scripts for running genHyBR (https://github.com/juliannechung/genHyBR). Those scripts are difficult to implement with many chemical transport models like the GEOS-Chem model, and the codes in this repository are easier to implement with thesee chemical transport models. Specifically, at each iteration, Julianne's code first calls the adjoint of the forward model and then calls the forward model. Models like GEOS-Chem by contrast, will first run the forward model and then run the adjoint (the reverse order of Julianne's code). The scripts attached here have been reformulated to first run the forward and then the adjoint model at each iteration. These scripts are therefore easier to implement with models like GEOS-Chem.

**Implementation**
These scripts require the user to write two Matlab functions. The function afun.m should (1) Write a set of emissions to a file format that can be read by the forward/adjoint model of your choice, (2) should launch the forward and adjoint models, and (3) read in the results of those forward/adjoint model runs. This function is called at every iteration of the genHyBR algorithm, except the first iteration. The first iteration of geHyBR requires a function called afunadj.m. That function should (1) write a set of observations to a file format that can be read by the adjoint model of your choice, (2) run that adjoint model, and (3) read in the results of that adjoint model run. No forward model run is required for the function afunadj.m.

Leyang Feng (lfeng13@jhu.edu) has written "afun" and "afunadj" functions for the GEOS-Chem methane adjoint models. These scripts from Leyang could be used as a template for creating your own "afun" and "afunadj" functions for GEOS-Chem.

**References**
Chung, J., & Saibaba, A. K. (2017). Generalized hybrid iterative methods for large-scale Bayesian inverse problems. SIAM Journal on Scientific Computing, 39(5), S24-S46, doi:10.1137/16M1081968.

Cho, T., Chung, J., Miller, S. M., and Saibaba, A. K. (2022). Computationally efficient methods for large-scale atmospheric inverse modeling, Geosci. Model Dev., 15, 5547–5565, https://doi.org/10.5194/gmd-15-5547-2022. (Note: There is a mistake in the eqs. for u and v in Sect. 3.1.1 of this manuscript. Refer to Chung and Saibaba (above) for the correct version of these equations.)
