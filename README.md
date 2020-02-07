<h2>The Integrated Nested Laplace Approximation within Monte Carlo methods</h2>

------

* Project Thesis delivered January 20, 2020
* Project Title: Adaptive Multiple Importance Sampling with the Integrated Nested Laplace Approximation
* Author: Martin Outzen Berild

------

This repository contains the code used in the project thesis. The different methodologies was applied to four examples; biviariate linear model, Bayesian Lasso model, linear model with missing covariates, and Spatial Econometric Model. The examples are equivalelent to the ones presented in <a href="#mcmcwinla">Gomez-Rubio and Rue (2018)</a>.

<h3 id="abstract">Abstract</h3>

The Integrated Nested Laplace Approximation (INLA) is a deterministic approach for Bayesian inference on latent Gaussian models (LGMs) introduced by <a href="#inla">Rue et al. (2009)</a>. 
INLA focuses on approximating fast and accurate posterior marginals of the parameters in the model. 
We have explored existing methods that combine INLA and Monte Carlo simulation to extend the set of models that can be fitted with INLA to those that can be expressed as conditional LGMs. 
More specifically, we describe the INLA within Metropolis-Hastings algorithm introduced by 
<a href="#mcmcwinla">Gomez-Rubio and Rue (2018)</a> and the importance sampling (IS) with INLA approach proposed by <a href="#iswinla">Gomez-Rubio (2019)</a>.
In light of these methods and the recent development of more robust importance sampling schemes, we present a novel approach that combines INLA and the adaptive multiple importance sampling (AMIS; <a href="#amis">Cornuet et al. (2012)</a>) algorithm. 

The aim is to assess these methods on a series of applications with simulated and existing datasets, and compare their performance based on accuracy and efficiency.

<h3 id="examples">Examples</h3>

<details>
  <summary id="blm" style ="cursor: pointer; font-size: 1.5em;">Bivariate Linear (click to view)</summary>
  
To apply the combined methods on the bivariate linear model, run the <a href="https://github.com/berild/project-thesis/blob/master/linreg/linreg.R">linreg.R</a> script. The functions for each algorithm is presented in <a href="https://github.com/berild/project-thesis/blob/master/linreg/linreg_amis_w_inla.R">linreg_amis_w_inla.R</a>, <a href="https://github.com/berild/project-thesis/blob/master/linreg/linreg_is_w_inla.R">linreg_is_w_inla.R</a>, and <a href="https://github.com/berild/project-thesis/blob/master/linreg/linreg_mcmc_w_inla.R">linreg_mcmc_w_inla.R</a>. General functions used in the algorithms, plotting, and evaluation is given in <a href="https://github.com/berild/project-thesis/blob/master/linreg/linreg_general_functions.R">linreg_general_functions.R</a>. Result of our simulation is given with <a href="https://git-lfs.github.com/">git-lfs</a> in <a href="https://github.com/berild/project-thesis/tree/master/linreg/sims">sims</a> and use <a href="https://github.com/berild/project-thesis/blob/master/linreg/plot_linreg.R">plot_linreg.R</a> to replicate our plots.

  <h4>Result</h4>
  
<img src="https://i.imgur.com/592pwcu.png"
     alt="Bivariate Linear Regression"
     style="width: 70%; display: block; margin-left: auto; margin-right: auto;" /> 
</details>

<details>
  <summary id="bl" style ="cursor: pointer; font-size: 1.5em;">Bayesian Lasso (click to view)</summary>
  
To apply the combined methods on the Bayesian lasso model, run the <a href="https://github.com/berild/project-thesis/blob/master/lasso/lasso.R">lasso.R</a> script. The functions for each algorithm is presented in <a href="https://github.com/berild/project-thesis/blob/master/lasso/lasso_amis_w_inla.R">lasso_amis_w_inla.R</a>, <a href="https://github.com/berild/project-thesis/blob/master/lasso/lasso_is_w_inla.R">lasso_is_w_inla.R</a>, and <a href="https://github.com/berild/project-thesis/blob/master/lasso/lasso_mcmc_w_inla.R">lasso_mcmc_w_inla.R</a>. General functions used in the algorithms, plotting, and evaluation is given in <a href="https://github.com/berild/project-thesis/blob/master/lasso/lasso_general_functions.R">lasso_general_functions.R</a>. Result of our simulation is given with <a href="https://git-lfs.github.com/">git-lfs</a> in <a href="https://github.com/berild/project-thesis/tree/master/lasso/sims">sims</a> and use <a href="https://github.com/berild/project-thesis/blob/master/lasso/plot_lasso.R">plot_lasso.R</a> to replicate our plots.
  
  <h4>Result</h4>
  
  
  | Parameter | Lasso |  MCMC w/INLA  |   IS w/INLA   | AMIS w/INLA   |
|:---------:|:-----:|:-------------:|:-------------:|---------------|
|   AtBat   | 0.000 | -0.007(0.006) | -0.088(0.002) | -0.008(0.007) |
|    Hits   | 0.182 |  0.166(0.013) |  0.260(0.004) | 0.165(0.013)  |
|   HmRun   | 0.000 |  0.026(0.004) |  0.082(0.001) | 0.024(0.004)  |
|    Runs   | 0.000 |  0.074(0.007) |  0.099(0.001) | 0.072(0.007)  |
|    RBI    | 0.232 |  0.204(0.012) |  0.015(0.002) | 0.209(0.012)  |
  
<img src="https://i.imgur.com/P35eM7f.png"
     alt="Bayesian Lasso Regression"
     style="width: 70%; display: block; margin-left: auto; margin-right: auto;" /> 
</details>


<details>
  <summary id="mc" style ="cursor: pointer; font-size: 1.5em;">Missing Covariates (click to view)</summary>
    
 To apply the combined methods on the linear model with missing covariates, run the <a href="https://github.com/berild/project-thesis/blob/master/missing/missing.R">lasso.R</a> script. The functions for each algorithm is presented in <a href="https://github.com/berild/project-thesis/blob/master/missing/missing_amis_w_inla.R">missing_amis_w_inla.R</a>, <a href="https://github.com/berild/project-thesis/blob/master/missing/missing_is_w_inla.R">missing_is_w_inla.R</a>, and <a href="https://github.com/berild/project-thesis/blob/master/missing/missing_mcmc_w_inla.R">missing_mcmc_w_inla.R</a>. General functions used in the algorithms, plotting, and evaluation is given in <a href="https://github.com/berild/project-thesis/blob/master/missing/missing_general_functions.R">missing_general_functions.R</a>. Result of our simulation is given with <a href="https://git-lfs.github.com/">git-lfs</a> in <a href="https://github.com/berild/project-thesis/tree/master/missing/sims">sims</a> and use <a href="https://github.com/berild/project-thesis/blob/master/missing/plot_missing.R">plot_missing.R</a> to replicate our plots.   
  
  <h4>Result</h4>
  
  
| Parameter |  MCMC w/INLA  |   IS w/INLA   |  AMIS w/INLA  |
|:---------:|:-------------:|:-------------:|:-------------:|
|  &beta;0  |  43.69(62.20) |  44.21(62.27) |  42.87(62.02  |
|  &beta;1  |   4.86(2.19)  |   4.84(2.19)  |   4.89(2.19)  |
|  &beta;2  |  29.37(17.87) |  29.40(17.90) |  29.49(17.84) |
|  &beta;3  |  49.61(23.12) |  49.49(23.24) |  49.83(23.14) |
|   &tau;   | 0.001(0.0005) | 0.001(0.0005) | 0.001(0.0005) |
  
<img src="https://i.imgur.com/R1CtyWh.png"
     alt="Missing Covariates Parameter Marginals"
     style="width: 70%; display: block; margin-left: auto; margin-right: auto;" /> 
<img src="https://i.imgur.com/UBF1GEf.png"
     alt="Missing Covariates Observation Marginals"
     style="width: 70%; display: block; margin-left: auto; margin-right: auto;" /> 
</details>


<details>
  <summary id="sem" style ="cursor: pointer; font-size: 1.5em;">Spatial Econometric Model (click to view)</summary>
  
   To apply the combined methods on the spatial lag model, run the <a href="https://github.com/berild/project-thesis/blob/master/sem/sem.R">sem.R</a> script. The functions for each algorithm is presented in <a href="https://github.com/berild/project-thesis/blob/master/sem/sem_amis_w_inla.R">sem_amis_w_inla.R</a>, <a href="https://github.com/berild/project-thesis/blob/master/sem/sem_is_w_inla.R">sem_is_w_inla.R</a>, and <a href="https://github.com/berild/project-thesis/blob/master/sem/sem_mcmc_w_inla.R">sem_mcmc_w_inla.R</a>. General functions used in the algorithms, plotting, and evaluation is given in <a href="https://github.com/berild/project-thesis/blob/master/sem/sem_general_functions.R">sem_general_functions.R</a>. Result of our simulation is given with <a href="https://git-lfs.github.com/">git-lfs</a> in <a href="https://github.com/berild/project-thesis/tree/master/sem/sims">sims</a> and use <a href="https://github.com/berild/project-thesis/blob/master/sem/plot_sem.R">plot_sem.R</a> to replicate our plots.   

  <h4>Result</h4>

| Parameter |  Max.Like.  | MCMC w/INLA |  IS w/INLA  | AMIS w/INLA |
|:---------:|:-----------:|:-----------:|:-----------:|:-----------:|
| Intercept | 61.05(5.31) | 60.65(6.23) | 60.62(6.23) | 60.66(6.23) |
|    INC    | -1.00(0.34) | -0.97(0.38) | -0.97(0.37) | -0.96(0.37) |
|   HOVAL   | -0.31(0.09) | -0.31(0.09) | -0.31(0.09) | -0.30(0.09) |
|   &rho;   |  0.52(0.14) |  0.55(0.14) |  0.55(0.13) |  0.55(0.13) |
|   &tau;   |   0.01(-)   | 0.01(0.002) | 0.01(0.002) | 0.01(0.001) |
 
<img src="https://i.imgur.com/2QK6jUA.png"
     alt="Spatial Econometric Model"
     style="width: 70%; display: block; margin-left: auto; margin-right: auto;" /> 
</details>

<h3>Further Work</h3>

<ul>
  <li>We will create a R package for the methodologies which will be linked here later.</li>
  <li>We will apply the AMIS with INLA method in a master thesis and the repository will be linked here later.</li>
  <li>We will link to the project thesis paper, when the grading is completed.</li>
</ul>


<h3 id="references">References</h3>

<a href="https://link.springer.com/article/10.1007/s11222-017-9778-y" id ="mcmcwinla">Gomez-Rubio, V., & Rue, H. (2018). Markov chain Monte Carlo with the Integrated Nested Laplace Approximation. *Statistics and Computing, 28*(5), 1033.</a>


<a href="https://www.jstor.org/stable/40247579?seq=1" id ="inla">Rue, H., Martino, S., & Chopin, N. (2009). Approximate Bayesian inference for latent Gaussian models by using integrated nested Laplace approximations. *Journal of the Royal Statistical Society: Series B (Statistical Methodology), 71*(2), 319-392.</a>


<a id ="iswinla" href="">Gomez-Rubio, V. (2019). Importance Sampling with the Integrated Nested Laplace Aproximation. *BISP 2019 conference,Madrid, Spain, 12-14 June*, Poster</a>

<a href = "https://www.jstor.org/stable/23357226?seq=1" id = "amis">Cornuet, J., Marin, J., Mira, A., & Robert, C. (2012). Adaptive Multiple Importance Sampling. *Scandinavian Journal of Statistics, 39*(4), 798-812. </a>