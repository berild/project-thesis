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
  <summary id="blm" style ="cursor: pointer; font-size: 1.5em;">Bivariate Linear</summary>
  
To apply the combined methods on the bivariate linear model, run the <a href="https://github.com/berild/project-thesis/blob/master/linreg/linreg.R">linreg.R</a> script. The functions for each algorithm is presented in <a href="https://github.com/berild/project-thesis/blob/master/linreg/linreg_amis_w_inla.R">linreg_amis_w_inla.R</a>, <a href="https://github.com/berild/project-thesis/blob/master/linreg/linreg_is_w_inla.R">linreg_is_w_inla.R</a>, and <a href="https://github.com/berild/project-thesis/blob/master/linreg/linreg_mcmc_w_inla.R">linreg_mcmc_w_inla.R</a>. General functions used in the algorithms, plotting, and evaluation is given in <a href="https://github.com/berild/project-thesis/blob/master/linreg/linreg_general_functions.R">linreg_general_functions.R</a>. Result of our simulation is given with <a href="https://git-lfs.github.com/"><img src="https://avatars0.githubusercontent.com/u/20246716?s=280&v=4" alt="git-lfs logo" style="width: 1.5em;" />git-lfs</a> in <a href="https://github.com/berild/project-thesis/tree/master/linreg/sims">sims</a> and use <a href="https://github.com/berild/project-thesis/blob/master/linreg/plot_linreg.R">plot_linreg.R</a> to replicate our plots.
  
  <h4>Model</h4>
  
  In this example, we simulate our own data using a linear model with the Gaussian noise term $\epsilon_i \sim \mathcal{N}(0,1/\tau)$, with zero mean and precision $\tau$. 
  The response $y_i$ is modelled by two covariates $x_{1i}$ and $x_{2i}$ generated from a uniform distribution between zero and one, i.e. $x_{1i},x_{2i} \sim \mathcal{U}(0,1)$.
  We generate 100 such samples and calculate the response by 
  $$
  y_i = \alpha + \beta_1x_{1i} + \beta_2x_{2i} + \epsilon_i,
  $$
  using $\alpha = 3$, $\beta_1 = 2$, $\beta_2 = -2$, and $\tau = 1$. 
  
  <h4>Result</h4>
  
<img src="https://i.imgur.com/592pwcu.png"
     alt="Bivariate Linear Regression"
     style="width: 70%; display: block; margin-left: auto; margin-right: auto;" /> 
</details>

<details>
  <summary id="bl" style ="cursor: pointer; font-size: 1.5em;">Bayesian Lasso</summary>
  
To apply the combined methods on the Bayesian lasso model, run the <a href="https://github.com/berild/project-thesis/blob/master/lasso/lasso.R">lasso.R</a> script. The functions for each algorithm is presented in <a href="https://github.com/berild/project-thesis/blob/master/lasso/lasso_amis_w_inla.R">lasso_amis_w_inla.R</a>, <a href="https://github.com/berild/project-thesis/blob/master/lasso/lasso_is_w_inla.R">lasso_is_w_inla.R</a>, and <a href="https://github.com/berild/project-thesis/blob/master/lasso/lasso_mcmc_w_inla.R">lasso_mcmc_w_inla.R</a>. General functions used in the algorithms, plotting, and evaluation is given in <a href="https://github.com/berild/project-thesis/blob/master/lasso/lasso_general_functions.R">lasso_general_functions.R</a>. Result of our simulation is given with <a href="https://git-lfs.github.com/"><img src="https://avatars0.githubusercontent.com/u/20246716?s=280&v=4" alt="git-lfs logo" style="width: 1.5em;" />git-lfs</a> in <a href="https://github.com/berild/project-thesis/tree/master/lasso/sims">sims</a> and use <a href="https://github.com/berild/project-thesis/blob/master/lasso/plot_lasso.R">plot_lasso.R</a> to replicate our plots.
  
  <h4>Model</h4>
  
  The data used in this example is found in the Hitters dataset in the R package ISLR. It holds Major League Baseball data from 1986 and 1987. 
  We want to model the 1987 Salary based on AtBat, Hits, HmRun, Runs, and RBI from the 1986 season. 
  Lasso regression is achieved by solving the minimization problem 
  $$
  \sum_{i=1}^N \left(y_i - \alpha - \sum_{j=1}^{n_\beta} \beta_j x_{ji}\right)^2 + \lambda\sum_{j=1}^{n_\beta}|\beta_j|,
  $$
  where $\lambda$ controls the shrinkage of the $\beta_j$ parameters. We what to achieve inference about the parameters in the Bayesian Lasso model.
  
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
  <summary id="mc" style ="cursor: pointer; font-size: 1.5em;">Missing Covariates</summary>
    
 To apply the combined methods on the linear model with missing covariates, run the <a href="https://github.com/berild/project-thesis/blob/master/missing/missing.R">lasso.R</a> script. The functions for each algorithm is presented in <a href="https://github.com/berild/project-thesis/blob/master/missing/missing_amis_w_inla.R">missing_amis_w_inla.R</a>, <a href="https://github.com/berild/project-thesis/blob/master/missing/missing_is_w_inla.R">missing_is_w_inla.R</a>, and <a href="https://github.com/berild/project-thesis/blob/master/missing/missing_mcmc_w_inla.R">missing_mcmc_w_inla.R</a>. General functions used in the algorithms, plotting, and evaluation is given in <a href="https://github.com/berild/project-thesis/blob/master/missing/missing_general_functions.R">missing_general_functions.R</a>. Result of our simulation is given with <a href="https://git-lfs.github.com/"><img src="https://avatars0.githubusercontent.com/u/20246716?s=280&v=4" alt="git-lfs logo" style="width: 1.5em;" />git-lfs</a> in <a href="https://github.com/berild/project-thesis/tree/master/missing/sims">sims</a> and use <a href="https://github.com/berild/project-thesis/blob/master/missing/plot_missing.R">plot_missing.R</a> to replicate our plots.   
    
  <h4>Model</h4>
  
  In this example, we will apply the combined methods on the nhanes dataset provided by the R package mice. We will model the cholesterol level by BMI and age using a linear model with Gaussian noise term:
  $$
  y_i = \beta_0 + \beta_1 x_{1i} + \beta_2 x_{2i} + \beta_3 x_{3i} + \epsilon_i.
  $$
  The BMI contains missing values which we impute with the Monte Carlo method, in addition, a predictive distribution of the missing values in cholesterol level is provide by INLA. 
  
  <h4>Result</h4>
  | Parameter |  MCMC w/INLA  |   IS w/INLA   |  AMIS w/INLA  |
|:---------:|:-------------:|:-------------:|:-------------:|
| $\beta_0$ |  43.69(62.20) |  44.21(62.27) |  42.87(62.02  |
| $\beta_1$ |   4.86(2.19)  |   4.84(2.19)  |   4.89(2.19)  |
| $\beta_2$ |  29.37(17.87) |  29.40(17.90) |  29.49(17.84) |
| $\beta_3$ |  49.61(23.12) |  49.49(23.24) |  49.83(23.14) |
|   $\tau$  | 0.001(0.0005) | 0.001(0.0005) | 0.001(0.0005) |
  
<img src="https://i.imgur.com/R1CtyWh.png"
     alt="Missing Covariates Parameter Marginals"
     style="width: 70%; display: block; margin-left: auto; margin-right: auto;" /> 
<img src="https://i.imgur.com/UBF1GEf.png"
     alt="Missing Covariates Observation Marginals"
     style="width: 70%; display: block; margin-left: auto; margin-right: auto;" /> 
</details>


<details>
  <summary id="sem" style ="cursor: pointer; font-size: 1.5em;">Spatial Econometric Model</summary>
  
   To apply the combined methods on the spatial lag model, run the <a href="https://github.com/berild/project-thesis/blob/master/sem/sem.R">sem.R</a> script. The functions for each algorithm is presented in <a href="https://github.com/berild/project-thesis/blob/master/sem/sem_amis_w_inla.R">sem_amis_w_inla.R</a>, <a href="https://github.com/berild/project-thesis/blob/master/sem/sem_is_w_inla.R">sem_is_w_inla.R</a>, and <a href="https://github.com/berild/project-thesis/blob/master/sem/sem_mcmc_w_inla.R">sem_mcmc_w_inla.R</a>. General functions used in the algorithms, plotting, and evaluation is given in <a href="https://github.com/berild/project-thesis/blob/master/sem/sem_general_functions.R">sem_general_functions.R</a>. Result of our simulation is given with <a href="https://git-lfs.github.com/"><img src="https://avatars0.githubusercontent.com/u/20246716?s=280&v=4" alt="git-lfs logo" style="width: 1.5em;" />git-lfs</a> in <a href="https://github.com/berild/project-thesis/tree/master/sem/sims">sims</a> and use <a href="https://github.com/berild/project-thesis/blob/master/sem/plot_sem.R">plot_sem.R</a> to replicate our plots.   
    
  <h4>Model</h4>
  In this example, we will apply the combined methods on a spatial lag model. Furthermore, we will consider Columbus dataset provided by the R package spData, and we will model the crime rates in Columbus based on the houshold income (INC) and the houseing value (HOVAL). 
  The spatial lag model is expressed as
  $$
  \mathbf{y} = \rho\mathbf{W}\mathbf{y} + \mathbf{X}\beta + \mathbf{u},
  $$
  with the Gaussian noise $\mathbf{u} \sim \mathcal{N}(0,\tau)$. This expression is generally rewritten as
  $$
  \mathbf{y} = (\mathbf{I}_n - \rho\mathbf{W})^{-1}\mathbf{X}\beta + \epsilon,
  $$
  where the noise term is given by
  $$
  \epsilon \sim \mathcal{N}\left(0,\frac{1}{\tau}\left[(\mathbf{I}_n - \rho\mathbf{W}^T)(\mathbf{I}_n - \rho\mathbf{W})\right]^{-1}\right).
  $$
  
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
* We will create a R package for the methodologies which will be linked here later.
* We will apply the AMIS with INLA method in a master thesis and the repository will be linked here later.
* We will link to the project thesis paper, when the grading is completed.

<h3 id="references">References</h3>

<a href="https://link.springer.com/article/10.1007/s11222-017-9778-y" id ="mcmcwinla">Gomez-Rubio, V., & Rue, H. (2018). Markov chain Monte Carlo with the Integrated Nested Laplace Approximation. *Statistics and Computing, 28*(5), 1033.</a>


<a href="https://www.jstor.org/stable/40247579?seq=1" id ="inla">Rue, H., Martino, S., & Chopin, N. (2009). Approximate Bayesian inference for latent Gaussian models by using integrated nested Laplace approximations. *Journal of the Royal Statistical Society: Series B (Statistical Methodology), 71*(2), 319-392.</a>


<a id ="iswinla" href="">Gomez-Rubio, V. (2019). Importance Sampling with the Integrated Nested Laplace Aproximation. *BISP 2019 conference,Madrid, Spain, 12-14 June*, Poster</a>

<a href = "https://www.jstor.org/stable/23357226?seq=1" id = "amis">Cornuet, J., Marin, J., Mira, A., & Robert, C. (2012). Adaptive Multiple Importance Sampling. *Scandinavian Journal of Statistics, 39*(4), 798-812. </a>