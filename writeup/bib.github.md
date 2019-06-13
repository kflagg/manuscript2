Illian, Sørbye, and Rue (2012) \*
=================================

-   Presents a methodological "toolbox" for using INLA to fit spatial marked Markov processes driven by latent Gaussian process intensity functions
-   Discretize window into cells and model cell counts via Poisson GLMM with log link
-   Linear predictor includes *p* spatial functions and a nonlinear function of a constructed covariate
    -   Spatial functions modeled as random walks constrained to sum to zero
    -   Constructed covariate is a summary describing small-scale interaction among points
    -   Constructed covariate allows interaction-like behavior in LGCP, easier to fit than Gibbs process because interaction acts directly on intensity
    -   Relationship between constructed covariate and intensity modeled as random walk constrained to sum to zero
    -   Constructed covariate used in the paper is distance from grid cell center to nearest point outside grid cell
-   Also includes a residual term which can be used to assess model adequacy
-   Choice of prior for spatial effect is tied to choice of grid size
-   Suggest a grid "not finer than that given by the data"
-   Model fit to a simulated repulsive Strauss process and posterior predictive check of *L*-function appeared reasonable
-   Model fit to a simulated Thomas cluster process did not have strong enough clustering
-   Example with many small clusters and a spatial trend in the intensity of cluster centers described the trend well but underestimated the strength of clustering
-   Fit to homogeneous Poisson process does not include spurious spatial effects
-   Discussion of other possible constructed covariates and difficulties of identifying different effects at same scale
-   Can include covariates observed at different locations that the points being modeled
-   Covariate term in linear predictor for intensity is prediction from spatial model for covariate
-   Also illustrate model for marks with spatial GLMM for marks

Simpson et al. (2016) \*
========================

-   Includes a clear look at the math of the likelihood and approximation
-   Fine grid method of Illian, Sørbye, and Rue (2012) is accurate but computationally wasteful because it induces a dense covariance matrix
    -   CAR alleviates this but requires regular lattice
    -   Lattice needs to be fine enough to accurately represent the locations of points which is finer than needed to accurately represent the latent spatial field
    -   Keeping constant resolution of the lattice (espicially in unsurveyed regions) is particularly wasteful
-   Better approximation represents spatial field as a linear combination of a finite number of basis functions
-   Explains the **A** matrix
-   Ends up using a Baddeley and Turner style approximation to the likelihood
    -   Looks like independent Poisson random variables at point locations (with mean 1) and at numerical integration nodes (with mean equal to the appoximated intensity at the node)
    -   Requires a good numerical integration scheme

Banerjee et al. (2008) \*
=========================

-   Simplify high-dimensional spatial and spatiotemporal predictions by projecting the spatial process into a lower-dimensional subset observed at "knots"
-   Smaller covariance matrix to invert
-   MCMC requires "reasonably informative priors"
-   Useful with tens of thousands of prediction points, needs hundreds of knots
-   "\[F\]aster alternatives that avoid MCMC sampling can be employed" (cites the INLA paper)

Andrieu, Doucet, and Holenstein (2010) \*
=========================================

-   PMCMC combines SMC and MCMC to handle high-dimensional models
-   SMC used to generate proposals
-   Focus on state space models
-   Particle degeneracy less of a problem than with pure SMC
-   Efficient exploration wthout a lot of MCMC tuning but slower than INLA

Girolami and Calderhead (2011)
==============================

-   Exploits Riemann manifold structure of parameter space for efficient MC sampling without tuning in high-dimensional situations
-   LGCP example, notes that INLA could speed up computation
-   Langevin diffusion/Metropolis adjusted Langevin algorithm(MALA)?
-   Extend HMC and MALA, stochastic differential equations are involved

Marin et al. (2012) \*
======================

-   Review paper

Murray, Adams, and MacKay (2010) \*
===================================

-   MCMC for models with multivariate Gaussian priors
-   Slower but more general than approximations
-   No tuning parameters
-   Better than Gibbs sampling when latent variables are correlated
-   Sort of a modified MH with proposals on an ellipse that passes through the current state

Bagchi et al. (2014)
====================

-   Experiment applying insecticide and fungicide in a rainforest to see the effects on plant development
-   Weekly seed and seedling counts over 17 weeks
-   Hierarchical model with negative binomial response and latent normal variables fit in R-INLA
-   Priors not stated

Diggle, Menezes, and Su (2010) \*
=================================

-   Geostatistical process model expanded by conditioning on a point process model for the observation locations
-   Needed when the process of selecting observation locations depends on the response
-   Describes ML inference by conditional simulation
-   Discussion response by Rue, Martino, Mondal, and Chopin mentions INLA as applicable

Schrödle and Held (2011)
========================

-   Use INLA to fit hierarchical model with Poisson observations and a latent spatiotemporal Gaussian process to annual region-level counts of cattle salmonellosis in Switzerland
-   Describe how to fit an Intrinsic Gaussian Markov random field model (with linear constraints to ensure parameter identifiability) in INLA

Miller et al. (2013)
====================

-   Discusses mapping dolphin population density using distance sampling along transects and a GAM for Poisson process intensity
-   INLA mentioned as a way to quickly fit Bayesian models

Wang and Blei (2013) \*
=======================

-   Discusses Laplace variational inference and delta method variational inference as generic approaches for non-conjugate exponential family models
-   Describes Laplace approximation as Taylor expansion around MAP point, resulting in Gaussian approximation
-   Conclude that Laplace approximation is more accurate than delta method

Cameletti et al. (2013)
=======================

-   Discusses Gaussian Markov random field as lower-dimensional alternative to Gaussian process with dense covariance for spatiotemporal prediction
-   Can be computed quickly using INLA and stochastic PDEs
-   SPDE approach uses finite element approximation (linear combination of basis functions) of spatial field on triangular grid
-   Illustration of particulate matter modeling fit in R-INLA

Martins et al. (2013)
=====================

-   Presents improvements implemented in R-INLA in the four years after the original INLA paper
-   Reviews a variety of applications where INLA had been previously shown to do well
-   Summarizes INLA
-   Adds some flexibility to allow more complicated dependence structures

Noor et al. (2014)
==================

-   Spatialtemporal model used to map prevalence of a malaria-spreading parasite across Africa
-   Fit using R-INLA and the SPDE approach (Cameletti et al. 2013)
-   Possibly questionable data collection and variable selection

Fong, Rue, and Wakefield (2010)
===============================

-   Compares PQL and INLA for several biostatistical examples
-   Points out that INLA approximation can be inaccurate for binomial GLMMs with small denominators

Scott et al. (2016) \*
======================

-   Setting where datasets are so large that communication within a computing cluster is cost-prohibitive
-   "Consensus Monte Carlo operates by running a separate Monte Carlo algorithm on each machine, and then averaging individual Monte Carlo draws across machines."
-   Authors suggest this is useful when the parameter space is too high-dimensional for INLA to be efficient or when joint posteriors are needed
-   Data partitioned into "shards" and MC samplers run independently on each shard
-   "Consensus posterior" constructed from draws which are weighted averages of draws from each sampler

References
==========

Andrieu, Christophe, Arnaud Doucet, and Roman Holenstein. 2010. “Particle Markov Chain Monte Carlo Methods.” *Journal of the Royal Statistical Society: Series B (Statistical Methodology)* 72 (3): 269–342.

Bagchi, Robert, Rachel E Gallery, Sofia Gripenberg, Sarah J Gurr, Lakshmi Narayan, Claire E Addis, Robert P Freckleton, and Owen T Lewis. 2014. “Pathogens and Insect Herbivores Drive Rainforest Plant Diversity and Composition.” *Nature* 506 (7486): 85.

Banerjee, Sudipto, Alan E Gelfand, Andrew O Finley, and Huiyan Sang. 2008. “Gaussian Predictive Process Models for Large Spatial Data Sets.” *Journal of the Royal Statistical Society: Series B (Statistical Methodology)* 70 (4): 825–48.

Cameletti, Michela, Finn Lindgren, Daniel Simpson, and Håvard Rue. 2013. “Spatio-Temporal Modeling of Particulate Matter Concentration Through the SPDE Approach.” *AStA Advances in Statistical Analysis* 97 (2): 109–31.

Diggle, Peter J, Raquel Menezes, and Ting-li Su. 2010. “Geostatistical Inference Under Preferential Sampling.” *Journal of the Royal Statistical Society: Series C (Applied Statistics)* 59 (2): 191–232.

Fong, Youyi, Håvard Rue, and Jon Wakefield. 2010. “Bayesian Inference for Generalized Linear Mixed Models.” *Biostatistics* 11 (3): 397–412.

Girolami, Mark, and Ben Calderhead. 2011. “Riemann Manifold Langevin and Hamiltonian Monte Carlo Methods.” *Journal of the Royal Statistical Society: Series B (Statistical Methodology)* 73 (2): 123–214.

Illian, Janine B, Sigrunn H Sørbye, and Håvard Rue. 2012. “A Toolbox for Fitting Complex Spatial Point Process Models Using Integrated Nested Laplace Approximation (Inla).” *The Annals of Applied Statistics*, 1499–1530.

Marin, Jean-Michel, Pierre Pudlo, Christian P Robert, and Robin J Ryder. 2012. “Approximate Bayesian Computational Methods.” *Statistics and Computing* 22 (6): 1167–80.

Martins, Thiago G, Daniel Simpson, Finn Lindgren, and Håvard Rue. 2013. “Bayesian Computing with Inla: New Features.” *Computational Statistics & Data Analysis* 67: 68–83.

Miller, David L, M Louise Burt, Eric A Rexstad, and Len Thomas. 2013. “Spatial Models for Distance Sampling Data: Recent Developments and Future Directions.” *Methods in Ecology and Evolution* 4 (11): 1001–10.

Murray, Iain, Ryan Prescott Adams, and David JC MacKay. 2010. “Elliptical Slice Sampling.”

Noor, Abdisalan M, Damaris K Kinyoki, Clara W Mundia, Caroline W Kabaria, Jonesmus W Mutua, Victor A Alegana, Ibrahima Socé Fall, and Robert W Snow. 2014. “The Changing Risk of Plasmodium Falciparum Malaria Infection in Africa: 2000–10: A Spatial and Temporal Analysis of Transmission Intensity.” *The Lancet* 383 (9930): 1739–47.

Schrödle, Birgit, and Leonhard Held. 2011. “Spatio-Temporal Disease Mapping Using Inla.” *Environmetrics* 22 (6): 725–34.

Scott, Steven L, Alexander W Blocker, Fernando V Bonassi, Hugh A Chipman, Edward I George, and Robert E McCulloch. 2016. “Bayes and Big Data: The Consensus Monte Carlo Algorithm.” *International Journal of Management Science and Engineering Management* 11 (2): 78–88.

Simpson, Daniel, Janine B Illian, Finn Lindgren, Sigrunn H Sørbye, and Havard Rue. 2016. “Going Off Grid: Computationally Efficient Inference for Log-Gaussian Cox Processes.” *Biometrika* 103 (1): 49–70.

Wang, Chong, and David M Blei. 2013. “Variational Inference in Nonconjugate Models.” *Journal of Machine Learning Research* 14 (Apr): 1005–31.


