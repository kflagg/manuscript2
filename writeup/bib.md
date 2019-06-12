---
bibliography: ../references.bib
---


# @illianetal \*

- Presents a methodological "toolbox" for using INLA to fit spatial marked Markov processes driven by latent Gaussian process intensity functions
- Discretize window into cells and model cell counts via Poisson GLMM with log link
- Linear predictor includes _p_ spatial functions and a nonlinear function of a constructed covariate
    - Spatial functions modeled as random walks constrained to sum to zero
    - Constructed covariate is a summary describing small-scale interaction among points
    - Constructed covariate allows interaction-like behavior in LGCP, easier to fit than Gibbs process because interaction acts directly on intensity
    - Relationship between constructed covariate and intensity modeled as random walk constrained to sum to zero
    - Constructed covariate used in the paper is distance from grid cell center to nearest point outside grid cell
- Also includes a residual term which can be used to assess model adequacy
- Choice of prior for spatial effect is tied to choice of grid size
- Suggest a grid "not finer than that given by the data"
- Model fit to a simulated repulsive Strauss process and posterior predictive check of _L_-function appeared reasonable
- Model fit to a simulated Thomas cluster process did not have strong enough clustering
- Example with many small clusters and a spatial trend in the intensity of cluster centers described the trend well but underestimated the strength of clustering
- Fit to homogeneous Poisson process does not include spurious spatial effects
- Discussion of other possible constructed covariates and difficulties of identifying different effects at same scale
- Can include covariates observed at different locations that the points being modeled
- Covariate term in linear predictor for intensity is prediction from spatial model for covariate
- Also illustrate model for marks with spatial GLMM for marks


# @simpsonetal \*


# @banerjeegelfandfinley \*

- Simplify high-dimensional spatial and spatiotemporal predictions by projecting the spatial process into a lower-dimensional subset observed at "knots"
- Smaller covariance matrix to invert
- MCMC requires "reasonably informative priors"
- Useful with tens of thousands of prediction points, needs hundreds of knots
- "[F]aster alternatives that avoid MCMC sampling can be employed" (cites the INLA paper)


# @andrieudoucetholenstein \*

- PMCMC combines SMC and MCMC to handle high-dimensional models
- SMC used to generate proposals
- Focus on state space models
- Particle degeneracy less of a problem than with pure SMC
- Efficient exploration wthout a lot of MCMC tuning but slower than INLA


# @girolamicalderhead

- Exploits Riemann manifold structure of parameter space for efficient MC sampling without tuning in high-dimensional situations
- LGCP example, notes that INLA could speed up computation
- Langevin diffusion/Metropolis adjusted Langevin algorithm(MALA)?
- Extend HMC and MALA, stochastic differential equations are involved


# @marinetal \*

- Review paper


# @murrayadamsmackayl \*

- MCMC for models with multivariate Gaussian priors
- Slower but more general than approximations
- No tuning parameters
- Better than Gibbs sampling when latent variables are correlated
- Sort of a modified MH with proposals on an ellipse that passes through the current state


# @bagchietal

- Experiment applying insecticide and fungicide in a rainforest to see the effects on plant development
- Weekly seed and seedling counts over 17 weeks
- Hierarchical model with negative binomial response and latent normal variables fit in R-INLA
- Priors not stated

# @digglemenezessu \*

- Geostatistical process model expanded by conditioning on a point process model for the observation locations
- Needed when the process of selecting observation locations depends on the response
- Describes ML inference by conditional simulation
- Discussion response by Rue, Martino, Mondal, and Chopin mentions INLA as applicable


# @schroedleheld

- Use INLA to fit hierarchical model with Poisson observations and a latent spatiotemporal Gaussian process to annual region-level counts of cattle salmonellosis in Switzerland
- Describe how to fit an Intrinsic Gaussian Markov random field model (with linear constraints to ensure parameter identifiability) in INLA


# @milleretal

- Discusses mapping dolphin population density using distance sampling along transects and a GAM for Poisson process intensity
- INLA mentioned as a way to quickly fit Bayesian models


# @wangblei \*


- Discusses Laplace variational inference and delta method variational inference as generic approaches for non-conjugate exponential family models
- Describes Laplace approximation as Taylor expansion around MAP point, resulting in Gaussian approximation
- Conclude that Laplace approximation is more accurate than delta method


# @camelettilindgrensimpsonrue

- Discusses Gaussian Markov random field as lower-dimensional alternative to Gaussian process with dense covariance for spatiotemporal prediction
- Can be computed quickly using INLA and stochastic PDEs
- SPDE approach uses finite element approximation (linear combination of basis functions) of spatial field on triangular grid
- Illustration of particulate matter modeling fit in R-INLA


# @martinssimpsonlingrenrue

- Presents improvements implemented in R-INLA in the four years after the original INLA paper
- Reviews a variety of applications where INLA had been previously shown to do well
- Summarizes INLA
- Adds some flexibility to allow more complicated dependence structures


# @nooretal

- Spatialtemporal model used to map prevalence of a malaria-spreading parasite across Africa
- Fit using R-INLA and the SPDE approach [@camelettilindgrensimpsonrue]
- Possibly questionable data collection and variable selection


# @fongruewakefield

- Compares PQL and INLA for several biostatistical examples
- Points out that INLA approximation can be inaccurate for binomial GLMMs with small denominators


# @scottetal \*

- Setting where datasets are so large that communication within a computing cluster is cost-prohibitive
- "Consensus Monte Carlo operates by running a separate Monte Carlo algorithm on each machine, and then averaging individual Monte Carlo draws across machines."
- Authors suggest this is useful when the parameter space is too high-dimensional for INLA to be efficient or when joint posteriors are needed
- Data partitioned into "shards" and MC samplers run independently on each shard
- "Consensus posterior" constructed from draws which are weighted averages of draws from each sampler


# References
