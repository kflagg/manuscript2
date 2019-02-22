---
bibliography: ../references.bib
---

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

# @wangblei

# @camelettilindgrensimpsonrue

# @martinssimpsonlingrenrue

# @nooretal

# @fongruewakefield

# @scottetal

# References
