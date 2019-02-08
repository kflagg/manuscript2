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

Miller et al. (2013)
====================

Wang and Blei (n.d.)
====================

Cameletti et al. (2013)
=======================

Martins et al. (2013)
=====================

Noor et al. (2014)
==================

Fong, Rue, and Wakefield (2010)
===============================

Scott et al. (2016)
===================

References
==========

Andrieu, Christophe, Arnaud Doucet, and Roman Holenstein. 2010. “Particle Markov Chain Monte Carlo Methods.” *Journal of the Royal Statistical Society: Series B (Statistical Methodology)* 72 (3): 269–342.

Bagchi, Robert, Rachel E Gallery, Sofia Gripenberg, Sarah J Gurr, Lakshmi Narayan, Claire E Addis, Robert P Freckleton, and Owen T Lewis. 2014. “Pathogens and Insect Herbivores Drive Rainforest Plant Diversity and Composition.” *Nature* 506 (7486): 85.

Banerjee, Sudipto, Alan E Gelfand, Andrew O Finley, and Huiyan Sang. 2008. “Gaussian Predictive Process Models for Large Spatial Data Sets.” *Journal of the Royal Statistical Society: Series B (Statistical Methodology)* 70 (4): 825–48.

Cameletti, Michela, Finn Lindgren, Daniel Simpson, and Håvard Rue. 2013. “Spatio-Temporal Modeling of Particulate Matter Concentration Through the Spde Approach.” *AStA Advances in Statistical Analysis* 97 (2): 109–31.

Diggle, Peter J, Raquel Menezes, and Ting-li Su. 2010. “Geostatistical Inference Under Preferential Sampling.” *Journal of the Royal Statistical Society: Series C (Applied Statistics)* 59 (2): 191–232.

Fong, Youyi, Håvard Rue, and Jon Wakefield. 2010. “Bayesian Inference for Generalized Linear Mixed Models.” *Biostatistics* 11 (3): 397–412.

Girolami, Mark, and Ben Calderhead. 2011. “Riemann Manifold Langevin and Hamiltonian Monte Carlo Methods.” *Journal of the Royal Statistical Society: Series B (Statistical Methodology)* 73 (2): 123–214.

Marin, Jean-Michel, Pierre Pudlo, Christian P Robert, and Robin J Ryder. 2012. “Approximate Bayesian Computational Methods.” *Statistics and Computing* 22 (6): 1167–80.

Martins, Thiago G, Daniel Simpson, Finn Lindgren, and Håvard Rue. 2013. “Bayesian Computing with Inla: New Features.” *Computational Statistics & Data Analysis* 67: 68–83.

Miller, David L, M Louise Burt, Eric A Rexstad, and Len Thomas. 2013. “Spatial Models for Distance Sampling Data: Recent Developments and Future Directions.” *Methods in Ecology and Evolution* 4 (11): 1001–10.

Murray, Iain, Ryan Prescott Adams, and David JC MacKay. 2010. “Elliptical Slice Sampling.”

Noor, Abdisalan M, Damaris K Kinyoki, Clara W Mundia, Caroline W Kabaria, Jonesmus W Mutua, Victor A Alegana, Ibrahima Socé Fall, and Robert W Snow. 2014. “The Changing Risk of Plasmodium Falciparum Malaria Infection in Africa: 2000–10: A Spatial and Temporal Analysis of Transmission Intensity.” *The Lancet* 383 (9930): 1739–47.

Schrödle, Birgit, and Leonhard Held. 2011. “Spatio-Temporal Disease Mapping Using Inla.” *Environmetrics* 22 (6): 725–34.

Scott, Steven L, Alexander W Blocker, Fernando V Bonassi, Hugh A Chipman, Edward I George, and Robert E McCulloch. 2016. “Bayes and Big Data: The Consensus Monte Carlo Algorithm.” *International Journal of Management Science and Engineering Management* 11 (2): 78–88.

Wang, Chong, and David M Blei. n.d. “Variational Inference in Nonconjugate Models.” *Journal of Machine Learning Research* 14 (Apr): 1005–31.


