In this repository, one can find the R routines used in the article [Estimating the Efficiency of Brazilian Electricity Distribution Utilities](https://doi.org/10.1080/02664763.2021.1890000). It is a joint work with Ralph S. Silva, Mario Jorge Mendonça and Amaro Olimpio Pereira Jr., published in the Journal of Applied Statistics.

The article proposes an alternative methodology from the Brazilian Electricity Regulatory Agency on the efficiency estimation for the Brazilian electricity distribution sector. Our proposal combines robust state-space models and stochastic frontier analysis to measure the operational cost efficiency in a panel data set from 60 Brazilian electricity distribution utilities. The modeling joins the main literature in energy economics with advanced econometric and statistic techniques in order to estimate the efficiencies. Moreover, the suggested model is able to deal with changes in the inefficiencies across time whilst the Bayesian paradigm – through Markov chain Monte Carlo techniques – facilitates the inference on all unknowns. The method enables a significant degree of flexibility in the resultant efficiencies and a complete photography of the distribution sector.

This repo includes:

- **MCMC_Code.R** : Main file containing the MCMC routine; 
- **Full_Conditionals.R** : Auxiliary file with all full conditional distributions;
- **base.RData** : real data set.
