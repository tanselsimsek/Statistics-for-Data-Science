# Statistics-for-Data-Science
---

This is a compulsory course given by Sapienza University of Rome to review the basics of statistics & probability (+ going a bit further!) with the aim at providing the fundamentals tools for:

- setting up a suitable probabilistic model;

- understanding the principles behind the main inferential problems: estimation, testing, model checking and forecasting, uncertainty quantification;

- implementing inference on observed data through the likelihood function using both optimization and simulation-based approximations (e.g. resampling, bootstrap, Monte Carlo ecc.)

- developing statistical computations within a suitable software environment (mainly R)

## Content
- **Masci-Simsek_part_A.Rmd** : A suitable simulation study in R to double-check the performance of the Randomized Max-Cut
Algorithm.

- **Masci-Simsek_part_B.Rmd** : A program in R to simulate the preferential attachment process, starting with 4 pages linked together as a directed cycle on 4 vertices, adding pages each with one outlink until there are 1 million pages, and using γ = 0.5 as the probability
a link is to a page chosen uniformly at random and 1 − γ = 0.5 as the probability a link is copied from existing links. Additionally, analyze the complimentary cumulative degree distribution.

- **final_version_hw2.R** : Starting with the idea of learning The Bart from data using a MoG, design a simulation study in order to compare the performance of different model selection techniques. 

- **masci_simsek.Rmd** : Graphically represent all the estimated graphs and try to draw some conclusion: are there clear co-activation differences between the two groups.
