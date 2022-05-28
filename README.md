# Determination-of-the-optimal-number-of-strata-for-propensity-score-subclassification
This repository stores a R code to determine the optimal number of strata for propensity score subclassification proposed by Orihara and Hamada (2021).
Also, a toy example is prepared in the R code.

The R code needs to prepare the following variables:
- TT: treatment variable, PS: propensity score, YY: outcome variable

Unfortunately, the proposed method has some limitations:
- Treatment variables are only dichotomous.
- The maximum number of strata is necessary to use the method.

Also please note the following when use the R code:
- The programs has not validated correctly such as double programming (there is only self check).
- The programs are freely available; however, there are no responsibility.
