
# RCMsize

<!-- badges: start -->
<!-- badges: end -->

The goal of RCMsize is to calculate the sample size required for studies where the main outcome measure is the seroconversion rate (SCR).
It provides tools to compute the probability of an individual being seropositive, given specific parameters such as the seroconversion rate (SCR) , the seroreversion rate (SRR), and the individual's age, using a reversible catalytic model.
Additionally, the package allows for the calculation of seroprevalence (SP) and its corresponding confidence interval, as well as the confidence interval for the seroconversion rate (SCR).

## Installation

To install the RCMsize package from GitHub, you first need to have the devtools package installed. If you haven't done so yet, you can easily install it using the following command:
```{r, eval=FALSE}
install.packages("devtools")
```
Once you have devtools installed, you can proceed to install the RCMsize package with the following command:
```{r, eval=FALSE}
devtools::install_github("https://github.com/marciagraca/RCMsize")
```
After the installation is complete, you can load the RCMsize package into your R session by running:

```{r setup}
library(RCMsize)
```


## Example

In this example we calculate the required sample size so that the relative width of the confidence interval for the seroconversion rate (SCR) is equal to a specified value (`RL`). The function calculates the necessary sample size by iteratively adjusting the sample size until the confidence interval for SCR meets the desired width criteria. This calculation is based on the seroprevalence, confidence intervals for seroprevalence, and the seroconversion rate.

### Parameters:
- `SCR`: The seroconversion rate.
- `RL`: The desired relative width for the confidence interval width for the seroconversion rate.
- `SRR`: The seroreversion rate.
- `ages`: A vector representing the distribution of ages in the population.
- `A_max`: The maximum age in the population.
- `limits`: The lower and upper limits for the `SCR`.
- `max_iter`: The maximum number of iterations for adjusting the sample size (default is 10000).
- `conf.level`: The confidence level for the confidence interval (default is 0.95).
- `method`: The method for calculating the confidence interval for seroprevalence, default is `"asymptotic"`.


### Example:

```{r}
A_max <- 80
age_distribution <- rep(1 / A_max, A_max)
sample_s(0.03, 1, 0.01, age_distribution, A_max, limits = c(0, 1))
```

# References

For more information on reversible catalytic models, please refer to the following [article](https://rdcu.be/d0lp5).


