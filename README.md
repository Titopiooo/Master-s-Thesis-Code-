# Regional inequalities in health outcomes in Italy

This repository contains the Stata code used for the empirical analysis of regional inequalities in health outcomes in Italy over the period 2010–2019.

## Project overview

The analysis focuses on regional differences in standardized mortality and examines whether mortality is disproportionately concentrated among poorer regions. The empirical strategy includes:

- computation of the concentration index of mortality
- decomposition of the concentration index
- yearly robustness checks of the concentration index
- indirect standardization for need variables
- robustness regressions with year fixed effects and region-clustered standard errors
- instrumental variable analysis using lagged healthcare expenditure

## Unit of analysis

Region-year observations for Italian regions over the period 2010–2019.

## Main variables

- `mortstd`: standardized mortality rate
- `gdp`: GDP per capita
- `share65plus`: share of population aged 65 and over
- `pct_degree`: share of population with tertiary education
- `chroniccondition`: proxy for health need / chronic conditions prevalence
- `LEAscore`: regional LEA monitoring score
- `piano_rientro`: dummy for recovery plan status
- `healthcareexpend`: public healthcare expenditure per capita

## Files

- `Codefintemplate.do`: Stata do-file used for the main analysis

## Data availability

The dataset used for the analysis (`fintemplate.dta`) is not publicly shared in this repository. The code is provided for transparency and to document the analytical workflow used in the thesis.

## Software

The analysis was conducted in Stata.
