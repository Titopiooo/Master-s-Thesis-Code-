/****************************************************************************************
Project: Regional inequalities in health outcomes in Italy
Purpose: Replication file for concentration index, decomposition, indirect standardization,
         and robustness / IV analysis
Dataset: fintemplate.dta
Unit of analysis: Region-year
Period: 2010–2019
Software: Stata

Description:
This do-file computes the concentration index of mortality, decomposes inequality into
its observable determinants, performs yearly robustness checks, indirectly standardizes
mortality for need variables, and estimates robustness and IV models for healthcare
expenditure.

Note:
The dataset is not publicly included in the repository. Users should place
"fintemplate.dta" in the working directory before running the file.
****************************************************************************************/

clear all
set more off

* Optional: set your own working directory before running the file
* cd "your/path/here"

*------------------------------------------------------------------------------
* 0) Load data
*------------------------------------------------------------------------------
use "fintemplate.dta", clear

* Quick sanity checks
describe
summ year mortstd gdp chroniccondition share65plus pct_degree LEAscore piano_rientro healthcareexpend
count
* Expected: 210 obs

*------------------------------------------------------------------------------
* 1) Create GDP rank and fractional rank (needed for CI and covariance formula)
*------------------------------------------------------------------------------
* IMPORTANT:
* - rank() assigns rank 1..N (poorest to richest if gdp sorted ascending)
* - fractional rank: (rank - 0.5) / N

sort gdp
egen rank = rank(gdp)
gen fracrank = (rank - 0.5) / _N

label var rank     "rank of (gdp)"
label var fracrank "fractional rank by GDP (rank-0.5)/N"

*------------------------------------------------------------------------------
* 2) Concentration Index (CI) of the outcome: mortstd
*    CI = 2 * cov(y, fracrank) / mean(y)
*------------------------------------------------------------------------------
summ mortstd
scalar mu_y = r(mean)

corr mortstd fracrank, covariance
scalar cov_y = r(cov_12)

scalar CI_y = 2 * cov_y / mu_y
display "Mean(y) mu_y = " mu_y
display "cov(y, fracrank) = " cov_y
display "CI_y (mortstd) = " CI_y

*------------------------------------------------------------------------------
* 3) Regression model (pooled 2010–2019)
*    NOTE: Used for decomposition coefficients; not interpreted as causal.
*------------------------------------------------------------------------------
reg mortstd share65plus chroniccondition pct_degree LEAscore piano_rientro healthcareexpend

* Optional: multicollinearity check
vif

*------------------------------------------------------------------------------
* 4) Store regression coefficients for the decomposition
*------------------------------------------------------------------------------
scalar b_share65 = _b[share65plus]
scalar b_chronic = _b[chroniccondition]
scalar b_degree  = _b[pct_degree]
scalar b_lea     = _b[LEAscore]
scalar b_pr      = _b[piano_rientro]
scalar b_exp     = _b[healthcareexpend]

*------------------------------------------------------------------------------
* 5) Decomposition pieces for each x:
*    - mu_x = mean(x)
*    - CI_x = 2 * cov(x, fracrank) / mu_x
*    - elasticity_x = (beta_x * mu_x) / mu_y
*    - contribution_x = elasticity_x * CI_x
*    - percent_x = (contribution_x / CI_total) * 100
*------------------------------------------------------------------------------
* Helper procedure repeated for each variable (explicit, to keep it transparent)

*************************************
* 5.1 healthcareexpend
*************************************
summ healthcareexpend
scalar mu_exp = r(mean)

corr healthcareexpend fracrank, covariance
scalar cov_exp = r(cov_12)

scalar CI_exp = 2 * cov_exp / mu_exp
scalar elas_exp = (b_exp * mu_exp) / mu_y
scalar contr_exp = elas_exp * CI_exp
scalar perc_exp = (contr_exp / CI_y) * 100

display "---- healthcareexpend ----"
display "mu_exp=" mu_exp " CI_exp=" CI_exp " elasticity=" elas_exp " contr=" contr_exp " %=" perc_exp

*************************************
* 5.2 pct_degree
*************************************
summ pct_degree
scalar mu_deg = r(mean)

corr pct_degree fracrank, covariance
scalar cov_deg = r(cov_12)

scalar CI_deg = 2 * cov_deg / mu_deg
scalar elas_deg = (b_degree * mu_deg) / mu_y
scalar contr_deg = elas_deg * CI_deg
scalar perc_deg = (contr_deg / CI_y) * 100

display "---- pct_degree ----"
display "mu_deg=" mu_deg " CI_deg=" CI_deg " elasticity=" elas_deg " contr=" contr_deg " %=" perc_deg

*************************************
* 5.3 LEAscore
*************************************
summ LEAscore
scalar mu_lea = r(mean)

corr LEAscore fracrank, covariance
scalar cov_lea = r(cov_12)

scalar CI_lea = 2 * cov_lea / mu_lea
scalar elas_lea = (b_lea * mu_lea) / mu_y
scalar contr_lea = elas_lea * CI_lea
scalar perc_lea = (contr_lea / CI_y) * 100

display "---- LEAscore ----"
display "mu_lea=" mu_lea " CI_lea=" CI_lea " elasticity=" elas_lea " contr=" contr_lea " %=" perc_lea

*************************************
* 5.4 piano_rientro (binary)
*************************************
summ piano_rientro
scalar mu_pr = r(mean)

corr piano_rientro fracrank, covariance
scalar cov_pr = r(cov_12)

scalar CI_pr = 2 * cov_pr / mu_pr
scalar elas_pr = (b_pr * mu_pr) / mu_y
scalar contr_pr = elas_pr * CI_pr
scalar perc_pr = (contr_pr / CI_y) * 100

display "---- piano_rientro ----"
display "mu_pr=" mu_pr " CI_pr=" CI_pr " elasticity=" elas_pr " contr=" contr_pr " %=" perc_pr

*************************************
* 5.5 share65plus
*************************************
summ share65plus
scalar mu_65 = r(mean)

corr share65plus fracrank, covariance
scalar cov_65 = r(cov_12)

scalar CI_65 = 2 * cov_65 / mu_65
scalar elas_65 = (b_share65 * mu_65) / mu_y
scalar contr_65 = elas_65 * CI_65
scalar perc_65 = (contr_65 / CI_y) * 100

display "---- share65plus ----"
display "mu_65=" mu_65 " CI_65=" CI_65 " elasticity=" elas_65 " contr=" contr_65 " %=" perc_65

*************************************
* 5.6 chroniccondition
*************************************
summ chroniccondition
scalar mu_ch = r(mean)

corr chroniccondition fracrank, covariance
scalar cov_ch = r(cov_12)

scalar CI_ch = 2 * cov_ch / mu_ch
scalar elas_ch = (b_chronic * mu_ch) / mu_y
scalar contr_ch = elas_ch * CI_ch
scalar perc_ch = (contr_ch / CI_y) * 100

display "---- chroniccondition ----"
display "mu_ch=" mu_ch " CI_ch=" CI_ch " elasticity=" elas_ch " contr=" contr_ch " %=" perc_ch

*------------------------------------------------------------------------------
* 6) Check: sum of explained contributions vs CI_y
*    CI_y ≈ sum(contributions) + residual
*------------------------------------------------------------------------------
scalar sum_contrib = contr_exp + contr_deg + contr_lea + contr_pr + contr_65 + contr_ch
scalar residual = CI_y - sum_contrib

display "=============================="
display "CI_y           = " CI_y
display "sum_contrib    = " sum_contrib
display "residual (CI_y - sum) = " residual
display "=============================="

********************************************************************************
* DIAGNOSTIC: Is inequality driven by CI or elasticity?
********************************************************************************

scalar tot_CI = abs(CI_exp) + abs(CI_deg) + abs(CI_lea) + abs(CI_pr) + abs(CI_65) + abs(CI_ch)
scalar tot_elas = abs(elas_exp) + abs(elas_deg) + abs(elas_lea) + abs(elas_pr) + abs(elas_65) + abs(elas_ch)

preserve
clear
set obs 6

gen str20 factor = ""
gen CI_x = .
gen elasticity = .
gen contribution = .
gen share_CI = .
gen share_elas = .
gen str20 driver = ""

* HEALTH EXP
replace factor = "Health exp." in 1
replace CI_x = CI_exp in 1
replace elasticity = elas_exp in 1
replace contribution = contr_exp in 1

* EDUCATION
replace factor = "Education" in 2
replace CI_x = CI_deg in 2
replace elasticity = elas_deg in 2
replace contribution = contr_deg in 2

* LEA
replace factor = "LEA score" in 3
replace CI_x = CI_lea in 3
replace elasticity = elas_lea in 3
replace contribution = contr_lea in 3

* PIANO DI RIENTRO
replace factor = "Recovery plan" in 4
replace CI_x = CI_pr in 4
replace elasticity = elas_pr in 4
replace contribution = contr_pr in 4

* SHARE 65+
replace factor = "Share 65+" in 5
replace CI_x = CI_65 in 5
replace elasticity = elas_65 in 5
replace contribution = contr_65 in 5

* CHRONIC
replace factor = "Chronic cond." in 6
replace CI_x = CI_ch in 6
replace elasticity = elas_ch in 6
replace contribution = contr_ch in 6

* shares
replace share_CI = abs(CI_x) / tot_CI
replace share_elas = abs(elasticity) / tot_elas

* driver
replace driver = "More CI" if share_CI > share_elas
replace driver = "More elasticity" if share_elas >= share_CI

format CI_x elasticity contribution %9.4f
format share_CI share_elas %9.3f

list factor CI_x elasticity contribution share_CI share_elas driver, noobs sep(0)

restore

********************************************************************************
* Robustness check: Concentration Index of mortstd by year (2010, 2015, 2019)
********************************************************************************

* Store the pooled CI (optional, just for comparison)
summ mortstd, meanonly
scalar mu_y_pooled = r(mean)
corr mortstd fracrank, covariance
scalar cov_y_pooled = r(cov_12)
scalar CI_y_pooled = 2 * cov_y_pooled / mu_y_pooled
display "Pooled CI_y (mortstd) = " CI_y_pooled

* Loop over selected years
foreach t in 2010 2015 2019 {

    preserve
        keep if year == `t'
        count
        display "-----------------------------"
        display "Year = `t' | N = " r(N)

        sort gdp
        egen rank_y = rank(gdp)
        gen fracrank_y = (rank_y - 0.5) / _N

        summ mortstd, meanonly
        scalar mu_y_year = r(mean)

        corr mortstd fracrank_y, covariance
        scalar cov_y_year = r(cov_12)

        scalar CI_y_year = 2 * cov_y_year / mu_y_year
        display "CI_y (mortstd) in `t' = " CI_y_year

    restore
}

********************************************************************************
* ROBUSTNESS CHECK: Concentration Index of mortality by year (2010–2019)
********************************************************************************

preserve

* Create variable to store yearly concentration index
gen CI_year = .

* Get list of years
levelsof year, local(years)

* Loop over all years
foreach y of local years {

    * Sort observations by GDP (ranking variable)
    sort year gdp

    * Compute rank and fractional rank within each year
    by year: egen rank_temp = rank(gdp) if year == `y'
    by year: gen fracrank_temp = (rank_temp - 0.5) / _N if year == `y'

    * Compute mean mortality
    quietly summarize mortstd if year == `y', meanonly
    scalar mu = r(mean)

    * Compute covariance between mortality and fractional rank
    quietly corr mortstd fracrank_temp if year == `y', covariance
    scalar cov = r(cov_12)

    * Compute concentration index
    scalar CI_temp = 2 * cov / mu

    * Store result
    replace CI_year = CI_temp if year == `y'

    * Drop temporary variables
    drop rank_temp fracrank_temp
}

* Collapse dataset to one observation per year
collapse (mean) CI_year, by(year)

* Sort and display results
sort year
list, clean

********************************************************************************
* FIGURE: Yearly concentration index of mortality (Italy, 2010–2019)
********************************************************************************

twoway line CI_year year, ///
    xtitle("Year") ///
    ytitle("Concentration index of mortality") ///
    title("Yearly concentration index of mortality (Italy, 2010–2019)") ///
    xlabel(2010(1)2019) ///
    yline(0)

graph save "Concentration_index_graph.gph", replace

restore

************************************************
* Optional bar chart of percentage contributions from the decomposition
************************************************

preserve
clear
set obs 7

* Create variables for graph
gen str30 factor = ""
gen contribution = .

* Insert percentage contributions
replace factor = "Health exp." in 1
replace contribution = perc_exp in 1

replace factor = "Education" in 2
replace contribution = perc_deg in 2

replace factor = "LEA score" in 3
replace contribution = perc_lea in 3

replace factor = "Recovery plan" in 4
replace contribution = perc_pr in 4

replace factor = "Share 65+" in 5
replace contribution = perc_65 in 5

replace factor = "Chronic cond." in 6
replace contribution = perc_ch in 6

replace factor = "Residual" in 7
replace contribution = (residual/CI_y)*100 in 7

* Order bars by absolute contribution (largest first)
gen abs_contr = abs(contribution)
gsort -abs_contr

* Plot bar chart
graph bar contribution, over(factor, sort(1) label(angle(30))) ///
    ytitle("Contribution to inequality (%)") ///
    title("Decomposition of regional inequality in mortality (2010–2019)") ///
    legend(off)

* Save graph
graph save "Graphmigliorato.gph", replace

restore

********************************************************************************
* INDIRECT STANDARDIZATION OF MORTALITY FOR NEED VARIABLES
* Need variables: share65plus, chroniccondition
********************************************************************************

reg mortstd share65plus chroniccondition

predict yhat_need

summ mortstd, meanonly
gen mort_IS = mortstd - yhat_need + r(mean)

corr mort_IS fracrank, covariance
scalar cov_IS = r(cov_12)

summ mort_IS, meanonly
scalar mu_IS = r(mean)

scalar CI_IS = 2 * cov_IS / mu_IS

display "Original CI = " CI_y
display "Indirectly standardized CI (inequity) = " CI_IS

********************************************************************************
* ROBUSTNESS AND INSTRUMENTAL VARIABLE ANALYSIS
* This section reproduces the final regressions reported in Table 6 / Appendix A1
********************************************************************************

preserve

*------------------------------------------------------------------------------
* 1) Prepare panel structure and lagged expenditure variable
*------------------------------------------------------------------------------
capture confirm variable region_id
if _rc encode region, gen(region_id)

xtset region_id year

capture drop L_exp
gen L_exp = L.healthcareexpend

*------------------------------------------------------------------------------
* 2) Robust OLS - full sample (N = 210)
*    Includes year fixed effects and region-clustered standard errors
*------------------------------------------------------------------------------
reg mortstd healthcareexpend share65plus chroniccondition pct_degree LEAscore piano_rientro i.year, vce(cluster region_id)

*------------------------------------------------------------------------------
* 3) Robust OLS - IV sample only (N = 189)
*    Restrict sample to observations with non-missing lagged expenditure
*    This makes OLS directly comparable with the IV specification
*------------------------------------------------------------------------------
reg mortstd healthcareexpend share65plus chroniccondition pct_degree LEAscore piano_rientro i.year if !missing(L_exp), vce(cluster region_id)

*------------------------------------------------------------------------------
* 4) IV regression - healthcare expenditure instrumented with one-year lag
*    Same sample as column (3), with year fixed effects and clustered standard errors
*------------------------------------------------------------------------------
ivregress 2sls mortstd (healthcareexpend = L_exp) share65plus chroniccondition pct_degree LEAscore piano_rientro i.year if !missing(L_exp), vce(cluster region_id)

*------------------------------------------------------------------------------
* 5) First-stage diagnostics
*------------------------------------------------------------------------------
estat firststage

restore

********************************************************************************
* END OF ROBUSTNESS AND IV ANALYSIS
********************************************************************************

********************************************************************************
* End of do-file
********************************************************************************
