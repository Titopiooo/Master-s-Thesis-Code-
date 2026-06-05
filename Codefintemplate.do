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
For the map, the Eurostat NUTS shapefile is needed. If the shapefile is not available,
the map section is skipped without stopping the do-file.
****************************************************************************************/

clear all
set more off
* Optional: set your own working directory before running the file
* cd "your/path/here"

* Optional output folder
capture mkdir "outputs"

* Optional map settings
* Put the Eurostat NUTS shapefile in this folder if you want to reproduce the map.
local mapdir "NUTS_RG_20M_2021_3035-2"
local shapefile "NUTS_RG_20M_2021_3035.shp"

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
* 2.1) NEW: Inference for pooled concentration index of observed mortality
*      This produces Appendix Table A1, Panel A
*------------------------------------------------------------------------------
preserve

keep if !missing(mortstd, fracrank)

summ fracrank
scalar var_rank = r(Var)

summ mortstd, meanonly
scalar mean_y = r(mean)

gen lhs_ci = 2 * var_rank * (mortstd / mean_y)

reg lhs_ci fracrank, robust

scalar CI_obs_reg = _b[fracrank]
scalar SE_obs_reg = _se[fracrank]
scalar P_obs_reg  = 2 * ttail(e(df_r), abs(CI_obs_reg / SE_obs_reg))
scalar LB_obs_reg = CI_obs_reg - invttail(e(df_r), 0.025) * SE_obs_reg
scalar UB_obs_reg = CI_obs_reg + invttail(e(df_r), 0.025) * SE_obs_reg

display "Observed mortality CI with robust SE"
display "CI = " CI_obs_reg
display "SE = " SE_obs_reg
display "p-value = " P_obs_reg
display "Lower 95% CI = " LB_obs_reg
display "Upper 95% CI = " UB_obs_reg

restore

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

********************************************************************************
* NEW: APPENDIX TABLE A1, PANEL B
* Yearly concentration indices with robust standard errors and confidence intervals
********************************************************************************

tempfile yearlyci
tempname memhold

postfile `memhold' int year double ci se p lb ub using `yearlyci', replace

levelsof year, local(years)

foreach yy of local years {

    preserve

    keep if year == `yy' & !missing(mortstd, gdp)

    count
    local N = r(N)

    sort gdp
    egen rank_y = rank(gdp), unique
    gen fracrank_y = (rank_y - 0.5) / `N'

    summarize fracrank_y
    scalar var_rank_y = r(Var)

    summarize mortstd, meanonly
    scalar mean_y_year = r(mean)

    gen lhs_ci_y = 2 * var_rank_y * (mortstd / mean_y_year)

    quietly reg lhs_ci_y fracrank_y, robust

    scalar ci = _b[fracrank_y]
    scalar se = _se[fracrank_y]
    scalar p = 2 * ttail(e(df_r), abs(ci / se))
    scalar lb = ci - invttail(e(df_r), 0.025) * se
    scalar ub = ci + invttail(e(df_r), 0.025) * se

    post `memhold' (`yy') (ci) (se) (p) (lb) (ub)

    restore
}

postclose `memhold'

preserve

use `yearlyci', clear
format ci se p lb ub %9.5f

list, clean noobs

export excel using "outputs/yearly_CI_inference.xlsx", firstrow(variables) replace

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

*------------------------------------------------------------------------------
* NEW: Inference for concentration index of need-adjusted mortality
*      This completes Appendix Table A1, Panel A
*------------------------------------------------------------------------------
preserve

keep if !missing(mort_IS, fracrank)

summ fracrank
scalar var_rank_IS = r(Var)

summ mort_IS, meanonly
scalar mean_IS = r(mean)

gen lhs_ci_IS = 2 * var_rank_IS * (mort_IS / mean_IS)

reg lhs_ci_IS fracrank, robust

scalar CI_IS_reg = _b[fracrank]
scalar SE_IS_reg = _se[fracrank]
scalar P_IS_reg  = 2 * ttail(e(df_r), abs(CI_IS_reg / SE_IS_reg))
scalar LB_IS_reg = CI_IS_reg - invttail(e(df_r), 0.025) * SE_IS_reg
scalar UB_IS_reg = CI_IS_reg + invttail(e(df_r), 0.025) * SE_IS_reg

display "Need-adjusted mortality CI with robust SE"
display "CI = " CI_IS_reg
display "SE = " SE_IS_reg
display "p-value = " P_IS_reg
display "Lower 95% CI = " LB_IS_reg
display "Upper 95% CI = " UB_IS_reg

restore

* Export Appendix Table A1, Panel A
putexcel set "outputs/Appendix_Table_A1_Panel_A_main_CI.xlsx", replace
putexcel A1=("Measure") B1=("Concentration Index") C1=("Robust SE") D1=("p-value") E1=("95% CI lower") F1=("95% CI upper")
putexcel A2=("Observed mortality") B2=(CI_obs_reg) C2=(SE_obs_reg) D2=(P_obs_reg) E2=(LB_obs_reg) F2=(UB_obs_reg)
putexcel A3=("Need-adjusted mortality") B3=(CI_IS_reg) C3=(SE_IS_reg) D3=(P_IS_reg) E3=(LB_IS_reg) F3=(UB_IS_reg)

********************************************************************************
* NEW: TABLE 7
* Sensitivity check: alternative need definitions
********************************************************************************

preserve

capture drop yhat_need_age yhat_need_chronic yhat_need_both
capture drop mort_need_age mort_need_chronic mort_need_both

summ mortstd, meanonly
scalar mu_y_check = r(mean)

corr mortstd fracrank, covariance
scalar cov_y_check = r(cov_12)

scalar CI_y_check = 2 * cov_y_check / mu_y_check

summ mortstd, meanonly
scalar mean_mort = r(mean)

* Model A: need adjustment using age structure only
reg mortstd share65plus
predict yhat_need_age, xb
gen mort_need_age = mortstd - yhat_need_age + mean_mort

summ mort_need_age, meanonly
scalar mu_need_age = r(mean)

corr mort_need_age fracrank, covariance
scalar cov_need_age = r(cov_12)

scalar CI_need_age = 2 * cov_need_age / mu_need_age

* Model B: need adjustment using chronic conditions only
reg mortstd chroniccondition
predict yhat_need_chronic, xb
gen mort_need_chronic = mortstd - yhat_need_chronic + mean_mort

summ mort_need_chronic, meanonly
scalar mu_need_chronic = r(mean)

corr mort_need_chronic fracrank, covariance
scalar cov_need_chronic = r(cov_12)

scalar CI_need_chronic = 2 * cov_need_chronic / mu_need_chronic

* Model C: need adjustment using age structure and chronic conditions
reg mortstd share65plus chroniccondition
predict yhat_need_both, xb
gen mort_need_both = mortstd - yhat_need_both + mean_mort

summ mort_need_both, meanonly
scalar mu_need_both = r(mean)

corr mort_need_both fracrank, covariance
scalar cov_need_both = r(cov_12)

scalar CI_need_both = 2 * cov_need_both / mu_need_both

* Percentage reduction and remaining inequality
scalar red_age = (abs(CI_y_check) - abs(CI_need_age)) / abs(CI_y_check) * 100
scalar rem_age = abs(CI_need_age) / abs(CI_y_check) * 100

scalar red_chronic = (abs(CI_y_check) - abs(CI_need_chronic)) / abs(CI_y_check) * 100
scalar rem_chronic = abs(CI_need_chronic) / abs(CI_y_check) * 100

scalar red_both = (abs(CI_y_check) - abs(CI_need_both)) / abs(CI_y_check) * 100
scalar rem_both = abs(CI_need_both) / abs(CI_y_check) * 100

display "===================================================="
display "SENSITIVITY CHECK 1: ALTERNATIVE NEED DEFINITIONS"
display "Observed CI                         = " CI_y_check
display "Need-adjusted CI, age only           = " CI_need_age
display "Need-adjusted CI, chronic only       = " CI_need_chronic
display "Need-adjusted CI, age + chronic      = " CI_need_both
display "----------------------------------------------------"
display "Reduction age only (%)               = " red_age
display "Remaining age only (%)               = " rem_age
display "Reduction chronic only (%)           = " red_chronic
display "Remaining chronic only (%)           = " rem_chronic
display "Reduction age + chronic (%)          = " red_both
display "Remaining age + chronic (%)          = " rem_both
display "===================================================="

putexcel set "outputs/need_standardization_sensitivity.xlsx", replace
putexcel A1=("Specification") B1=("CI") C1=("Reduction_pct") D1=("Remaining_pct")
putexcel A2=("Observed mortality") B2=(CI_y_check) C2=(0) D2=(100)
putexcel A3=("Need-adjusted: age only") B3=(CI_need_age) C3=(red_age) D3=(rem_age)
putexcel A4=("Need-adjusted: chronic conditions only") B4=(CI_need_chronic) C4=(red_chronic) D4=(rem_chronic)
putexcel A5=("Need-adjusted: age + chronic conditions") B5=(CI_need_both) C5=(red_both) D5=(rem_both)

restore

********************************************************************************
* NEW: TABLE 8
* Decomposition-based horizontal inequity check
********************************************************************************

scalar need_contribution = contr_65 + contr_ch
scalar HI_decomp = CI_y - need_contribution

scalar need_share_pct = (need_contribution / CI_y) * 100
scalar HI_decomp_pct = (HI_decomp / CI_y) * 100

display "===================================================="
display "DECOMPOSITION-BASED HORIZONTAL INEQUITY CHECK"
display "Observed CI                         = " CI_y
display "Contribution of share65plus          = " contr_65
display "Contribution of chroniccondition     = " contr_ch
display "Total need contribution              = " need_contribution
display "HI = CI - need contribution          = " HI_decomp
display "Need contribution as % of CI         = " need_share_pct
display "HI as % of observed CI               = " HI_decomp_pct
display "===================================================="

putexcel set "outputs/decomposition_based_HI.xlsx", replace
putexcel A1=("Measure") B1=("Value")
putexcel A2=("Observed concentration index") B2=(CI_y)
putexcel A3=("Contribution of age structure") B3=(contr_65)
putexcel A4=("Contribution of chronic conditions") B4=(contr_ch)
putexcel A5=("Total need contribution") B5=(need_contribution)
putexcel A6=("Decomposition-based horizontal inequity index") B6=(HI_decomp)

********************************************************************************
* NEW: TABLE 9
* Leave-one-region-out sensitivity check
********************************************************************************

tempfile loo_results
tempname handle

postfile `handle' str40 excluded_region double CI_observed CI_need_adjusted using `loo_results', replace

levelsof region, local(regions)

foreach r of local regions {

    preserve

        keep if region != `"`r'"'

        * Recalculate GDP fractional rank after excluding one region
        capture drop rank_loo fracrank_loo
        sort gdp
        egen rank_loo = rank(gdp)
        gen fracrank_loo = (rank_loo - 0.5) / _N

        * Observed concentration index
        quietly summarize mortstd, meanonly
        scalar mu_obs_loo = r(mean)

        quietly corr mortstd fracrank_loo, covariance
        scalar cov_obs_loo = r(cov_12)

        scalar CI_obs_loo = 2 * cov_obs_loo / mu_obs_loo

        * Need-adjusted concentration index: age structure + chronic conditions
        capture drop yhat_loo mort_need_loo

        quietly regress mortstd share65plus chroniccondition
        predict yhat_loo, xb

        quietly summarize mortstd, meanonly
        scalar mean_mort_loo = r(mean)

        gen mort_need_loo = mortstd - yhat_loo + mean_mort_loo

        quietly summarize mort_need_loo, meanonly
        scalar mu_need_loo = r(mean)

        quietly corr mort_need_loo fracrank_loo, covariance
        scalar cov_need_loo = r(cov_12)

        scalar CI_need_loo = 2 * cov_need_loo / mu_need_loo

        post `handle' (`"`r'"') (CI_obs_loo) (CI_need_loo)

    restore
}

postclose `handle'

preserve

use `loo_results', clear

format CI_observed CI_need_adjusted %9.5f

list, clean noobs

summarize CI_observed CI_need_adjusted

export excel using "outputs/leave_one_region_out_CI.xlsx", firstrow(variables) replace

restore

********************************************************************************
* NEW: TABLE 3
* Descriptive statistics by macroarea
********************************************************************************

preserve

capture drop macroarea
gen str20 macroarea = ""

* North
replace macroarea = "North" if inlist(region, "PIEMONTE", "VALLE D'AOSTA", "LIGURIA", "LOMBARDIA", "BOLZANO")
replace macroarea = "North" if inlist(region, "TRENTO", "VENETO", "FRIULI-VENEZIA GIULIA", "EMILIA-ROMAGNA")

* Centre
replace macroarea = "Centre" if inlist(region, "TOSCANA", "UMBRIA", "MARCHE", "LAZIO")

* South and Islands
replace macroarea = "South and Islands" if inlist(region, "ABRUZZO", "MOLISE", "CAMPANIA", "PUGLIA", "BASILICATA")
replace macroarea = "South and Islands" if inlist(region, "CALABRIA", "SICILIA", "SARDEGNA")

tab macroarea
list region if macroarea == ""

assert macroarea != ""

collapse (mean) mortstd gdp healthcareexpend LEAscore piano_rientro, by(macroarea)

gen recovery_plan_share = piano_rientro * 100
drop piano_rientro

format mortstd %9.2f
format gdp %9.0f
format healthcareexpend %9.0f
format LEAscore %9.1f
format recovery_plan_share %9.1f

list macroarea mortstd gdp healthcareexpend LEAscore recovery_plan_share, clean noobs

export excel using "outputs/macroarea_descriptive_table.xlsx", firstrow(variables) replace

restore

********************************************************************************
* NEW: FIGURE 1
* Map of mean standardized mortality by region
********************************************************************************

preserve

collapse (mean) mortstd LEAscore healthcareexpend gdp, by(region)

save "outputs/mapvars_mean_region.dta", replace

restore

capture confirm file "`mapdir'/`shapefile'"

if _rc {
    display "Map shapefile not found. Map section skipped."
    display "To reproduce the map, place the Eurostat NUTS shapefile in the folder indicated by local mapdir."
}
else {

    preserve

    use "outputs/mapvars_mean_region.dta", clear

    gen str4 NUTS_ID = ""

    replace NUTS_ID = "ITF1" if region == "ABRUZZO"
    replace NUTS_ID = "ITF5" if region == "BASILICATA"
    replace NUTS_ID = "ITH1" if region == "BOLZANO"
    replace NUTS_ID = "ITF6" if region == "CALABRIA"
    replace NUTS_ID = "ITF3" if region == "CAMPANIA"
    replace NUTS_ID = "ITH5" if region == "EMILIA-ROMAGNA"
    replace NUTS_ID = "ITH4" if region == "FRIULI-VENEZIA GIULIA"
    replace NUTS_ID = "ITI4" if region == "LAZIO"
    replace NUTS_ID = "ITC3" if region == "LIGURIA"
    replace NUTS_ID = "ITC4" if region == "LOMBARDIA"
    replace NUTS_ID = "ITI3" if region == "MARCHE"
    replace NUTS_ID = "ITF2" if region == "MOLISE"
    replace NUTS_ID = "ITC1" if region == "PIEMONTE"
    replace NUTS_ID = "ITF4" if region == "PUGLIA"
    replace NUTS_ID = "ITG2" if region == "SARDEGNA"
    replace NUTS_ID = "ITG1" if region == "SICILIA"
    replace NUTS_ID = "ITI1" if region == "TOSCANA"
    replace NUTS_ID = "ITH2" if region == "TRENTO"
    replace NUTS_ID = "ITI2" if region == "UMBRIA"
    replace NUTS_ID = "ITC2" if region == "VALLE D'AOSTA"
    replace NUTS_ID = "ITH3" if region == "VENETO"

    list region NUTS_ID, clean noobs

    count if missing(NUTS_ID) | NUTS_ID == ""
    assert r(N) == 0

    save "outputs/mapvars_mean_region_nuts.dta", replace

    restore

    preserve

    cd "`mapdir'"

    spshape2dta "`shapefile'", replace

    local shpbase = subinstr("`shapefile'", ".shp", "", .)

    use "`shpbase'.dta", clear

    keep if CNTR_CODE == "IT" & LEVL_CODE == 2
    count
    list NUTS_ID NAME_LATN in 1/21, clean

    save "../outputs/italy_nuts2_attr.dta", replace

    restore

    preserve

    use "outputs/italy_nuts2_attr.dta", clear

    merge 1:1 NUTS_ID using "outputs/mapvars_mean_region_nuts.dta"

    tab _merge
    keep if _merge == 3
    drop _merge

    save "outputs/italy_map_merged.dta", replace

    capture which spmap
    if _rc {
        ssc install spmap, replace
    }

    spmap mortstd using "`mapdir'/`shpbase'_shp.dta", id(_ID) ///
        clmethod(custom) clbreaks(77 82 86 90 95 102) ///
        fcolor(eltblue ebblue edkblue dknavy navy) ///
        ocolor(white white white white white) ///
        osize(vthin vthin vthin vthin vthin) ///
        title("") note("") legend(pos(7) ring(0) size(small))

    graph save "outputs/Italycolored.gph", replace
    graph export "outputs/italy_mortality_map_blue.png", width(3000) replace

    restore
}


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

* Manual first-stage F and partial R-squared used for Appendix Table A2
quietly reg healthcareexpend L_exp share65plus chroniccondition pct_degree LEAscore piano_rientro i.year if !missing(L_exp), vce(cluster region_id)
test L_exp
scalar firststageF_manual = r(F)

quietly reg healthcareexpend share65plus chroniccondition pct_degree LEAscore piano_rientro i.year if !missing(L_exp)
scalar rss_reduced = e(rss)

quietly reg healthcareexpend L_exp share65plus chroniccondition pct_degree LEAscore piano_rientro i.year if !missing(L_exp)
scalar rss_full = e(rss)

scalar partialR2_manual = (rss_reduced - rss_full) / rss_reduced

display "First-stage F-statistic = " firststageF_manual
display "Partial R-squared = " partialR2_manual

estimates restore m3
estadd scalar firststageF = firststageF_manual, replace
estadd scalar partialR2 = partialR2_manual, replace

esttab m1 m2 m3 using "outputs/Appendix_Table_A2.rtf", replace label ///
    mtitles("(1) Robust OLS full sample" "(2) Robust OLS IV sample" "(3) IV") ///
    cells("b(star fmt(4)) p(par fmt(4))") ///
    stats(N r2 yearfe clusterse firststageF partialR2, ///
    labels("Observations" "R-squared" "Year fixed effects" "Region-clustered standard errors" "First-stage F-statistic" "Partial R-squared")) ///
    star(* 0.10 ** 0.05 *** 0.01)

restore

********************************************************************************
* NEW: TABLE 11
* Wild-cluster bootstrap inference for selected coefficients
********************************************************************************

capture which boottest
if _rc ssc install boottest, replace

capture confirm variable region_id
if _rc encode region, gen(region_id)

xtset region_id year

reg mortstd healthcareexpend share65plus chroniccondition pct_degree LEAscore piano_rientro i.year, vce(cluster region_id)

boottest healthcareexpend, cluster(region_id) reps(9999) seed(12345)

boottest LEAscore, cluster(region_id) reps(9999) seed(12345)

boottest piano_rientro, cluster(region_id) reps(9999) seed(12345)

********************************************************************************
* END OF ROBUSTNESS AND IV ANALYSIS
********************************************************************************

********************************************************************************
* End of do-file
********************************************************************************

