use uhc_ssa, clear

tab dead

tab nUHC

tab country

tab dead nUHC, row

bysort country: tab dead nUHC, row

 mhodds dead nUHC , by(country)

/*STEPS in conducting Binary Data meta-analysis, i.e. mortality (yes or no), hypertensive (yes or no)

STEP 1 - LOAD DATA
STEP 2 - DECLARE, UPDATE & DESCRIBE meta data
STEP 3 - SUMMARIZE meta data by using a TABLE or a FOREST PLOT
STEP 4- EXPLORE HETEROGENEITY - SUB-GROUP and META-REGRESSION analysis
STEP 5- EXPLORE and ADDRESS SMALL-STUDY EFFECTS*/

Binary-outcome summaries        # of successes (treated)
                                # of failures (treated)
                                # of successes (controls)
                                # of failures (controls)

Binary type	Description
lnoratio	log odds-ratio; the default
lnrratio	log risk-ratio (also known as log rate ratio and log relative risk
rdiff	risk difference
lnorpeto	Peto’s log odds-ratio

*Experimental arm
gen dead1   = (nUHC == 1) & (dead == 1)
gen nodead1 = (nUHC == 1) & (dead == 0)

*Control arm
gen dead0   = (nUHC == 0) & (dead == 1)
gen nodead0 = (nUHC == 0) & (dead == 0)

preserve

collapse (sum) dead1 nodead1 dead0 nodead0, by(country)

list

meta esize dead1 nodead1 dead0 nodead0, studylabel(country) esize(lnrratio)

meta summarize, fixed

meta forestplot, random
graph display

#delimit ;
    meta forestplot, random nullrefline 
    columnopts(_data1, supertitle(No UHC))  
    columnopts(_data2, supertitle(UHC)) 
    columnopts(_a _c, title(Dead)) 
    columnopts(_b _d, title(Alive));
#delimit cr
graph display

use uhc_ssa, clear

global ind "kid_male mat_age noed mat_currwork mat_dec3 kid_bord kid_u5c mat_hinsur mat_fhead media_access mat_wealth"
global com "rural  com_poverty_hl com_uemp_hl com_illit_hl com_diversity_hl"
global country "country"

oaxaca dead  $ind $com , by(nUHC)   logit   pooled relax

coefplot (., keep(explained:*)),  bylabel("Explained") || ///
(., keep(unexplained:*)),  bylabel("Unexplained") || , ///
drop(*:_cons) recast(bar) barwidth(0.5) citop ciopts(recast(rcap) color(black)) ///
byopts(cols(1)) xline(0, lpattern(dash)) xlabel(, grid glstyle(minor_grid) glpattern(dash))


global factors "b4.country b1.kid_male b3.mat_cage b5.mat_wealth b2.mat_edu b1.rural b1.mat_currwork b3.media_access"

*PSM
#delimit ;
	psmatch2 CCI_non  $factors ,
	out(dead) logit ate  llr neighbor(5)
	;
#delimit cr

pstest, both  graph

graph display

gen pair1 = _id if _treated==0
replace pair1 = _n1 if _treated==1
gen pair2 = _id if _treated==0
replace pair2 = _n2 if _treated==1
gen pair3 = _id if _treated==0
replace pair3 = _n3 if _treated==1
gen pair4 = _id if _treated==0
replace pair4 = _n4 if _treated==1
gen pair5 = _id if _treated==0
replace pair5 = _n5 if _treated==1

bysort pair1: egen paircount1 = count(pair1)
bysort pair2: egen paircount2 = count(pair2)
bysort pair3: egen paircount3 = count(pair3)
bysort pair4: egen paircount4 = count(pair4)
bysort pair5: egen paircount5 = count(pair5)
egen byte paircount = anycount(paircount1 paircount2 paircount3 paircount4  paircount5), values(2)
drop if paircount==0

tab _treated

cs  dead _treated

global factors "i.country mat_hinsur i.kid_male i.mat_cage i.mat_wealth  i.rural i.mat_currwork i.media_access"

logit dead  nUHC $factors, or

regpar, at(nUHC=0)


Asymmetric 95% CIs for the untransformed proportions
under Scenario 0 and Scenario 1
and for the untransformed population attributable risk (PAR)
                Estimate     Minimum     Maximum 
  Scenario_0   .05863874   .04579201   .07480702 
  Scenario_1   .03765061   .02651294    .0532113 
         PAR   .02098813   .01134642   .03062594 

- In current (factual) scenario, the childhood mortality was 59 per 1000 
- but in the theoretical minimum scenario, the childhood mortality was 38 per 
- The proportion of childhood deaths that would have been avoided if all mother had access to UHC was 2.1%, 
suggesting that excess 21 childhood death could have been avoided for every 1000 children.


Asymmetric 95% CIs for the untransformed proportions
under Scenario 0 and Scenario 1
and for the untransformed population attributable risk (PAR)
                Estimate     Minimum     Maximum 
  Scenario_0   .05863874   .04579201   .07480702 
  Scenario_1   .03765061   .02651294    .0532113 
         PAR   .02098813   .01134642   .03062594 

We see that in the real world (Scenario 0), 5.9% of children are expected to die before the 5th birthday but 
that in the dream scenario where all mother have access to UHC (Scenario 1), 
only 3.8% of children are expected to die before the 5th birthday. 
The diﬀerence between these scenario percentages (PAR) is 2.1%, 
with conﬁdence limits from 3.2% to 13.5%. 

The PAR can be interpreted as the proportion of all babies that have low birthweight because 
they were born in scenario 0 instead of in scenario 1.

regpar, at(nUHC=0) subpop(if rural ==1)
regpar, at(nUHC=0) subpop(if rural ==0)

regpar, at(nUHC=0) subpop(if country == 1)
regpar, at(nUHC=0) subpop(if country == 2)
regpar, at(nUHC=0) subpop(if country == 3)
regpar, at(nUHC=0) subpop(if country == 4)

logit dead  nUHC `factor', or

punaf, at(nUHC=0) eform

*13.9% of the ‘disease burden’ of childhood mortality might be eliminated by providing access to UHC health service indicators to all child, with confidence limits from 8.9% to 18.5%. 

Scenario 0: (asobserved) _all
Scenario 1: nUHC=0
Confidence intervals for the means under Scenario 0 and Scenario 1
and for the population unattributable faction (PUF)
Total number of observations used: 962
------------------------------------------------------------------------------
             | Mean/Ratio   Std. Err.      z    P>|z|     [95% Conf. Interval]
-------------+----------------------------------------------------------------
  Scenario_0 |   .0582121   .0074087   -22.34   0.000     .0453606    .0747045
  Scenario_1 |   .0376884     .00675   -18.30   0.000     .0265312    .0535376
         PUF |   .6474336     .07635    -3.69   0.000     .5138252    .8157837
------------------------------------------------------------------------------

95% CI for the population attributable fraction (PAF)
                Estimate     Minimum     Maximum 
         PAF   .35256644    .1842163    .4861748 

- 35.2% of the ‘disease burden’ of childhood mortality might be eliminated by providing access to 
UHC health service indicators to all child, with confidence limits from 18.4% to 48.6%. 

punaf, at(nUHC=0) eform subpop(if rural ==1)
punaf, at(nUHC=0) eform subpop(if rural ==0)

punaf, at(nUHC=0) eform subpop(if country == 1)
punaf, at(nUHC=0) eform subpop(if country == 2)
punaf, at(nUHC=0) eform subpop(if country == 3)
punaf, at(nUHC=0) eform subpop(if country == 4)


