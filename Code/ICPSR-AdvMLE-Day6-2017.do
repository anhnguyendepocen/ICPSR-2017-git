* ICPSR Advanced MLE -- 2017
*
* Day Six Stata code
*
* Note: Results are very slightly different than when using
* R due to differences in algorithms.

* Stratification:

ssc install rnd, replace

clear

set obs 400
gen Z = invnorm(uniform())
gen X = 0
replace X = 1 in 201/400
gen ZG = exp(2 + 0.5*Z)
rndwei 400 1 ZG
rename xw T0
rndwei 400 0.5 ZG
rename xw T1
gen T = T0
replace T = T1 if X==1
gen C=1
drop T0 T1

* Plot

stset T, failure(C)

sts graph, by(X) legend(order(1 "p = 0.5" 2 "p = 1"))

* Models:

stcox X Z, nohr
predict NoStrat, basesurv

stcox Z, strata(X) nohr
predict Strat, basesurv

* Plot:

twoway (line NoStrat T, sort lcolor(black) lwidth(medthick)) /*
 */ (line Strat T if X==0, sort lcolor(blue) lwidth(medthick) lpattern(dash)) /*
 */ (line Strat T if X==1, sort lcolor(cranberry) lwidth(medthick) /*
 */ lpattern(dash)), ytitle(Predicted Survival) xtitle(Time) legend(off)

* Weibull:

streg Z, strata(X) dist(weibull)

***************************************
* Duration dependence:
*
* Weibull model on SCOTUS data:

insheet using "scotus.csv", clear
stset service, failure(retire) id(justice)

streg age pension pagree, nohr dist(weibull)

* Add -age- to the shape parameter:

streg age pension pagree, nohr dist(weibull) anc(age)
mat B = e(b)
svmat B, names(Bs)
replace Bs5 = Bs5[_n-1] if Bs5[_n-1]!=.
replace Bs6 = Bs6[_n-1] if Bs6[_n-1]!=.

* Subsequent plot:

range Agesim 32 91
gen Phat = exp(Bs6 + Bs5*Agesim)

twoway (line Phat Agesim, sort lcolor(black) lwidth(medthick)), /*
 */ ytitle(Estimate of p) yline(1, lwidth(medium) lpattern(dash) /*
 */ lcolor(gs8)) xtitle(Age) legend(off)



***************************************
* Competing risks

clear

insheet using "scotus2.csv"
summarize

* K-M plots:

stset service, en(t svcstart) failure(retire)
sts gen RetireKM = s RetUB = ub(s) RetLB = lb(s)

stset service, en(t svcstart) failure(death)
sts gen DeathKM = s DUB = ub(s) DLB = lb(s)

twoway (line RetireKM service, sort lcolor(black) lwidth(medthick) /*
 */ connect(stairstep)) (rline RetLB RetUB service, sort lpattern(dash) /*
 */ connect(stairstep)) (line DeathKM service, sort lcolor(cranberry) /*
 */ lwidth(medthick) lpattern(solid) connect(stairstep)) /*
 */ (rline DLB DUB service, sort lcolor(cranberry) lpattern(dash) /*
 */ connect(stairstep)), ytitle(Predicted Survival) xtitle(Time) legend(off)
 
* Models:

gen anyevent = retire + death
stset service, en(t svcstart) failure(anyevent)
stcox age chief south pension pagree, nohr efron

stset service, en(t svcstart) failure(retire)
stcox age chief south pension pagree, nohr efron

stset service, en(t svcstart) failure(death)
stcox age chief south pension pagree, nohr efron

* MNL:

gen lnT = ln(service+1)
mlogit threecat age chief south pension pagree lnT

**********************************
* Repeated events:

clear

use "OR.dta"

gen eventno=.
sort dyadid year
quietly by dyadid : replace eventno=sum(dispute)+1
replace eventno=eventno-1 if dispute==1
gen altduration=1
sort dyadid year
quietly by dyadid: replace altduration=altduration[_n-1]+1 /*
  */ if altduration[_n-1]~=.
gen altstart=altduration-1

* First events:

stset altduration, en(altstart) failure(dispute)
stcox allies contig capratio growth democ trade if eventno==1, /*
  */ nohr efron robust cluster(dyadid)

* Andersen-Gill:

stcox allies contig capratio growth democ trade, /*
  */ nohr efron robust cluster(dyadid)

* PWP-Elapsed Time:

stcox allies contig capratio growth democ trade, strata(eventno) /*
  */ nohr efron robust cluster(dyadid)

* PWP-Gap Time:

stset stop, en(start) failure(dispute)
stcox allies contig capratio growth democ trade, strata(eventno) /*
  */ nohr efron robust cluster(dyadid)

* Make WLW data (saving old):

save "OR.dta", replace
drop _*
expand 8
gen double dy=dyadid+(year/10000)
sort dyadid year
gen eventlabel = 1
sort dy
qui by dy : replace eventlabel = eventlabel[_n-1]+1 if eventlabel[_n-1]!=.
drop dy
sort dyadid year eventlabel
gen WLWdispute=0
replace WLWdispute=1 if dispute==1 & eventlabel==eventno
save "OR-WLW.dta"

* WLW model:

stset altduration, en(altstart) failure(WLWdispute)
stcox allies contig capratio growth democ trade, strata(eventno) /*
  */ nohr efron robust cluster(dyadid)

* Strata-by-covariate interaction:

clear
use "OR.dta"
stset stop, en(start) failure(dispute)
gen capXevent=capratio*eventno

stcox allies contig growth democ trade capratio capXevent, strata(eventno) /*
  */ nohr efron robust cluster(dyadid)
