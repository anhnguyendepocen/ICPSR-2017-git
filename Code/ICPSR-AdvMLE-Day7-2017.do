* ICPSR Advanced MLE -- 2017
*
* Day Seven Stata code
*
* Note: Results are very slightly different than when using
* R due to differences in algorithms.

* Simulated logits

ssc install rnd, replace

set seed 7222009
set obs 500
gen Z = invnorm(uniform())
gen W = invnorm(uniform())
gen XB = 0.2+0.5*W-0.5*Z
gen PY = (exp(XB)) / (1+exp(XB))
rndbin 500 PY 1
rename xb Y
gen X = 0
replace X = 1 in 251/500
replace X=0 if Y==0

logit Y W Z X

* Pets...

insheet using "Pets.csv"

replace partyid="DNR" if partyid==""
encode female, gen(Female)
encode married, gen(Married)
encode partyid, gen(PartyID)
encode education, gen(Education)

logit petfamily i.Female i.Married i.PartyID i.Education

xi: logit petfamily I.Female*I.Married I.PartyID I.Education

tab petfamily Married if Female==1
tab petfamily Married if Female==2

* Firth:

firthlogit petfamily _*


********************
* Survival
*
* Sims:

clear

set obs 200
set seed 7222009
gen X = 0
replace X=1 in 101/200
gen Scale = exp(0 + 0.2*X)
rndwei 200 1.2 Scale
rename xw T
rndbin 200 0.2 1
rename xb C
replace C=0 if X==0

stset T, failure(C)

stcox X, nohr 

streg X, nohr dist(weibull)

* The Heinze et al. correction has not yet been implemented in Stata...
* for this reason, I'm skipping the LHR analysis here.

