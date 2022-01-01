use NLSY97/NLSY97.dta, clear

*** Highest grade of the previous generation
gen HGC_FAM = max(R1302400, R1302500, R1302600, R1302700)
gen HDC_FAM = max(Z0502200, Z0502300, Z0502500, Z0502800)
label value HGC_FAM vlR1302400
label value HDC_FAM vlZ0502200

drop if HDC_FAM < 0
drop if Z9083900 < 0
drop if Z9083800 < 0

*** Generate the probability of getting a high-paid job
/*High-paid jobs are those paying at the 75th percentile while Low-paid jobs are those paying at 25th percentile*/
gen prob1 = 0
su U2857200 if U2857200 > 0, det
scalar lower25 = r(p25)
scalar upper25 = r(p75)
replace prob1 = 1 if U2857200 >= upper25
bys HDC_FAM: su prob1 if (U2857200 > 0 & U2857200 <= lower25) | U2857200 >= upper25 /*HDC: group with lowest probability and at least 20 observations is 2th (0.25), highest is 8th (0.816)*/
su prob1 if HDC_FAM == 2 & ((U2857200 > 0 & U2857200 <= lower25) | U2857200 >= upper25)
scalar minp1 = r(mean)
su prob1 if HDC_FAM == 8 & ((U2857200 > 0 & U2857200 <= lower25) | U2857200 >= upper25)
scalar maxp1 = r(mean)
scalar maxh1 = maxp1/minp1

su U2857200 if U2857200 >= upper25
scalar H1 = r(mean) /*H in the model as average of 25% highest-paid job in 2017*/
su U2857200  if U2857200 <= lower25 & U2857200 > 0
scalar L1 = r(mean) /*L in the model as average of 25% lowest-paid job in 2017*/



*** Find the amount of investment made
/*Investment = total amount of tuition fees + personal and federal loans over the years*/
foreach var in S2473400 S2473500 S4214800 S4214900 S5807000 S5807100 S7895900 S7896000 T0308700 T0308800 T2325300 T2325400 T3851600 T3851700 T5477300 T5477400 T6889500 T6889600 T8384400 T8384500 U0250800 U0250900 U2168900 U2169000 Z0414700 Z0415300 Z0519100 Z0519700 {
    replace `var' = . if `var' < 0
}

egen INV = rowtotal(S2473400 S2473500 S4214800 S4214900 S5807000 S5807100 S7895900 S7896000 T0308700 T0308800 T2325300 T2325400 T3851600 T3851700 T5477300 T5477400 T6889500 T6889600 T8384400 T8384500 U0250800 U0250900 U2168900 U2169000 Z0414700 Z0415300 Z0519100 Z0519700)

gen YEARLY_INV = INV/Z9083800 /*Dividing the total investment by the number of year of schooling (treating highest grade achieved as years of schooling)*/
su YEARLY_INV if YEARLY_INV > 0 & ((U2857200 > 0 & U2857200 <= lower25) | U2857200 >= upper25), det
scalar invrange1 = r(p95) - r(min) /*take the 95th percentile rather than maximum as the upper limit to exclude outliers*/

***Determine the coefficient of intergenerational linkage 
gen lnFAM = log(HDC_FAM)
gen lnSELF = log(Z9083900)
reghdfe lnSELF lnFAM if (U2857200 > 0 & U2857200 <= lower25) | U2857200 >= upper25, a(birth_year R1482600) /*to estimate alpha*/
scalar alpha1 = _b[lnFAM]

*Output summary table 1
estpost tabstat Z9083900 HDC_FAM YEARLY_INV U2857200 prob1 if (U2857200 > 0 & U2857200 <= lower25) | U2857200 >= upper25, stat(n mean sd min max) c(stat) by(prob1)
esttab, cells("count mean(fmt(%13.2fc)) sd(fmt(%13.2fc)) min max") nonumber nomtitle nonote noobs label collabels("Observations" "Mean" "SD" "Min" "Max")
*esttab using "summary1.tex", replace cells("count mean(fmt(%13.2fc)) sd(fmt(%13.2fc)) min max") nonumber nomtitle nonote noobs label booktabs f collabels("Observations" "Mean" "SD" "Min" "Max")


*************************************************************
************************ CALIBRATION ************************
*************************************************************
preserve
clear 
set obs 10000

****************************************
/*Discontinuity of investment decision*/
****************************************
scalar gamma1 = minp1
scalar barI1 = 2*maxp1/gamma1-2*(maxh1)^alpha1 /*associate the highest probability of success with maximum amount of investment and highest level of parental education*/
if barI1 > (maxp1/minp1 - 1) {
	scalar barI1 = maxh1 - 1
}
scalar scaled_H1 = H1*barI1/invrange1
scalar scaled_L1 = L1*barI1/invrange1
scalar lambda1 = 0.23
scalar threshold1 = 2.4317 /*estimated via Matlab*/
/*Matlab code
*pre-insert values of available parameters
func1 = @(h) -exp(-lambda1*scaled_L1)*(exp(lambda1*(maxh1-h^alpha1))-1)+gamma1/2*(exp(-lambda1*scaled_L1)-exp(-lambda1*scaled_H1))*((exp(lambda1*(maxh1-h^alpha1))*(maxh1+h^alpha1))-2*h^alpha1);
roots1 = fzero(func1, 2);
roots1
*/
scalar list lambda1 alpha1 gamma1 barI1 scaled_H1 scaled_L1 maxp1 minp1 maxh1 threshold1

*Expected utility
forvalues i = 1/6 {
	scalar h`i' = 1 + 2.1*`i'/6
	gen I`i' = (_n-1)*(maxh1 - h`i'^alpha1)/(_N-1)
	gen U`i' = -exp(-lambda1*(scaled_L1-I`i'))+(gamma1/2)*(2*h`i'^alpha1+I`i')*(exp(-lambda1*(scaled_L1-I`i'))-exp(-lambda1*(scaled_H1-I`i')))
	gen p`i' = gamma1/2*(2*h`i'^alpha1+I`i')	
}
forvalues i = 1/6 {
	gen U`i'_re = U`i'+ 1 + (U6[1]-U`i'[1])*(94+`i')/100 /*add a constant to make it easier to graph*/
}

*Expected wealth
gen h = 1 + (maxh1-1)*_n/_N
gen W = L1 + gamma1*h^alpha1*(H1-L1)
replace W = L1 + gamma1/2*(maxh1 + h^alpha1)*(H1-L1) if h > threshold1

*Graph 1A(a)
line U2_re I2, lpattern(-.) legend(col(1) symx(6) ring(0) position(8) label(1 "`=ustrunescape("h\u0302")'=1.7") label(2 "`=ustrunescape("h\u0302")'=2.05") label(3 "`=ustrunescape("h\u0302")'=2.4") label(4 "`=ustrunescape("h\u0302")'=2.75") label(5 "`=ustrunescape("h\u0302")'=3.1")) xtitle("I") lcolor(blue) graphregion(color(white)) xscale(titlegap(3)) lwidth(medthick) || line U3_re I3, lpattern(_) lcolor(blue) lwidth(medthick) || line U4_re I4, lpattern(-) lcolor(blue) lwidth(medthick) || line U5_re I5, lpattern(._) lcolor(red) lwidth(medthick) || line U6_re I6, lpattern(l) lcolor(red) lwidth(medthick)
*graph save "Figure1A(a)", replace
*graph export "Figure1A(a).png", as(png) replace

*Graph 1A(b)
line U2_re p2, lpattern(-.) legend(col(1) symx(6) ring(0) position(8) label(1 "`=ustrunescape("h\u0302")'=1.7") label(2 "`=ustrunescape("h\u0302")'=2.05") label(3 "`=ustrunescape("h\u0302")'=2.4") label(4 "`=ustrunescape("h\u0302")'=2.75") label(5 "`=ustrunescape("h\u0302")'=3.1")) xtitle("Probability") xscale(titlegap(3)) lcolor(blue) graphregion(color(white)) lwidth(medthick) || line U3_re p3, lpattern(_) lcolor(blue) lwidth(medthick) || line U4_re p4, lpattern(-) lcolor(blue) lwidth(medthick) || line U5_re p5, lpattern(._) lcolor(red) lwidth(medthick) || line U6_re p6, lpattern(l) lcolor(red) lwidth(medthick)
*graph save "Figure1A(b)", replace
*graph export "Figure1A(b).png", as(png) replace

*Graph 1B
line W h if h < threshold1, xlabel(`=maxh1' 1(0.5)3, ax(1) format(%13.1gc)) xlabel(`=maxh1' 1(0.5)3, ax(2) format(%13.1gc)) xscale(r(1 3.3) ax(1)) xscale(r(1 3.3) ax(2)) lwidth(thick) lcolor(blue) xline(`=threshold1', lp(-)) ytitle("Expected Wealth (US$)") yscale(titlegap(2)) xtitle("`=ustrunescape("h\u0302")'") xmlabel(`=threshold1', axis(1) format(%13.1fc) angle(45)) xmlabel(`=threshold1', axis(2) format(%13.1fc) angle(-45)) ytitle("", axis(2)) xtitle("", axis(2)) graphregion(color(white)) legend(off) xaxis(1 2) yaxis(1 2) xscale(titlegap(3)) || line W h if h > threshold1, lcolor(blue) lwidth(thick) 
*graph save "Figure1B", replace
*graph export "Figure1B.png", as(png) replace

restore


****************************
/*The evolution of bequest*/
****************************
*** Restrict the sample now to 15% upper and 15% lower
centile U2857200 if U2857200 > 0, centile(15 85)
scalar lower15 = r(c_1)
scalar upper15 = r(c_2)
su U2857200 if U2857200 >= upper15
scalar H2 = r(mean) /*H in the model as average of 15% highest-paid job in 2017*/
su U2857200  if U2857200 <= lower15 & U2857200 > 0
scalar L2 = r(mean) /*L in the model as average of 15% lowest-paid job in 2017*/
gen prob2 = 0
replace prob2 = 1 if U2857200 >= upper15
bys HDC_FAM: su prob2 if (U2857200 > 0 & U2857200 <= lower15) | U2857200 >= upper15 
/*HDC: upper and lower 15: lowest is 2nd (0.273) highest is 8th (0.786)*/

su prob2 if HDC_FAM == 1 & ((U2857200 > 0 & U2857200 <= lower15) | U2857200 >= upper15)
scalar minp2 = r(mean)
su prob2 if HDC_FAM == 8 & ((U2857200 > 0 & U2857200 <= lower15) | U2857200 >= upper15)
scalar maxp2 = r(mean)
scalar maxh2 = maxp2/minp2

*Output summary table 2
estpost tabstat Z9083900 HDC_FAM YEARLY_INV U2857200 prob2 if (U2857200 > 0 & U2857200 <= lower15) | U2857200 >= upper15, stat(n mean sd min max) c(stat) by(prob2)
esttab, cells("count mean(fmt(%13.2fc)) sd(fmt(%13.2fc)) min max") nonumber nomtitle nonote noobs label collabels("Observations" "Mean" "SD" "Min" "Max")
*esttab using "summary2.tex", replace cells("count mean(fmt(%13.2fc)) sd(fmt(%13.2fc)) min max") nonumber nomtitle nonote noobs label booktabs f collabels("Observations" "Mean" "SD" "Min" "Max")

reghdfe lnSELF lnFAM if (U2857200 > 0 & U2857200 <= lower15) | U2857200 >= upper15, a(birth_year R1482600)
scalar alpha2 = _b[lnFAM]
su YEARLY_INV if YEARLY_INV > 0 & ((U2857200 > 0 & U2857200 <= lower15) | U2857200 >= upper15), det
scalar invrange2 = r(p95) - r(min)


***Calibrated parameters (separating equilibrium)
preserve
clear
set obs 10
scalar R = 1.051
scalar b0 = 1
scalar beta = 0.3
scalar lambda2 = 0.24 /*pooling equilibrium for lambda >= 0.25*/
scalar gamma2 = minp2
scalar barI2 = 2*maxp2/gamma2-2*(maxh2)^alpha2 /*associate the highest probability of success with maximum amount of investment and highest level of parental education*/
if barI2 > (maxp2/minp2 - 1) {
	scalar barI2 = maxh2 - 1
}
scalar scaled_H2 = H2*barI2/invrange2
scalar scaled_L2 = L2*barI2/invrange2
scalar threshold2 = 1.9747 /*estimated via Matlab, similar code as above*/
scalar list lambda2 alpha2 gamma2 barI2 scaled_H2 scaled_L2 maxp2 minp2 maxh2 threshold2 R beta 

gen t = _n-1
tsset t

*Household 1: invest in education
scalar h01 = 2
gen h1 = h01 if t == 0
gen b1 = b0 if t == 0
gen I1 = 0 if t == 0

forvalues i = 1/`=_N' {
	scalar temp = -exp(-lambda2*scaled_L2)*(exp(lambda2*(maxh2-h1^alpha2)*R)-1)+gamma2/2*(exp(-lambda2*scaled_L2)-exp(-lambda2*scaled_H2))*((exp(lambda2*(maxh2-h1^alpha2)*R)*(maxh2+h1^alpha2))-2*h1^alpha2)
	replace I1 = maxh2-L.h1^alpha2 if L.h1 > threshold2 & t > 0 & temp > 0
	replace h1 = (maxh2 + L.h1^alpha2)/2  if L.h1 > threshold2 & t > 0 & temp > 0
	replace b1 = beta*((L.b1-maxh2+L.h1^alpha2)*R + scaled_L2 + (scaled_H2-scaled_L2)*gamma2*h1)  if t > 0 & temp > 0

	replace h1 = L.h1^alpha2  if (L.h1 < threshold2 | temp < 0) & t > 0 
	replace b1 = beta*(L.b1*R + scaled_L2 + (scaled_H2-scaled_L2)*gamma2*h1)  if (L.h1 < threshold2 | temp < 0) & t > 0 
	replace I1 = 0 if (L.h1 < threshold2 | temp < 0) & t > 0 
}

*Household 2: invest in the capital market
scalar h02 = 1.9
gen h2 = h02 if t == 0
gen b2 = b0 if t == 0
replace h2 = L.h2^alpha2 if t > 0
replace b2 = beta*(L.b2*R + scaled_L2 + (scaled_H2-scaled_L2)*gamma2*h2) if t > 0
gen I2 = 0

*Graph 2(a)
line b1 t, xlabel(0(1)9) graphregion(color(white)) ytitle("b{subscript:t}", size(large)) ylabel(,labsize(11pt)) xtitle(,size(large)) legend(col(1) symx(15) ring(0) position(4) label(1 "b{subscript:1,t}") label(2 "b{subscript:2,t}") size(18pt)) xlabel(,labsize(11pt)) lwidth(medthick) xscale(titlegap(3)) lcolor(blue) || line b2 t, lpattern(-) lwidth(thick) lcolor(black)
*graph save "Figure2(a)", replace
*graph export "Figure2(a).png", as(png) replace

*Graph 2(b)
line h1 t, xlabel(0(1)9) graphregion(color(white)) yline(`=threshold2', lp(_)) ttext(`=threshold2-0.08' 8 "`=ustrunescape("h\u0303")'=1.975", size(large)) ytitle("h{subscript:t}", size(large)) xtitle(,size(large)) ylabel(, labsize(11pt)) xlabel(, labsize(11pt)) legend(col(1) symx(15) ring(0) position(3) label(1 "h{subscript:1,t}") label(2 "h{subscript:2,t}") size(18pt)) lwidth(medthick) xscale(titlegap(3)) lcolor(blue) || line h2 t, lpattern(-) lwidth(thick) lcolor(black)
*graph save "Figure2(b)", replace
*graph export "Figure2(b).png", as(png) replace

*Graph 2(c)
line I1 t, xlabel(0(1)9, labsize(11pt)) ylabel(,labsize(11pt)) graphregion(color(white)) yline(`=barI2', lp(_)) ytitle("I{subscript:t}", size(large)) lwidth(medthick) xscale(titlegap(3)) yscale(titlegap(2)) xtitle(,size(large)) legend(col(1) symx(15) ring(0) position(3) label(1 "I{subscript:1,t}") label(2 "I{subscript:2,t}") size(18pt)) lcolor(blue) || line I2 t, lp(-) lwidth(thick) lcolor(black)
*graph save "Figure2(c)", replace
*graph export "Figure2(c).png", as(png) replace

restore

***Calibrated parameters (pooling equilibrium)
preserve
clear
set obs 10
scalar lambda3 = 0.246 /*pooling equilibrium for lambda >= 0.25*/
scalar threshold3 = 2.1658 

gen t = _n-1
tsset t

*Household 1: invest in education
scalar h01 = 2.5
gen h1 = h01 if t == 0
gen b1 = b0 if t == 0
gen I1 = 0 if t == 0

forvalues i = 1/`=_N' {
	scalar temp = -exp(-lambda3*scaled_L2)*(exp(lambda3*(maxh2-h1^alpha2)*R)-1)+gamma2/2*(exp(-lambda3*scaled_L2)-exp(-lambda3*scaled_H2))*((exp(lambda3*(maxh2-h1^alpha2)*R)*(maxh2+h1^alpha2))-2*h1^alpha2)
	replace I1 = maxh2-L.h1^alpha2 if L.h1 > threshold3 & t > 0 & temp > 0
	replace h1 = (maxh2 + L.h1^alpha2)/2  if L.h1 > threshold3 & t > 0 & temp > 0
	replace b1 = beta*((L.b1-maxh2+L.h1^alpha2)*R + scaled_L2 + (scaled_H2-scaled_L2)*gamma2*h1)  if t > 0 & temp > 0

	replace h1 = L.h1^alpha2  if (L.h1 < threshold3 | temp < 0) & t > 0 
	replace b1 = beta*(L.b1*R + scaled_L2 + (scaled_H2-scaled_L2)*gamma2*h1)  if (L.h1 < threshold3 | temp < 0) & t > 0 
	replace I1 = 0 if (L.h1 < threshold3 | temp < 0) & t > 0 
}

*Household 2: invest in the capital market
scalar h02 = 1.9
gen h2 = h02 if t == 0
gen b2 = b0 if t == 0
replace h2 = L.h2^alpha2 if t > 0
replace b2 = beta*(L.b2*R + scaled_L2 + (scaled_H2-scaled_L2)*gamma2*h2) if t > 0
gen I2 = 0

*Graph 3(a)
line b1 t, xlabel(0(1)9) graphregion(color(white)) ytitle("b{subscript:t}", size(large)) ylabel(,labsize(11pt)) xtitle(,size(large)) legend(col(1) symx(15) ring(0) position(4) label(1 "b{subscript:1,t}") label(2 "b{subscript:2,t}") size(18pt)) xlabel(,labsize(11pt)) lwidth(medthick) xscale(titlegap(3)) lcolor(blue) || line b2 t, lpattern(-) lwidth(thick) lcolor(black)
*graph save "Figure3(a)", replace
*graph export "Figure3(a).png", as(png) replace

*Graph 3(b)
line h1 t, xlabel(0(1)9) graphregion(color(white)) yline(`=threshold3', lp(_)) ttext(`=threshold3+0.08' 8 "`=ustrunescape("h\u0303")'=2.165", size(large)) ytitle("h{subscript:t}", size(large)) xtitle(,size(large)) ylabel(, labsize(11pt)) xlabel(, labsize(11pt)) legend(col(1) symx(15) ring(0) position(3) label(1 "h{subscript:1,t}") label(2 "h{subscript:2,t}") size(18pt)) lwidth(medthick) xscale(titlegap(3)) lcolor(blue) || line h2 t, lpattern(-) lwidth(thick) lcolor(black)
*graph save "Figure3(b)", replace
*graph export "Figure3(b).png", as(png) replace

*Graph 3(c)
line I1 t, xlabel(0(1)9, labsize(11pt)) ylabel(,labsize(11pt)) graphregion(color(white)) yline(`=barI2', lp(_)) ytitle("I{subscript:t}", size(large)) lwidth(medthick) xscale(titlegap(3)) yscale(titlegap(2)) xtitle(,size(large)) legend(col(1) symx(15) ring(0) position(3) label(1 "I{subscript:1,t}") label(2 "I{subscript:2,t}") size(18pt)) lcolor(blue) || line I2 t, lp(-) lwidth(thick) lcolor(black)
*graph save "Figure3(c)", replace
*graph export "Figure3(c).png", as(png) replace

restore



*********************
/*CRRA calibration*/
*********************
preserve
clear 
set obs 500


/*Discontinuity of investment decision*/
scalar sigma1 = 1
scalar b0 = 1.4
scalar list alpha2 gamma2 scaled_H2 scaled_L2 maxp2 minp2 maxh2 sigma1 R b0

*Expected utility
forvalues i = 1/6 {
	scalar h`i' = 1 + (maxh2 -1)*`i'/6
	gen I`i' = (_n-1)*(maxh2 - h`i'^alpha2)/(_N-1)
	gen p`i' = gamma2/2*(2*h`i'^alpha2+I`i')	
	gen U`i' = p`i'*((scaled_H2+(b0-I`i')*R)^(1-sigma1)-1)/(1-sigma1)+(1-p`i')*((scaled_L2+(b0-I`i')*R)^(1-sigma1)-1)/(1-sigma1) if sigma1 != 1
	replace U`i' = p`i'*log(scaled_H2+(b0-I`i')*R)+(1-p`i')*log(scaled_L2+(b0-I`i')*R) if sigma1 == 1
}

forvalues i = 1/6 {
	gen U`i'_re = U`i' + 1 + (U6[1]-U`i'[1])*(85+`i')/100 /*add a constant to make it easier to graph*/
	su U`i'_re
	if `i' == 6 {
		scalar maxU = r(max)
	}
}

gen h = 1 + (maxh2-1)*(_n-1)/(_N-1)
gen barI = maxh2 - h^alpha2
gen tempI = .

/*WARNING: THIS NEXT PIECE OF MATA CODE USES ESTIMATION TECHNIQUE WITH 5,000 ITERATIONS IN TOTAL (RUN TIME ~5 MINS)*/
clear mata
mata

gamma2 = st_numscalar("gamma2")
scaled_L2 = st_numscalar("scaled_L2")
scaled_H2 = st_numscalar("scaled_H2")
alpha2 = st_numscalar("alpha2")
R = st_numscalar("R")
b0 = st_numscalar("b0")
sigma1 = st_numscalar("sigma1")
maxh2 = st_numscalar("maxh2")

void myeval(todo, tempI, real scalar i, U, g, H)
{
	external gamma2, scaled_L2, scaled_H2, alpha2, R, b0, sigma1, maxh2
	real scalar temph
	temph = st_data(i, "h")
	real scalar I
	I = (maxh2 - temph^alpha2)*invlogit(tempI)
	real scalar p 
	p = gamma2/2*(2*temph^alpha2+I)
	if (sigma1 != 1) {
		U = p*((scaled_H2+(b0-I)*R)^(1-sigma1)-1)/(1-sigma1)+(1-p)*((scaled_L2+(b0-I)*R)^(1-sigma1)-1)/(1-sigma1)
	}
	else {
		U = p*log(scaled_H2+(b0-I)*R)+(1-p)*log(scaled_L2+(b0-I)*R)
	}
}

for (i=1;i<=st_nobs();i++) {
	S = optimize_init()
	optimize_init_evaluator(S, &myeval())
	optimize_init_which(S, "max")
	optimize_init_params(S,0)
	optimize_init_conv_maxiter(S, 10)
	optimize_init_conv_warning(S, "off")
	optimize_init_tracelevel(S, "none")
	optimize_init_argument(S,1,i)
	tempI = optimize(S)
	st_store(i, "tempI", tempI)
}
end
/*ENDS MATA CODE*/

gen I = (maxh2 - h^alpha2)*invlogit(tempI)
replace I = 0 if I < 1e-6
su h if I > 0
scalar threshold4 = r(min)
gen realI = I*invrange2/barI2
scalar realb0 = b0*invrange2/barI2
gen b = beta*(L2 + (realb0-realI)*R + gamma2/2*(2*h^alpha2+I)*(H2-L2))
gen b1 = b if I == 0
replace b1 = beta*(L2 + (realb0)*R + gamma2/2*(2*h^alpha2)*(H2-L2)) if I > 0
gen p = gamma2/2*(2*h^alpha2+I)


*Graph C1(a)
line U1_re I1, lpattern(-..-) lcolor(blue) lwidth(medthick) ylabel(`=maxU' 2.75(0.05)2.85, format(%13.3gc)) xlabel(`=maxh2-h1^alpha2' 0(0.3)1.5, format(%13.1gc)) || line U2_re I2, lpattern(-.) legend(col(1) symx(8) ring(0) position(8) label(1 "`=ustrunescape("h\u0302")'=1.3") label(2 "`=ustrunescape("h\u0302")'=1.6") label(3 "`=ustrunescape("h\u0302")'=1.9") label(4 "`=ustrunescape("h\u0302")'=2.2") label(5 "`=ustrunescape("h\u0302")'=2.5") label(6 "`=ustrunescape("h\u0302")'=2.86")) xtitle("I") lcolor(blue) graphregion(color(white)) xscale(titlegap(3)) lwidth(medthick) || line U3_re I3, lpattern(_) lcolor(blue) lwidth(medthick) || line U4_re I4, lpattern(-) lcolor(red) lwidth(medthick) || line U5_re I5, lpattern(._) lcolor(red) lwidth(medthick) || line U6_re I6, lpattern(l) lcolor(red) lwidth(medthick)
*graph save "FigureD1(a)", replace
*graph export "FigureD1(a).png", as(png) replace


*Graph C1(b)
line U1_re p1, lpattern(-..-) lcolor(blue) lwidth(medthick) ylabel(`=maxU' 2.75(0.05)2.85, format(%13.3gc)) || line U2_re p2, lpattern(-.) legend(col(1) symx(8) ring(0) position(8) label(1 "`=ustrunescape("h\u0302")'=1.3") label(2 "`=ustrunescape("h\u0302")'=1.6") label(3 "`=ustrunescape("h\u0302")'=1.9") label(4 "`=ustrunescape("h\u0302")'=2.2") label(5 "`=ustrunescape("h\u0302")'=2.5") label(6 "`=ustrunescape("h\u0302")'=2.86")) xtitle("Probability") xscale(titlegap(3)) lcolor(blue) graphregion(color(white)) lwidth(medthick) || line U3_re p3, lpattern(_) lcolor(blue) lwidth(medthick) || line U4_re p4, lpattern(-) lcolor(red) lwidth(medthick) || line U5_re p5, lpattern(._) lcolor(red) lwidth(medthick) || line U6_re p6, lpattern(l) lcolor(red) lwidth(medthick)
*graph save "FigureD1(b)", replace
*graph export "FigureD1(b).png", as(png) replace

*Graph C2(a)
line b1 h, color(black) lp(_.) ytitle("Expected Bequest (US$)") xline(`=threshold4', lp(-)) xmlabel(`=threshold4', ax(1) angle(45) format(%13.3gc)) xmlabel(`=threshold4', ax(2) angle(-45) format(%13.3gc)) yscale(titlegap(2)) xtitle("`=ustrunescape("h\u0302")'") ytitle("", axis(2)) xtitle("", axis(2)) graphregion(color(white)) legend(off) xaxis(1 2) yaxis(1 2) xscale(titlegap(3)) xlabel(`=maxh2' 1(0.5)2.5, ax(1) format(%13.1gc)) xlabel(`=maxh2' 1(0.5)2.5, ax(2) format(%13.1gc)) || line b h, lwidth(thick) lcolor(blue)   
*graph save "FigureD2(a)", replace
*graph export "FigureD2(a).png", as(png) replace


*Graph C2(b)
line realI h, lwidth(thick) lcolor(blue) xline(`=threshold4', lp(-)) xmlabel(`=threshold4', ax(1) angle(45) format(%13.3gc)) xmlabel(`=threshold4', ax(2) angle(-45) format(%13.3gc)) ytitle("Expected Investment (US$)") yscale(titlegap(2)) xtitle("`=ustrunescape("h\u0302")'")  ytitle("", axis(2)) xtitle("", axis(2)) graphregion(color(white)) legend(off) xaxis(1 2) yaxis(1 2) xscale(titlegap(3)) xlabel(`=maxh2' 1(0.5)2.5, ax(1) format(%13.1gc)) xlabel(`=maxh2' 1(0.5)2.5, ax(2) format(%13.1gc))
*graph save "FigureD2(b)", replace
*graph export "FigureD2(b).png", as(png) replace

restore

/*Intergenerational evolution of bequest and human capital*/
preserve
clear
set obs 10
gen t = _n-1

gen h1 = 2 if t == 0 	/*baseline h & b, average household*/
gen b1 = 1.4 if t == 0 	
gen h2 = 2.4 if t == 0 	/*poorer but more educated*/
gen b2 = 1.2 if t == 0 
gen h3 = 1.5 if t == 0 	/*richer but less educated*/
gen b3 = 1.7 if t == 0 
gen h4 = 2.8 if t == 0 /*more educated but very poor*/
gen b4 = 0.3 if t == 0
gen h5 = 2.8 if t == 0
gen b5 = 0.4 if t == 0

forvalues i = 1/5 {
	gen I`i' = .
}

clear mata
mata

gamma2 = st_numscalar("gamma2")
scaled_L2 = st_numscalar("scaled_L2")
scaled_H2 = st_numscalar("scaled_H2")
alpha2 = st_numscalar("alpha2")
R = st_numscalar("R")
sigma1 = st_numscalar("sigma1")
maxh2 = st_numscalar("maxh2")
beta = 0.2

void myeval(todo, tempI, U, g, H)
{
	external gamma2, scaled_L2, scaled_H2, alpha2, R, sigma1, maxh2, h, b, beta
	real scalar I
	I = (maxh2 - h^alpha2)*invlogit(tempI)
	real scalar p 
	p = gamma2/2*(2*h^alpha2+I)
	if (sigma1 != 1) {
		U = p*((scaled_H2+(b-I)*R)^(1-sigma1)-1)/(1-sigma1)+(1-p)*((scaled_L2+(b-I)*R)^(1-sigma1)-1)/(1-sigma1)
	}
	else {
		U = p*log(scaled_H2+(b-I)*R)+(1-p)*log(scaled_L2+(b-I)*R)
	}
}

void routine(real scalar i, real scalar j)
{ /*j is row, i is column*/
	string v1
	v1 = sprintf("h%g", i)
	string v2
	v2 = sprintf("b%g", i)
	string v3
	v3 = sprintf("I%g", i)
	scalar h
	h = st_data(j, v1)
	scalar b
	b = st_data(j, v2)
	
	pointer() scalar temph
	temph = crexternal("h")
	*temph = h
	pointer() scalar tempb
	tempb = crexternal("b")
	*tempb = b
	
	external maxh2, alpha2, myeval(), beta, gamma2, scaled_H2, scaled_L2, R
	
	transmorphic S
	S = optimize_init()
	optimize_init_evaluator(S, &myeval())
	optimize_init_which(S, "max")
	optimize_init_params(S,0)
	optimize_init_conv_maxiter(S,10)
	optimize_init_conv_warning(S, "off")
	optimize_init_tracelevel(S, "none")
	scalar tempI
	tempI = optimize(S)
	scalar I
	I = (maxh2 - h^alpha2)*invlogit(tempI)
	j=j+1
	st_store(j, v3, I)	
	scalar outh
	if (I != 0) {
		outh = h^alpha2 + I/2
	}
	else {
		outh = h^alpha2
	}
	st_store(j, v1, outh)
	scalar outb
	outb = beta*(scaled_L2 + (b-I)*R + gamma2/2*(2*h^alpha2+I)*(scaled_H2-scaled_L2))
	st_store(j, v2, outb)
	
	rmexternal("b")
	rmexternal("h")
}

end

/*Optimize for each entry of each household*/
forvalues i = 1/5 {
	forvalues j = 1/9 {
		mata: routine(`i',`j')
	}
}


*Graph C3(a)
line b1 t, lpattern(l) lcolor(red) lwidth(medthick) xlabel(0(1)9) || line b2 t, lpattern(-) legend(col(1) symx(8) ring(0) position(4) label(1 "b{subscript:1,0}= `: di %-4.1f `=b1'', h{subscript:1,0}=`: di %-4.1f `=h1''") label(2 "b{subscript:2,0}= `: di %-4.1f `=b2'', h{subscript:2,0}=`: di %-4.1f `=h2''") label(3 "b{subscript:3,0}= `: di %-4.1g `=b3'', h{subscript:3,0}=`: di %-4.1f `=h3''") label(4 "b{subscript:4,0}= `: di %-4.1f `=b4'', h{subscript:4,0}=`: di %-4.1f `=h4''") label(5 "b{subscript:5,0}= `: di %-4.1f `=b5'', h{subscript:5,0}=`: di %-4.1f `=h5''")) yscale(titlegap(3)) ytitle("b{subscript:t}") xscale(titlegap(3)) lcolor(red) graphregion(color(white)) lwidth(medthick) || line b3 t, lpattern(_) lcolor(blue) lwidth(medthick) || line b4 t, lpattern(-.) lcolor(blue) lwidth(medthick) || line b5 t, lp(l) lcolor(black) lwidth(medthick)
*graph save "FigureD3(a)", replace
*graph export "FigureD3(a).png", replace as(png)

*Graph C3(b)
line h1 t, lpattern(l) lcolor(red) lwidth(medthick) xlabel(0(1)9) || line h2 t, lpattern(-) legend(col(1) symx(8) ring(0) position(2) label(1 "b{subscript:1,0}= `: di %-4.1f `=b1'', h{subscript:1,0}=`: di %-4.1f `=h1''") label(2 "b{subscript:2,0}= `: di %-4.1f `=b2'', h{subscript:2,0}=`: di %-4.1f `=h2''") label(3 "b{subscript:3,0}= `: di %-4.1g `=b3'', h{subscript:3,0}=`: di %-4.1f `=h3''") label(4 "b{subscript:4,0}= `: di %-4.1f `=b4'', h{subscript:4,0}=`: di %-4.1f `=h4''") label(5 "b{subscript:5,0}= `: di %-4.1f `=b5'', h{subscript:5,0}=`: di %-4.1f `=h5''")) yscale(titlegap(3)) ytitle("h{subscript:t}") xscale(titlegap(3)) lcolor(red) graphregion(color(white)) lwidth(medthick) || line h3 t, lpattern(_) lcolor(blue) lwidth(medthick) || line h4 t, lpattern(-.) lcolor(blue) lwidth(medthick) || line h5 t, lp(l) lcolor(black) lwidth(medthick)
*graph save "FigureD3(b)", replace
*graph export "FigureD3(b).png", replace as(png)

*Graph C3(c)
line I1 t, lpattern(l) lcolor(red) lwidth(medthick) xlabel(0(1)9) || line I2 t, lpattern(-) legend(col(1) symx(8) ring(0) position(3) label(1 "b{subscript:1,0}= `: di %-4.1f `=b1'', h{subscript:1,0}=`: di %-4.1f `=h1''") label(2 "b{subscript:2,0}= `: di %-4.1f `=b2'', h{subscript:2,0}=`: di %-4.1f `=h2''") label(3 "b{subscript:3,0}= `: di %-4.1g `=b3'', h{subscript:3,0}=`: di %-4.1f `=h3''") label(4 "b{subscript:4,0}= `: di %-4.1f `=b4'', h{subscript:4,0}=`: di %-4.1f `=h4''") label(5 "b{subscript:5,0}= `: di %-4.1f `=b5'', h{subscript:5,0}=`: di %-4.1f `=h5''")) yscale(titlegap(3)) ytitle("I{subscript:t}") xscale(titlegap(3)) lcolor(red) graphregion(color(white)) lwidth(medthick) || line I3 t, lpattern(_) lcolor(blue) lwidth(medthick) || line I4 t, lpattern(-.) lcolor(blue) lwidth(medthick) || line I5 t, lp(l) lcolor(black) lwidth(medthick)
*graph save "FigureD3(c)", replace
*graph export "FigureD3(c).png", replace as(png)

restore