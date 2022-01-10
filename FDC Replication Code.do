*****************************************************************************************************	
***************************************************************************************************** 
*   Replication Code for 
*
*   “One pill, once a day”: the effect of simplified treatment regimens on retention in HIV care
*   Kluberg, Bor, LaValley, Evans, Hirasen, Maskew, Long, Fox
*
*****************************************************************************************************
*****************************************************************************************************


*** Notes on RDD specification used in all analyses
	* See discussion in https://cattaneo.princeton.edu/books/Cattaneo-Idrobo-Titiunik_2018_CUP-Vol1.pdf
	* We use the rdrobust defaults: local linear regression with a triangular kernel, MSE-optimal BW w regularization parameter to improve stability when bias is small, and robust bias-corrected CI with the same MSE-optimal bandwidth; in sensitivity analyses, use robust bias-corrected CI using coverage-error-rate CER optimal bandwidth.
	* outputs include: control mean, treated mean, MSE-optimal point estimate; bias-corrected point estimate; robust bias-corrected CI using MSE-optimal bandwidth; robust bias-corrected CI using CER-optimal bandwidth


	
	

*** 0. Preliminaries

* clear workspace
clear

* Install rdrobust command
ssc install rdrobust

* Set user directory
cd "/Users/jbor/Dropbox/Kluberg_FDC"

* Read in data
use "FDC_final.dta"

	

	
	
*** 1. Cohort Characteristics

* a) Total Sample (described in text)
count


* b) Descriptive characteristics of sample within +/- 180 days and at threshold (Table 1)

* -180 to 180 days
su male age_at_initiation cd4_baseline3 stage_high anemia fdc if abs(date_init_ctr) <=180
tab1 male stage_high anemia fdc if abs(date_init_ctr) <=180

* -180 to -1 days
su male age_at_initiation cd4_baseline3 stage_high anemia fdc if date_init_ctr >=-180 & date_init_ctr<0
tab1 male stage_high anemia fdc if date_init_ctr >=-180 & date_init_ctr<0

* 0 to 180 days
su male age_at_initiation cd4_baseline3 stage_high anemia fdc if date_init_ctr >=0 & date_init_ctr<=180
tab1 male stage_high anemia fdc if date_init_ctr >=0 & date_init_ctr<=180

* Differences at threshold, with MSE-optimal bw
foreach v of varlist male age_at_initiation cd4_baseline3 stage_high anemia fdc {
	rdrobust `v' date_init_ctr , p(1) kernel(triangular) bwselect(mserd) all
	}

* c) Regimen type (described in text)
*regimen if FDC
tab first_regimen if fdc==1
*regimen if not FDC
tab first_regimen if fdc==0

	
 
*** 2. Main RDD analysis: ITT and CACE results 

gen T = date_init_ctr>=0
forvalues o = 1/4 {
	
	* a) ITT MSE-optimal BW, point estimate, and robust bias-adj CI (table 2)
		rdrobust outcome`o' date_init_ctr , p(1) kernel(triangular) bwselect(mserd) all
		local ITT_h = e(h_l)
		
		* ITT Model re-run by hand to get control and treated predictions (table 2)
		cap drop ITT_w
		gen ITT_w = .
		replace ITT_w = (1 - abs(date_init_ctr / `ITT_h')) if date_init_ctr < 0 & date_init_ctr >= -`ITT_h'
		replace ITT_w = (1 - abs(date_init_ctr / `ITT_h')) if date_init_ctr >= 0 & date_init_ctr <= `ITT_h'
		regress outcome`o' c.date_init_ctr##i.T if abs(date_init_ctr) < `ITT_h' [aw=ITT_w], r 
		dis "ITT Control mean = " _b[_cons]
		dis "ITT Treated mean = " _b[_cons] + _b[1.T]
		
		* CACE MSE-optimal BW, point estimate, and robust bias-adj CI 
		rdrobust outcome`o' date_init_ctr , fuzzy(fdc) p(1) kernel(triangular)  bwselect(mserd)
		local CACE_h = e(h_l)
	
	* b) CACE Model re-run by hand to get ITT and First Stage based on CACE BW (appendix 7)
		cap drop CACE_w
		gen CACE_w = .
		replace CACE_w = (1 - abs(date_init_ctr / `CACE_h')) if date_init_ctr < 0 & date_init_ctr >= -`CACE_h'
		replace CACE_w = (1 - abs(date_init_ctr / `CACE_h')) if date_init_ctr >= 0 & date_init_ctr <= `CACE_h'
		regress fdc c.date_init_ctr##i.T if abs(date_init_ctr) < `CACE_h' [aw=CACE_w], r
		local b_FS = _b[1.T]
		regress outcome`o' c.date_init_ctr##i.T if abs(date_init_ctr) < `CACE_h' [aw=CACE_w], r
		local b_ITT = _b[1.T]
		dis "First Stage in fuzzy RDD = " `b_FS'
		dis "ITT in fuzzy RDD = " `b_ITT'
		dis "CACE by hand = " `b_ITT' / `b_FS'
	
	* c) CACE robustness checks (appendix 9)
		local twice_h = 2*`CACE_h'
		local half_h = 0.5*`CACE_h'
		rdrobust outcome`o' date_init_ctr , fuzzy(fdc) p(1) kernel(triangular) h(`half_h')
		rdrobust outcome`o' date_init_ctr , fuzzy(fdc) p(1) kernel(triangular) h(`twice_h')
	
	}


*** 3. Subgroup Analyses (appendix 8)

forvalues o = 1/1 {
local o = 1
	rdrobust outcome`o' date_init_ctr if male == 1, fuzzy(fdc) p(1) kernel(triangular)  bwselect(mserd)
	rdrobust outcome`o' date_init_ctr if male == 0, fuzzy(fdc) p(1) kernel(triangular)  bwselect(mserd)
	rdrobust outcome`o' date_init_ctr if age_at >= 16 & age_at <30, fuzzy(fdc) p(1) kernel(triangular)  bwselect(mserd)
	rdrobust outcome`o' date_init_ctr if age_at >= 30 & age_at <40, fuzzy(fdc) p(1) kernel(triangular)  bwselect(mserd)
	rdrobust outcome`o' date_init_ctr if age_at >= 40 & age_at <50, fuzzy(fdc) p(1) kernel(triangular)  bwselect(mserd)
	rdrobust outcome`o' date_init_ctr if age_at >= 50 & age_at <., fuzzy(fdc) p(1) kernel(triangular)  bwselect(mserd)
	rdrobust outcome`o' date_init_ctr if anemia == 1, fuzzy(fdc) p(1) kernel(triangular)  bwselect(mserd)
	rdrobust outcome`o' date_init_ctr if anemia == 0, fuzzy(fdc) p(1) kernel(triangular)  bwselect(mserd)
	rdrobust outcome`o' date_init_ctr if stage_high == 1, fuzzy(fdc) p(1) kernel(triangular)  bwselect(mserd)
	rdrobust outcome`o' date_init_ctr if stage_high == 0, fuzzy(fdc) p(1) kernel(triangular)  bwselect(mserd)
	rdrobust outcome`o' date_init_ctr if cd4_baseline >= 0 & cd4_baseline <=200, fuzzy(fdc) p(1) kernel(triangular)  bwselect(mserd)
	rdrobust outcome`o' date_init_ctr if cd4_baseline > 200 & cd4_baseline <., fuzzy(fdc) p(1) kernel(triangular)  bwselect(mserd)
	}

* interaction effects / test that effects are different
gen highCD4 = cd4_baseline>200
gen no_anemia = 1 - anemia
gen stage_low = 1 - stage_high
reg outcome1 c.date_init_ctr##i.post##i.highCD4 if abs(date_init_ctr) <150 & cd4_baseline <., r
reg outcome1 c.date_init_ctr##i.post##i.no_anemia if abs(date_init_ctr) <150, r
reg outcome1 c.date_init_ctr##i.post##i.stage_low if abs(date_init_ctr) <150, r
reg outcome1 c.date_init_ctr##i.post##i.male if abs(date_init_ctr) <150, r


	
*** 4. Figures

* set bandwidth to equal the MSE-optimal bandwidth in the ITT models
	
* FDC (First Stage) (fig 1)
	
	* create binned averages
	cap drop month_ctr
	gen month_ctr = 30*floor(date_init_ctr/30)+15
	bys month_ctr: egen mean_fdc = mean(fdc)
	bys month_ctr: replace mean_fdc = . if _n>1
	
	* plot
	rdrobust fdc date_init_ctr , p(1) kernel(triangular) bwselect(mserd) all
	local ITT_h = e(h_l)

	twoway (scatter mean_fdc month_ctr, mc(gs3) msymbol(Oh)) ///
		(lpoly fdc date_init_ctr if date_init_ctr <0, bw(`ITT_h') lw(med) lc(black) k(tri) degree(1)) ///
		(lpoly fdc date_init_ctr if date_init_ctr >0, bw(`ITT_h') lw(med) lc(black) k(tri) degree(1)), ///
		legend(off) graphregion(color(white)) ///
		xtitle(No. of Days Since Guideline Change, size(medlarge) margin(med)) ///
		ytitle(Proportion Starting FDC, size(medlarge) margin(med)) ///
		ylabel(0(.2)1, labsize(medlarge) format(%02.1f) angle(0) nogrid) ///
		xlabel(-547 -365 -182 0 182 365 547, labsize(medlarge)) legend(off) xline(0, lc(black) lp(dash) lw(thin))
				

				
				
* Outcome 1 (fig 2)
	
	* create binned averages
	bys month_ctr: egen mean_outcome1 = mean(outcome1)
	bys month_ctr: replace mean_outcome1 = . if _n>1
	
	* plot
	rdrobust outcome1 date_init_ctr , p(1) kernel(triangular) bwselect(mserd) all
	local ITT_h = e(h_l)

	twoway (scatter mean_outcome1 month_ctr,  mc(gs3) msymbol(Oh)) ///
		(lpoly outcome1 date_init_ctr if date_init_ctr <0, bw(`ITT_h') lw(med) lc(black) k(tri) degree(1)) ///
		(lpoly outcome1 date_init_ctr if date_init_ctr >0, bw(`ITT_h') lw(med) lc(black) k(tri) degree(1)), ///
		legend(off) graphregion(color(white)) ///
		xtitle(No. of Days Since Guideline Change, size(medlarge) margin(med)) ///
		ytitle(Proportion with a Gap in Care, size(medlarge) margin(med)) ///
		ylabel(0(.1).5, labsize(medlarge) format(%02.1f) angle(0) nogrid) ///
		xlabel(-547 -365 -182 0 182 365 547, labsize(medlarge)) legend(off) xline(0, lc(black) lp(dash) lw(thin))


* Outcome 2 (fig 3)

	* create binned averages
	bys month_ctr: egen mean_outcome2 = mean(outcome2)
	bys month_ctr: replace mean_outcome2 = . if _n>1
	
	* plot=
	rdrobust outcome2 date_init_ctr , p(1) kernel(triangular) bwselect(mserd) all
	local ITT_h = e(h_l)
	
	twoway (scatter mean_outcome2 month_ctr, mc(gs3) msymbol(Oh)) ///
		(lpoly outcome2 date_init_ctr if date_init_ctr <0, bw(`ITT_h')  lc(black) k(tri)  degree(1)) ///
		(lpoly outcome2 date_init_ctr if date_init_ctr >0, bw(`ITT_h')  lc(black) k(tri)  degree(1)), ///
		legend(off) graphregion(color(white)) ///
		xtitle(No. of Days Since Guideline Change, size(medlarge) margin(med)) ///
		ytitle(Proportion Not in Care at 12 Months, size(medlarge) margin(med)) ///
		ylabel(0(.1).5, labsize(medlarge) format(%02.1f) angle(0) nogrid) ///
		xlabel(-547 -365 -182 0 182 365 547, labsize(medlarge)) legend(off) xline(0, lc(black) lw(thin) lp(dash)) ///
		title( "A)", span pos(10) col(black)) t1title(" ", margin(small))


* Outcome 3 (fig 3)

	* create binned averages
	bys month_ctr: egen mean_outcome3 = mean(outcome3)
	bys month_ctr: replace mean_outcome3 = . if _n>1
	
	* plot=
	rdrobust outcome3 date_init_ctr , p(1) kernel(triangular) bwselect(mserd) all
	local ITT_h = e(h_l)
	
	twoway (scatter mean_outcome3 month_ctr, mc(gs3) msymbol(Oh)) ///
		(lpoly outcome3 date_init_ctr if date_init_ctr <0, bw(`ITT_h') lc(black) k(tri) degree(1)) ///
		(lpoly outcome3 date_init_ctr if date_init_ctr >0, bw(`ITT_h') lc(black) k(tri) degree(1)), ///
		legend(off) graphregion(color(white)) ///
		xtitle(No. of Days Since Guideline Change, size(medlarge) margin(med)) ///
		ytitle(Proportion with Long-Term Attrition, size(medlarge) margin(med)) ///
		ylabel(0(.1).5, labsize(medlarge) format(%02.1f) angle(0) nogrid) ///
		xlabel(-547 -365 -182 0 182 365 547, labsize(medlarge)) legend(off) xline(0, lc(black)  lw(thin) lp(dash)) ///
		title( "B)", span pos(10) col(black)) t1title(" ", margin(small))
	

* Outcome 4 (fig 3)
	
	* create binned averages
	bys month_ctr: egen mean_outcome4 = mean(outcome4)
	bys month_ctr: replace mean_outcome4 = . if _n>1
	
	* plot
	rdrobust outcome4 date_init_ctr , p(1) kernel(triangular) bwselect(mserd) all
	local ITT_h = e(h_l)
	
	twoway (scatter mean_outcome4 month_ctr, mc(gs3) msymbol(Oh)) ///
		(lpoly outcome4 date_init_ctr if date_init_ctr <0, bw(`ITT_h') lc(black) k(tri) degree(1)) ///
		(lpoly outcome4 date_init_ctr if date_init_ctr >0, bw(`ITT_h') lc(black) k(tri) degree(1)), ///
		legend(off) graphregion(color(white)) ///
		xtitle(No. of Days Since Guideline Change, size(medlarge) margin(med)) ///
		ytitle(Proportion with No Six-Month Viral Load, size(medlarge) margin(med)) ///
		ylabel(0(.1).5, labsize(medlarge) format(%02.1f) angle(0) nogrid) ///
		xlabel(-547 -365 -182 0 182 365 547, labsize(medlarge)) legend(off) xline(0, lc(black) lw(thin) lp(dash)) ///
		title( "C)", span pos(10) col(black)) t1title(" ", margin(small))

	


