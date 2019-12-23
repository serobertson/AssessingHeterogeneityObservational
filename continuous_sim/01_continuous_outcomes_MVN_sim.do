

clear all

set seed 12345678


local N_runs = 1000

timer on 1


		
qui{
						
forvalues scenario  = 1/ 16 {							

	foreach N_pt in  500 1000 5000 {	
	
		foreach prev in  25 50 {

clear all

import excel using "marg_cond_continuous_new.xlsx" , clear first

local sub1_p1 = 	`=sub1_p1[`scenario']'
local sub1_p0 = 	`=sub1_p0[`scenario']'
local sub0_p0 = 	`=sub0_p0[`scenario']'
local sub0_p1 = 	`=sub0_p1[`scenario']'

local sub1_rd = 	`=sub1_RD[`scenario']'
local sub0_rd = 	`=sub0_RD[`scenario']'		
local dif_rd = `=dif_RD[`scenario']'


tempname saved_results 
postfile `saved_results' double(scenario sim_run N_pt prev sub1_rd sub0_rd dif_rd  ///
								sub1_p1 sub1_p0  sub0_p1 sub0_p0 ///
								dif_RD_es1  /// /*ATE (RD=risk difference)*/
								dif_RD_es2  ///  
								dif_RD_es3  ///   	
								dif_RD_es4  ///
								dif_RD_es5  ///
								dif_RD_es6 ///
								dif_RD_es7 ///
								dif_RD_es8 ///
								dif_RD_es9 ///
								dif_RD_es10 ///
								dif_RD_es11 ///
								dif_RD_es12 ///
								dif_RD_es13 ///
								mu_S1_A1_es1 mu_S1_A0_es1 /// /*potential outcome mean in each subgroup */
								mu_S0_A1_es1 mu_S0_A0_es1 ///
								mu_S1_A1_es2 mu_S1_A0_es2 ///
								mu_S0_A1_es2 mu_S0_A0_es2 ///
								mu_S1_A1_es3 mu_S1_A0_es3 ///
								mu_S0_A1_es3 mu_S0_A0_es3 ///
								mu_S1_A1_es4 mu_S1_A0_es4 ///
								mu_S0_A1_es4 mu_S0_A0_es4 ///
								mu_S1_A1_es5 mu_S1_A0_es5 ///
								mu_S0_A1_es5 mu_S0_A0_es5 ///
								mu_S1_A1_es6 mu_S1_A0_es6 ///
								mu_S0_A1_es6 mu_S0_A0_es6 ///
								mu_S1_A1_es7 mu_S1_A0_es7 ///
								mu_S0_A1_es7 mu_S0_A0_es7 ///
								mu_S1_A1_es8 mu_S1_A0_es8 ///
								mu_S0_A1_es8 mu_S0_A0_es8 ///
								mu_S1_A1_es9 mu_S1_A0_es9 ///
								mu_S0_A1_es9 mu_S0_A0_es9 ///
								mu_S1_A1_es10 mu_S1_A0_es10 ///
								mu_S0_A1_es10 mu_S0_A0_es10 ///
								mu_S1_A1_es11 mu_S1_A0_es11 ///
								mu_S0_A1_es11 mu_S0_A0_es11 ///
								mu_S1_A1_es12 mu_S1_A0_es12 ///
								mu_S0_A1_es12 mu_S0_A0_es12 ///
								mu_S1_A1_es13 mu_S1_A0_es13 ///
								mu_S0_A1_es13 mu_S0_A0_es13) using ///
								cont_Y_MVN_rd_sim`scenario'_`N_pt'_`prev'.dta, replace 

forvalues i = 1/`N_runs' {

di in red "scenario = `scenario', prev = `prev', SS = `N_pt', MD1 = `sub1_rd', MD0 = `sub0_rd', simulation = `i'"

clear

set obs `N_pt'

local p = `prev'/100

generate S = rbinomial(1, `p')

/*
forvalues j = 1/8 {
	generate X`j' = rnormal()
}
*/

*Sample from a standard normal distribution with correlation matrix C
matrix c = ( 1, 0.2, -0.2, 0.2, -0.2, 0.2, -0.2, 0.1 \ ///
0.2, 1, 0.1, -0.1, 0.1, -0.1, 0.1, -0.1\ ///
-0.2, 0.1, 1, -0.1, 0.1, -0.1, 0.1, -0.1 \ ///
0.2, -0.1, -0.1, 1, 0.1, -0.1, 0.1, -0.1 \ ///
-0.2, 0.1, 0.1, 0.1, 1, -0.1, 0.1, -0.1 \ ///
0.2, -0.1, -0.1, -0.1, -0.1, 1, 0.1, -0.1 \ ///
-0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 1, -0.1 \ ///
0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, 1 )

drawnorm X1 X2 X3 X4 X5 X6 X7 X8, corr(c) 
		 
generate A = runiform() < invlogit(0 			+ log(1.5) * X1 + log(3) * X2 + log(2) * X3 + log(3) * X4   /// confounders
												+ log(1) * X5 + log(1) * X6 + log(1) * X7 + log(1) * X8  /// pure outcome predictors
											    ) if S == 1									

												
replace A = runiform() < invlogit(0 			+ log(3) * X1 + log(3) * X2 + log(2) * X3 + log(3) * X4   ///
												+ log(1) * X5 + log(1) * X6 + log(1) * X7 + log(1) * X8  	///
												) if S == 0

generate Y = 	    `sub1_p1'  +    3 * X1   + 1.5 * X2   + 1.5 * X3 ///
						+  X4 +   X5 +  X6 ///
						+ X7 + X8 +  rnormal(0, 1) if S == 1 & A==1
												
replace Y = 	    `sub1_p0'  +    3 * X1   + 1.5 * X2   + 1.5 * X3 ///
						+  X4 +   X5 +  X6 ///
						+ X7 + X8 +  rnormal(0, 1) if S == 1 & A==0
						
replace Y = 	    `sub0_p1'  +    1.5 * X1   + 1.5 * X2   + 1.5 * X3 ///
						+  X4 +   X5 +  X6 ///
						+ X7 + X8 +  rnormal(0, 1) if S == 0 & A==1
																								
replace Y = 	    `sub0_p0'  +    1.5 * X1   + 1.5 * X2   + 1.5 * X3 ///
						+  X4 +   X5 +  X6 ///
						+ X7 + X8 +  rnormal(0, 1) if S == 0 & A==0																						
																								
																																															
generate interX1_S = X1 * S
generate inter_A_S = A * S


* OM (outcome modeling)
regress Y X1 X2 X3 X4 X5 X6 X7 X8 A inter_A_S S interX1_S
	matrix estimatesOM = e(b)
	
	/* counterfactual treatment */
	generate tx = 1 
	generate inter_tx = tx * S

	generate g_OM_A1 = 		estimatesOM[1,1] * X1 ///
								+ estimatesOM[1,2] * X2 ///
								+ estimatesOM[1,3] * X3 ///
								+ estimatesOM[1,4] * X4 ///
								+ estimatesOM[1,5] * X5 ///
								+ estimatesOM[1,6] * X6 ///
								+ estimatesOM[1,7] * X7 ///
								+ estimatesOM[1,8] * X8 ///
								+ estimatesOM[1,9] * tx ///
								+ estimatesOM[1,10] * inter_tx ///
								+ estimatesOM[1,11] * S ///
								+ estimatesOM[1,12] * interX1_S ///
								+ estimatesOM[1,13] 
	generate g_OM_A0 = 		estimatesOM[1,1] * X1 ///
								+ estimatesOM[1,2] * X2 ///
								+ estimatesOM[1,3] * X3 ///
								+ estimatesOM[1,4] * X4 ///
								+ estimatesOM[1,5] * X5 ///
								+ estimatesOM[1,6] * X6 ///
								+ estimatesOM[1,7] * X7 ///
								+ estimatesOM[1,8] * X8 ///
								+ estimatesOM[1,11] * S ///
								+ estimatesOM[1,12] * interX1_S ///
								+ estimatesOM[1,13]  

	qui summ g_OM_A1 if S == 1
		local mu_S1_A1_es1 = r(mean)
	qui summ g_OM_A0 if S == 1
		local mu_S1_A0_es1 = r(mean)
		
	qui summ g_OM_A1 if S == 0
		local mu_S0_A1_es1 = r(mean)
	qui summ g_OM_A0 if S == 0
		local mu_S0_A0_es1 = r(mean)
	
		
	local dif_S1 =  `mu_S1_A1_es1' -  `mu_S1_A0_es1' 
	local dif_S0 =  `mu_S0_A1_es1' -  `mu_S0_A0_es1' 
	
	local dif_RD_es1 =  `dif_S1'  - `dif_S0' 
	di `dif_RD_es1'


*******************************************************************	
*******************************************************************
*******************************************************************

*propensity scores with pure outcome predictors (X5-X8) included
	
*IPW1, unnormalized

logit A X1 X2 X3 X4 X5 X6 X7 X8 S interX1_S
	predict ps1, pr
	generate w1 = A * (1/ps1) + (1 - A) * (1/(1-ps1))
	

generate summand_unnorm_S1 = S 
	summ summand_unnorm_S1
	local unnorm_S1 = r(sum)
	
generate summand_unnorm_S0 = (1-S)
	summ summand_unnorm_S0
	local unnorm_S0 = r(sum)

	
generate summand_IPW_S1_A1_term1 = w1 * S * A * Y 

generate summand_IPW_S1_A0_term1 = w1 * S * (1-A) * Y

generate summand_IPW_S0_A1_term1 = w1 * (1-S) * A * Y 

generate summand_IPW_S0_A0_term1 = w1 * (1-S) * (1 - A) * Y 

summ summand_IPW_S1_A1_term1
	local mu_S1_A1_term1 = r(sum)
	
summ summand_IPW_S1_A0_term1
	local mu_S1_A0_term1 = r(sum)
	
summ summand_IPW_S0_A1_term1
	local mu_S0_A1_term1 = r(sum)
	
summ summand_IPW_S0_A0_term1
	local mu_S0_A0_term1 = r(sum)
	
local mu_S1_A1_es2= `mu_S1_A1_term1'/`unnorm_S1' 
local mu_S1_A0_es2= `mu_S1_A0_term1'/`unnorm_S1' 
local mu_S0_A1_es2= `mu_S0_A1_term1'/`unnorm_S0' 
local mu_S0_A0_es2= `mu_S0_A0_term1'/`unnorm_S0'

local dif_S1= `mu_S1_A1_es2'- `mu_S1_A0_es2'
local dif_S0= `mu_S0_A1_es2'- `mu_S0_A0_es2'
local dif_RD_es2 = `dif_S1' - `dif_S0'


display `dif_RD_es2'
									
	
* IPW2
	
*norm term

generate summand_norm_S1_A1 = S * A * w1
	summ summand_norm_S1_A1
	local norm_S1_A1 = r(sum)
	
generate summand_norm_S1_A0 = S * (1-A) * w1
	summ summand_norm_S1_A0
	local norm_S1_A0 = r(sum)
	
generate summand_norm_S0_A1 = (1-S) * A * w1
	summ summand_norm_S0_A1
	local norm_S0_A1 = r(sum)
	
generate summand_norm_S0_A0 = (1-S) * (1-A) * w1
	summ summand_norm_S0_A0
	local norm_S0_A0 = r(sum)
	
local mu_S1_A1_es3= `mu_S1_A1_term1'/`norm_S1_A1' 
local mu_S1_A0_es3= `mu_S1_A0_term1'/`norm_S1_A0' 
local mu_S0_A1_es3= `mu_S0_A1_term1'/`norm_S0_A1' 
local mu_S0_A0_es3= `mu_S0_A0_term1'/`norm_S0_A0'

local dif_S1 = 	`mu_S1_A1_es3' -  `mu_S1_A0_es3'
local dif_S0 = 	`mu_S0_A1_es3' -  `mu_S0_A0_es3'
local dif_RD_es3 = `dif_S1'  - `dif_S0'


display `dif_RD_es3'
	
*DR1, unnormalized
	
generate summand_AIPW_S1_A1_term1 = w1 * S * A * (Y - g_OM_A1) 
generate summand_AIPW_S1_A1_term2= S * g_OM_A1

generate summand_AIPW_S1_A0_term1 = w1 * S * (1-A) * (Y - g_OM_A0) 
generate summand_AIPW_S1_A0_term2= S * g_OM_A0

generate summand_AIPW_S0_A1_term1 = w1 * (1-S) * A * (Y - g_OM_A1) 
generate summand_AIPW_S0_A1_term2= (1-S) * g_OM_A1

generate summand_AIPW_S0_A0_term1 = w1 * (1-S) * (1 - A) * (Y - g_OM_A0) 
generate summand_AIPW_S0_A0_term2= (1-S) * g_OM_A0

summ summand_AIPW_S1_A1_term1
	local mu_S1_A1_term1 = r(sum)
	
summ summand_AIPW_S1_A0_term1
	local mu_S1_A0_term1 = r(sum)
	
summ summand_AIPW_S0_A1_term1
	local mu_S0_A1_term1 = r(sum)
	
summ summand_AIPW_S0_A0_term1
	local mu_S0_A0_term1 = r(sum)
	
summ summand_AIPW_S1_A1_term2
	local mu_S1_A1_term2 = r(sum)
	
summ summand_AIPW_S1_A0_term2
	local mu_S1_A0_term2 = r(sum)
	
summ summand_AIPW_S0_A1_term2
	local mu_S0_A1_term2 = r(sum)
	
summ summand_AIPW_S0_A0_term2
	local mu_S0_A0_term2 = r(sum)
	
generate summand_S1_term2 = S
	summ summand_S1_term2
	local summand_S1_term2 = r(sum)

generate summand_S0_term2 =(1- S)
	summ summand_S0_term2
	local summand_S0_term2 = r(sum)
	
	
local mu_S1_A1_es4= `mu_S1_A1_term1'/`unnorm_S1' + `mu_S1_A1_term2'/`summand_S1_term2'
local mu_S1_A0_es4= `mu_S1_A0_term1'/`unnorm_S1' + `mu_S1_A0_term2'/`summand_S1_term2' 

local mu_S0_A1_es4= `mu_S0_A1_term1'/`unnorm_S0' + `mu_S0_A1_term2'/`summand_S0_term2' 
local mu_S0_A0_es4= `mu_S0_A0_term1'/`unnorm_S0' + `mu_S0_A0_term2'/`summand_S0_term2'

local dif_S1= 	`mu_S1_A1_es4'- `mu_S1_A0_es4'
local dif_S0= 	`mu_S0_A1_es4'- `mu_S0_A0_es4'
local dif_RD_es4 = `dif_S1' - `dif_S0'

/* DR2 */
	
local mu_S1_A1_es5= `mu_S1_A1_term1'/`norm_S1_A1' + `mu_S1_A1_term2'/`summand_S1_term2'
local mu_S1_A0_es5= `mu_S1_A0_term1'/`norm_S1_A0' + `mu_S1_A0_term2'/`summand_S1_term2' 

local mu_S0_A1_es5= `mu_S0_A1_term1'/`norm_S0_A1' + `mu_S0_A1_term2'/`summand_S0_term2' 
local mu_S0_A0_es5= `mu_S0_A0_term1'/`norm_S0_A0' + `mu_S0_A0_term2'/`summand_S0_term2'

local dif_S1= 	`mu_S1_A1_es5'- `mu_S1_A0_es5'
local dif_S0= 	`mu_S0_A1_es5'- `mu_S0_A0_es5'
local dif_RD_es5 = `dif_S1' - `dif_S0'

/* DR3*/

regress Y  X1 X2 X3 X4 X5 X6 X7 X8 A inter_A_S S interX1_S [pw = w1]
	matrix estimates_DR = e(b)
	

	generate g_DR_A1= 		estimates_DR[1,1] * X1 ///
								+ estimates_DR[1,2] * X2 ///
								+ estimates_DR[1,3] * X3 ///
								+ estimates_DR[1,4] * X4 ///
								+ estimates_DR[1,5] * X5 ///
								+ estimates_DR[1,6] * X6 ///
								+ estimates_DR[1,7] * X7 ///
								+ estimates_DR[1,8] * X8 ///
								+ estimates_DR[1,9] * tx ///
								+ estimates_DR[1,10] * inter_tx ///
								+ estimates_DR[1,11] * S ///
								+ estimates_DR[1,12] * interX1_S ///
								+ estimates_DR[1,13] 
	generate g_DR_A0 = 	estimates_DR[1,1] * X1 ///
								+ estimates_DR[1,2] * X2 ///
								+ estimates_DR[1,3] * X3 ///
								+ estimates_DR[1,4] * X4 ///
								+ estimates_DR[1,5] * X5 ///
								+ estimates_DR[1,6] * X6 ///
								+ estimates_DR[1,7] * X7 ///
								+ estimates_DR[1,8] * X8 ///
								+ estimates_DR[1,11] * S ///
								+ estimates_DR[1,12] * interX1_S ///
								+ estimates_DR[1,13]  
	qui summ g_DR_A1 if S == 1
		local mu_S1_A1_es6 = r(mean)
	qui summ g_DR_A0 if S == 1
		local mu_S1_A0_es6 = r(mean)	
		
	qui summ g_DR_A1 if S == 0
		local mu_S0_A1_es6 = r(mean)
	qui summ g_DR_A0 if S == 0
		local mu_S0_A0_es6= r(mean)
		
	local rd_1 =  `mu_S1_A1_es6' -  `mu_S1_A0_es6' 
	local rd_0 =  `mu_S0_A1_es6' -  `mu_S0_A0_es6' 
	
	local dif_RD_es6 =  `rd_1'  - `rd_0' 
	di `dif_RD_es6'

*MT, matching

psmatch2 A if S==1, outcome(Y) pscore(ps1) n(1) logit ate 
local psmatch_S1=r(ate)

*get potential outcome means 
 generate y_1_hat = Y if A==1 & S==1
 generate y_0_hat = Y if A==0 & S==1
 
 sort _id 
 replace y_1_hat = Y[_n1] if A==0 & S==1
 replace y_0_hat = Y[_n1] if A==1 & S==1
 
 summarize  y_1_hat if S==1
 local mu_S1_A1_es7 = r(mean)
  
 summarize y_0_hat if S==1
 local  mu_S1_A0_es7 = r(mean)
 
  
 display `mu_S1_A1_es7' - `mu_S1_A0_es7' /*check matches output of `ate' of psmatch2*/
 
 drop y_1_hat 
drop y_0_hat
drop _*

psmatch2 A if S==0, outcome(Y) pscore(ps1) n(1) logit ate 
local psmatch_S0=r(ate)

*get potential outcome means
 generate y_1_hat = Y if A==1 & S==0
 generate y_0_hat = Y if A==0 & S==0
 
 sort _id 
 replace y_1_hat = Y[_n1] if A==0 & S==0
 replace y_0_hat = Y[_n1] if A==1 & S==0
 
 summarize  y_1_hat if S==0
 local mu_S0_A1_es7 = r(mean)
  
 summarize y_0_hat if S==0
 local  mu_S0_A0_es7 = r(mean)
  
 display `mu_S0_A1_es7' - `mu_S0_A0_es7' 
drop y_1_hat 
drop y_0_hat


local dif_RD_es7=`psmatch_S1'-`psmatch_S0'

display `dif_RD_es7'

*******************************************************************	
*******************************************************************
*******************************************************************

*propensity scores without pure outcome predictors (X5-X8) included

drop summand*
drop g_DR*


*IPW1, unnormalized

logit A X1 X2 X3 X4 S interX1_S
	predict ps2, pr
	generate w2 = A * (1/ps2) + (1 - A) * (1/(1-ps2))

	
generate summand_IPW_S1_A1_term1 = w2 * S * A * Y 

generate summand_IPW_S1_A0_term1 = w2 * S * (1-A) * Y

generate summand_IPW_S0_A1_term1 = w2 * (1-S) * A * Y 

generate summand_IPW_S0_A0_term1 = w2 * (1-S) * (1 - A) * Y 

summ summand_IPW_S1_A1_term1
	local mu_S1_A1_term1 = r(sum)
	
summ summand_IPW_S1_A0_term1
	local mu_S1_A0_term1 = r(sum)
	
summ summand_IPW_S0_A1_term1
	local mu_S0_A1_term1 = r(sum)
	
summ summand_IPW_S0_A0_term1
	local mu_S0_A0_term1 = r(sum)
	
local mu_S1_A1_es8= `mu_S1_A1_term1'/`unnorm_S1' 
local mu_S1_A0_es8= `mu_S1_A0_term1'/`unnorm_S1' 
local mu_S0_A1_es8= `mu_S0_A1_term1'/`unnorm_S0' 
local mu_S0_A0_es8= `mu_S0_A0_term1'/`unnorm_S0'

local dif_S1= 	`mu_S1_A1_es8'- `mu_S1_A0_es8'
local dif_S0= 	`mu_S0_A1_es8'- `mu_S0_A0_es8'
local dif_RD_es8 = `dif_S1' - `dif_S0'


display `dif_RD_es8'
display `mu_S0_A1_es8'
display `mu_S0_A0_es8'
	*brrrr								
	
* IPW2

generate summand_norm_S1_A1 = S * A * w1
	summ summand_norm_S1_A1
	local norm_S1_A1 = r(sum)
	
generate summand_norm_S1_A0 = S * (1-A) * w1
	summ summand_norm_S1_A0
	local norm_S1_A0 = r(sum)
	
generate summand_norm_S0_A1 = (1-S) * A * w1
	summ summand_norm_S0_A1
	local norm_S0_A1 = r(sum)
	
generate summand_norm_S0_A0 = (1-S) * (1-A) * w1
	summ summand_norm_S0_A0
	local norm_S0_A0 = r(sum)
	
local mu_S1_A1_es9= `mu_S1_A1_term1'/`norm_S1_A1' 
local mu_S1_A0_es9= `mu_S1_A0_term1'/`norm_S1_A0' 
local mu_S0_A1_es9= `mu_S0_A1_term1'/`norm_S0_A1' 
local mu_S0_A0_es9= `mu_S0_A0_term1'/`norm_S0_A0'

local dif_S1 = 	`mu_S1_A1_es9' -  `mu_S1_A0_es9'
local dif_S0 = 	`mu_S0_A1_es9' -  `mu_S0_A0_es9'
local dif_RD_es9 = `dif_S1'  - `dif_S0'


display `dif_RD_es9'
	
	
*DR1, unnormalized

generate summand_unnorm_S1 = S 
	summ summand_unnorm_S1
	local unnorm_S1 = r(sum)
	
generate summand_unnorm_S0 = (1-S)
	summ summand_unnorm_S0
	local unnorm_S0 = r(sum)
	
	
generate summand_AIPW_S1_A1_term1 = w2 * S * A * (Y - g_OM_A1) 
generate summand_AIPW_S1_A1_term2= S * g_OM_A1

generate summand_AIPW_S1_A0_term1 = w2 * S * (1-A) * (Y - g_OM_A0) 
generate summand_AIPW_S1_A0_term2= S * g_OM_A0

generate summand_AIPW_S0_A1_term1 = w2 * (1-S) * A * (Y - g_OM_A1) 
generate summand_AIPW_S0_A1_term2= (1-S) * g_OM_A1

generate summand_AIPW_S0_A0_term1 = w2 * (1-S) * (1 - A) * (Y - g_OM_A0) 
generate summand_AIPW_S0_A0_term2= (1-S) * g_OM_A0

summ summand_AIPW_S1_A1_term1
	local mu_S1_A1_term1 = r(sum)
	
summ summand_AIPW_S1_A0_term1
	local mu_S1_A0_term1 = r(sum)
	
summ summand_AIPW_S0_A1_term1
	local mu_S0_A1_term1 = r(sum)
	
summ summand_AIPW_S0_A0_term1
	local mu_S0_A0_term1 = r(sum)
	
summ summand_AIPW_S1_A1_term2
	local mu_S1_A1_term2 = r(sum)
	
summ summand_AIPW_S1_A0_term2
	local mu_S1_A0_term2 = r(sum)
	
summ summand_AIPW_S0_A1_term2
	local mu_S0_A1_term2 = r(sum)
	
summ summand_AIPW_S0_A0_term2
	local mu_S0_A0_term2 = r(sum)
	
generate summand_S1_term2 = S
	summ summand_S1_term2
	local summand_S1_term2 = r(sum)

generate summand_S0_term2 =(1- S)
	summ summand_S0_term2
	local summand_S0_term2 = r(sum)
	
	
local mu_S1_A1_es10= `mu_S1_A1_term1'/`unnorm_S1' + `mu_S1_A1_term2'/`summand_S1_term2'
local mu_S1_A0_es10= `mu_S1_A0_term1'/`unnorm_S1' + `mu_S1_A0_term2'/`summand_S1_term2' 

local mu_S0_A1_es10= `mu_S0_A1_term1'/`unnorm_S0' + `mu_S0_A1_term2'/`summand_S0_term2' 
local mu_S0_A0_es10= `mu_S0_A0_term1'/`unnorm_S0' + `mu_S0_A0_term2'/`summand_S0_term2'

local dif_S1= 	`mu_S1_A1_es10'- `mu_S1_A0_es10'
local dif_S0= 	`mu_S0_A1_es10'- `mu_S0_A0_es10'
local dif_RD_es10 = `dif_S1' - `dif_S0'

/* DR2, normalized */
	
local mu_S1_A1_es11= `mu_S1_A1_term1'/`norm_S1_A1' + `mu_S1_A1_term2'/`summand_S1_term2'
local mu_S1_A0_es11= `mu_S1_A0_term1'/`norm_S1_A0' + `mu_S1_A0_term2'/`summand_S1_term2' 

local mu_S0_A1_es11= `mu_S0_A1_term1'/`norm_S0_A1' + `mu_S0_A1_term2'/`summand_S0_term2' 
local mu_S0_A0_es11= `mu_S0_A0_term1'/`norm_S0_A0' + `mu_S0_A0_term2'/`summand_S0_term2'

local dif_S1= 	`mu_S1_A1_es11'- `mu_S1_A0_es11'
local dif_S0= 	`mu_S0_A1_es11'- `mu_S0_A0_es11'
local dif_RD_es11 = `dif_S1' - `dif_S0'

/* DR3 */

regress Y  X1 X2 X3 X4 A inter_A_S S interX1_S [pw = w2]
	matrix estimates_M6 = e(b)
	

	generate g_DR_A1= 		estimates_DR[1,1] * X1 ///
								+ estimates_DR[1,2] * X2 ///
								+ estimates_DR[1,3] * X3 ///
								+ estimates_DR[1,4] * X4 ///
								+ estimates_DR[1,5] * X5 ///
								+ estimates_DR[1,6] * X6 ///
								+ estimates_DR[1,7] * X7 ///
								+ estimates_DR[1,8] * X8 ///
								+ estimates_DR[1,9] * tx ///
								+ estimates_DR[1,10] * inter_tx ///
								+ estimates_DR[1,11] * S ///
								+ estimates_DR[1,12] * interX1_S ///
								+ estimates_DR[1,13] 
	generate g_DR_A0 = 	estimates_DR[1,1] * X1 ///
								+ estimates_DR[1,2] * X2 ///
								+ estimates_DR[1,3] * X3 ///
								+ estimates_DR[1,4] * X4 ///
								+ estimates_DR[1,5] * X5 ///
								+ estimates_DR[1,6] * X6 ///
								+ estimates_DR[1,7] * X7 ///
								+ estimates_DR[1,8] * X8 ///
								+ estimates_DR[1,11] * S ///
								+ estimates_DR[1,12] * interX1_S ///
								+ estimates_DR[1,13]  
	qui summ g_DR_A1 if S == 1
		local mu_S1_A1_es12 = r(mean)
	qui summ g_DR_A0 if S == 1
		local mu_S1_A0_es12 = r(mean)	
		
	qui summ g_DR_A1 if S == 0
		local mu_S0_A1_es12 = r(mean)
	qui summ g_DR_A0 if S == 0
		local mu_S0_A0_es12= r(mean)
		
	local rd_1 =  `mu_S1_A1_es12' -  `mu_S1_A0_es12' 
	local rd_0 =  `mu_S0_A1_es12' -  `mu_S0_A0_es12' 
	
	local dif_RD_es12 =  `rd_1'  - `rd_0' 
	di `dif_RD_es12'
*MT, matching


psmatch2 A if S==1, outcome(Y) pscore(ps2) n(1) logit ate 
local psmatch_S1=r(ate)

*get potential outcome means 
 generate y_1_hat = Y if A==1 & S==1
 generate y_0_hat = Y if A==0 & S==1
 
 sort _id 
 replace y_1_hat = Y[_n1] if A==0 & S==1
 replace y_0_hat = Y[_n1] if A==1 & S==1
 
 summarize  y_1_hat if S==1
 local mu_S1_A1_es13 = r(mean)
  
 summarize y_0_hat if S==1
 local  mu_S1_A0_es13 = r(mean)
 
  
 display `mu_S1_A1_es13' - `mu_S1_A0_es13' /*check matches output of `ate' of psmatch2*/
 
 drop y_1_hat 
drop y_0_hat
drop _*

psmatch2 A if S==0, outcome(Y) pscore(ps2) n(1) logit ate 
local psmatch_S0=r(ate)

*get potential outcome means
 generate y_1_hat = Y if A==1 & S==0
 generate y_0_hat = Y if A==0 & S==0
 
 sort _id 
 replace y_1_hat = Y[_n1] if A==0 & S==0
 replace y_0_hat = Y[_n1] if A==1 & S==0
 
 summarize  y_1_hat if S==0
 local mu_S0_A1_es13 = r(mean)
  
 summarize y_0_hat if S==0
 local  mu_S0_A0_es13 = r(mean)
  
 display `mu_S0_A1_es13' - `mu_S0_A0_es13' 
drop y_1_hat 
drop y_0_hat


local dif_RD_es13=`psmatch_S1'-`psmatch_S0'

display `dif_RD_es13'

	
post `saved_results'  	(`scenario') (`i') (`N_pt') (`prev') (`sub1_rd') ///
						(`sub0_rd') (`dif_rd') ///
						(`sub1_p1') (`sub1_p0') (`sub0_p1') (`sub0_p0') ///
						(`dif_RD_es1')  ///
						(`dif_RD_es2')  ///
						(`dif_RD_es3')  ///
						(`dif_RD_es4')  ///
						(`dif_RD_es5')  ///
						(`dif_RD_es6') ///
						(`dif_RD_es7') ///
						(`dif_RD_es8') ///
						(`dif_RD_es9') ///
						(`dif_RD_es10') ///
						(`dif_RD_es11') ///
						(`dif_RD_es12') ///
						(`dif_RD_es13') ///
						(`mu_S1_A1_es1') (`mu_S1_A0_es1') ///
						(`mu_S0_A1_es1') (`mu_S0_A0_es1') ///
						(`mu_S1_A1_es2') (`mu_S1_A0_es2') ///
						(`mu_S0_A1_es2') (`mu_S0_A0_es2') ///
						(`mu_S1_A1_es3') (`mu_S1_A0_es3') ///
						(`mu_S0_A1_es3') (`mu_S0_A0_es3') ///
						(`mu_S1_A1_es4') (`mu_S1_A0_es4') ///
						(`mu_S0_A1_es4') (`mu_S0_A0_es4') ///
						(`mu_S1_A1_es5') (`mu_S1_A0_es5') ///
						(`mu_S0_A1_es5') (`mu_S0_A0_es5') ///
						(`mu_S1_A1_es6') (`mu_S1_A0_es6') ///
						(`mu_S0_A1_es6') (`mu_S0_A0_es6') ///
						(`mu_S1_A1_es7') (`mu_S1_A0_es7') ///
						(`mu_S0_A1_es7') (`mu_S0_A0_es7') ///
						(`mu_S1_A1_es8') (`mu_S1_A0_es8') ///
						(`mu_S0_A1_es8') (`mu_S0_A0_es8') ///
						(`mu_S1_A1_es9') (`mu_S1_A0_es9') ///
						(`mu_S0_A1_es9') (`mu_S0_A0_es9') ///
						(`mu_S1_A1_es10') (`mu_S1_A0_es10') ///
						(`mu_S0_A1_es10') (`mu_S0_A0_es10') ///
						(`mu_S1_A1_es11') (`mu_S1_A0_es11') ///
						(`mu_S0_A1_es11') (`mu_S0_A0_es11') ///
						(`mu_S1_A1_es12') (`mu_S1_A0_es12') ///
						(`mu_S0_A1_es12') (`mu_S0_A0_es12') ///
						(`mu_S1_A1_es13') (`mu_S1_A0_es13') ///
						(`mu_S0_A1_es13') (`mu_S0_A0_es13') 
								
						
				}
postclose `saved_results'
//
			}
//
		}
//
	}
//

}
//






 timer off 1
 
 
 timer list 1

 
 
 
