clear all


*****************************************************************************************
forvalues scenario = 1/16 {

	foreach N_pt in  500  1000  5000  {
		foreach prev in 25 50 {
		*append using "cont_Y_IID_rd_sim`scenario'_`N_pt'_`prev'.dta"
		append using "cont_Y_MVN_rd_sim`scenario'_`N_pt'_`prev'.dta"

}
//
}
//
}
//

forvalues i = 1/13 {

	*bias of diff_rd
	generate bias_diff_rd_es`i' = (dif_RD_es`i' - dif_rd)
	generate sq_err_diff_rd_es`i' = (dif_RD_es`i' - dif_rd)^2
	
	*bias of each potential outcome mean
	generate bias_mu_S1_A1_es`i' = (mu_S1_A1_es`i' - sub1_p1)
	generate sq_err_mu_S1_A1_es`i' = (mu_S1_A1_es`i' - sub1_p1)^2
	
	generate bias_mu_S1_A0_es`i' = (mu_S1_A0_es`i' - sub1_p0)
	generate sq_err_mu_S1_A0_es`i' = (mu_S1_A0_es`i' - sub1_p0)^2
	
	generate bias_mu_S0_A1_es`i' = (mu_S0_A1_es`i' - sub0_p1)
	generate sq_err_mu_S0_A1_es`i' = (mu_S0_A1_es`i' - sub0_p1)^2
	
	generate bias_mu_S0_A0_es`i' = (mu_S0_A0_es`i' - sub0_p0)
	generate sq_err_mu_S0_A0_es`i' = (mu_S0_A0_es`i' - sub0_p0)^2
}
//

forvalues i = 1/13 {
	foreach estimator_name in diff_rd mu_S1_A1 mu_S1_A0 mu_S0_A1 mu_S0_A0 {

	bysort scenario N_pt prev: egen mean_bias_`estimator_name'_es`i' = mean(bias_`estimator_name'_es`i')
	bysort scenario N_pt prev: egen mean_sq_err_`estimator_name'_es`i' = mean(sq_err_`estimator_name'_es`i')

}
}
//


duplicates drop scenario N_pt prev , force



forvalues i = 1/13 {
foreach estimator_name in diff_rd mu_S1_A1 mu_S1_A0 mu_S0_A1 mu_S0_A0 {

	generate mean_variance_`estimator_name'_es`i'  = mean_sq_err_`estimator_name'_es`i'  - (mean_bias_`estimator_name'_es`i' )^2

}
}
//

	
save "summaries_continuous.dta" ,replace

*****************************************************************************************

use "summaries_continuous.dta" , clear




generate prev_str = string(prev/100, "%9.2f")
generate prev_str2 = string(prev, "%9.0f")

generate N_pt_str = string(N_pt  /100, "%9.2f")

generate sub1_rd_str = string(sub1_rd, "%9.3f")

generate sub0_rd_str = string(sub0_rd, "%9.3f")

generate diff_rd_str = string(dif_rd, "%9.1f")


sort prev N_pt scenario



foreach measure in 	mean_bias ///
					mean_sq_err  ///
					mean_variance {
	local form = "%9.3f"	
							
forvalues i = 1/13 {	
foreach estimator_name in diff_rd mu_S1_A1 mu_S1_A0 mu_S0_A1 mu_S0_A0 {
				
	generate str_`measure'_`estimator_name'_es`i' = string(`measure'_`estimator_name'_es`i', "`form'" )
}

	}
}	

**********	
foreach prev in 25 50 { 
foreach measure in 	mean_bias ///
					mean_sq_err  ///
					mean_variance {
foreach estimator_name in diff_rd mu_S1_A1 mu_S1_A0 mu_S0_A1 mu_S0_A0 {
generate OM = str_`measure'_`estimator_name'_es1
generate IPW1 = str_`measure'_`estimator_name'_es2					
generate IPW2 = str_`measure'_`estimator_name'_es3					
generate DR1 = str_`measure'_`estimator_name'_es4					
generate DR2 = str_`measure'_`estimator_name'_es5					
generate DR3 = str_`measure'_`estimator_name'_es6
generate MT = str_`measure'_`estimator_name'_es7

generate IPW1_star = str_`measure'_`estimator_name'_es8				
generate IPW2_star = str_`measure'_`estimator_name'_es9				
generate DR1_star = str_`measure'_`estimator_name'_es10					
generate DR2_star = str_`measure'_`estimator_name'_es11					
generate DR3_star = str_`measure'_`estimator_name'_es12
generate MT_star = str_`measure'_`estimator_name'_es13						
					

export delimited scenario N_pt prev_str2 sub1_rd_str sub0_rd_str diff_rd_str ///
  OM IPW1 IPW2 DR1 DR2 DR3 MT ///
  IPW1_star IPW2_star DR1_star DR2_star DR3_star MT_star ///
  using formatted_`measure'_prev`prev'_`estimator_name'.csv if prev_str2 == "`prev'", replace 
  
  drop OM 
  drop IPW1
  drop IPW2
  drop DR1
  drop DR2
  drop DR3
  drop MT
  drop IPW1_star
  drop IPW2_star
  drop DR1_star
  drop DR2_star
  drop DR3_star
  drop MT_star
	
	}	
}	

}



