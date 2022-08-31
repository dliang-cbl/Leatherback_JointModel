# JOURNAL PUBLICATION CITATION: Liang et al., integrating telemetry and point observations to inform management and conservation of migratory marine species. Ecosphere.
# R code and example data demonstrating the data fusion. 
Author Name: Dong Liang
Affiliation: University of Maryland Center for Environmental Science, Chesapeake Biological Laboratory

Address: 146 Williams Street, PO Box 0038, Solomons, Maryland, 20688, USA

Email: dliang@umces.edu

# File list and description 
The code requires installation of the following R packages: raster, ctmcmove,mgcv,abind. For parallel processing, the code requires doParallel. The code was tested under R version 4.1.1 and R-Studio version 1.4.1103 under Windows 10.  Most recent versions of R are recommended for testing these codes.

	Main.R - R vignette for joint specie distribution modeling
	
	loc.csv - Simulated fisheries observation data
	track.csv - Simulated tracking data
  seascape.RData - An R list of raster stacks of seascape variables

	nmiss_2.R - Internal helper code
	s2inla_6.R - Internal helper code
	ctmc2glm_study_3.R – Convert tracks to GLM structure
	cellFromXY_nbr.R – Internal helper code
	ppm_st_utils.R – Internal helper code dynamic point process modeling
	ppm_st_glm_6f.R – Build pseudo-absence data from point patterns
	term_utils.R – Internal helper code
	sim2jam_study_1.R – Internal helper code
	all_group_1.R – Internal helper code
	predict_gam_term_2.R – Internal helper code
	sdm_bc_5d.R – Joint modeling from processed telemetry and fishery data
	predict_sdm_bc_5d.R – Predict intensity from joint model
	extract_1.R – Internal helper code
	intensity_2.R – Internal helper code
