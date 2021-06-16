#' """ Example code for MixSIAR using vignettes for introduction
#'     with bonefish example
#'     @author: Ryan James
#'     date: 6/16/21"""

# *** show rainbow parantheses, ctrl + shift + c, and soft wrap 

# set working directory

# load libraries 
library(MixSIAR)
library(tidyverse)
options(max.print = 6000000)

# this script combines multiple vignettes from the MixSIAR github page (https://brianstock.github.io/MixSIAR/) to show some of the common uses for MixSIAR for mixing models. The last example is with real data to show how you can apply real data         
# MixSIAR manual (https://github.com/brianstock/MixSIAR/blob/master/inst/mixsiar_manual_small.pdf)

# Example 1 wolves (simple) https://brianstock.github.io/MixSIAR/articles/wolves_ex.html -----

# load the mix data
mix.filename = system.file("extdata", "wolves_consumer.csv", package = "MixSIAR")

# look at the dataset 
wvs = read_csv(mix.filename)
wvs


# Load the mixture/consumer data
mix = load_mix_data(filename=mix.filename, # file of consumers iso
                    iso_names=c("d13C","d15N"), # iso column names
                    factors=c("Region"), # factor column names
                    fac_random=c(F), 
                    fac_nested=c(F), 
                    cont_effects=NULL)

# load the source data
source.filename = system.file("extdata", "wolves_sources.csv", package = "MixSIAR")

# look at source data 
swv = read_csv(source.filename)
swv

# Load the source data
source = load_source_data(filename=source.filename,
                          source_factors="Region", 
                          conc_dep=FALSE, 
                          data_type="means", 
                          mix)

# load trophic discrimination
discr.filename = system.file("extdata", "wolves_discrimination.csv", package = "MixSIAR")

# look at tdf 
tdw = read_csv(discr.filename)
tdw

# Load the discrimination/TDF data
discr = load_discr_data(filename=discr.filename, mix)

# plot biplot
# Make an isospace plot
plot_data(filename="isospace_plot", plot_save_pdf=F, plot_save_png=F, mix,source,discr) 
# check salmon 

# Write the JAGS model file
model_filename = "MixSIAR_model.txt"   # Name of the JAGS model file
resid_err = TRUE #variation in consumer assimilation (STock and Semmens 2016 for in depth explanation https://esajournals-onlinelibrary-wiley-com.ezproxy.fiu.edu/doi/10.1002/ecy.1517)
process_err = TRUE #variation in sampling of source values
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

# run model 
# | run ==  | Chain Length | Burn-in | Thin | # Chains |
#   | ------------- | ------------- | ------------- | ------------- | ------------- |
#   | "test" | 1,000 | 500 | 1 | 3 |
#   | "very short" | 10,000 | 5,000 | 5 | 3 |
#   | "short" | 50,000 | 25,000 | 25 | 3 |
#   | "normal" | 100,000 | 50,000 | 50 | 3 |
#   | "long" | 300,000 | 200,000 | 100 | 3 |
#   | "very long" | 1,000,000 | 500,000 | 500 | 3 |
#   | "extreme" | 3,000,000 | 1,500,000 | 500 | 3 |

jags.w = run_model(run="test", mix, source, discr, model_filename)


# set output options
output_options = list(summary_save = TRUE,
                      summary_name = "mm_results/wolves_summary_statistics",
                      sup_post = FALSE,
                      plot_post_save_pdf = F,
                      plot_post_name = "posterior_density",
                      sup_pairs = FALSE,
                      plot_pairs_save_pdf = F,
                      plot_pairs_name = "pairs_plot",
                      sup_xy = TRUE,
                      plot_xy_save_pdf = FALSE,
                      plot_xy_name = "xy_plot",
                      gelman = TRUE,
                      heidel = FALSE,
                      geweke = TRUE,
                      diag_save = T,
                      diag_name = "mm_results/wolves_diagnostics",
                      indiv_effect = FALSE,
                      plot_post_save_png = FALSE,
                      plot_pairs_save_png = FALSE,
                      plot_xy_save_png = FALSE,
                      diag_save_ggmcmc = FALSE)
# display outputs
output_JAGS(jags.w, mix, source, output_options)

# Example 2 Geese (Concentration dependence)----
# load consumer data
mix.filename = system.file("extdata", "geese_consumer.csv", package = "MixSIAR")

# look at data
gc = read_csv(mix.filename)
gc

mix = load_mix_data(filename=mix.filename,
                     iso_names=c("d13C","d15N"),
                     factors="Group",
                     fac_random=FALSE,
                     fac_nested=FALSE,
                     cont_effects=NULL)

# load source data 
source.filename = system.file("extdata", "geese_sources.csv", package = "MixSIAR")

# look at source data
gs = read_csv(source.filename)
gs

source = load_source_data(filename=source.filename,
                           source_factors=NULL,
                           conc_dep=TRUE,
                           data_type="means",
                           mix)

# Load tdf
discr.filename = system.file("extdata", "geese_discrimination.csv", package = "MixSIAR")

# look at tef data
gt = read_csv(discr.filename)
gt 

discr = load_discr_data(filename=discr.filename, mix)

# plot biplot
# Make an isospace plot
plot_data(filename="isospace_plot", plot_save_pdf=F, plot_save_png=F, mix,source,discr) 
# check salmon 

# Write the JAGS model file
model_filename = "MixSIAR_model.txt"   # Name of the JAGS model file
resid_err = TRUE #variation in consumer assimilation (STock and Semmens 2016 for in depth explanation https://esajournals-onlinelibrary-wiley-com.ezproxy.fiu.edu/doi/10.1002/ecy.1517)
process_err = TRUE #variation in sampling of source values
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

jags.g = run_model(run="test", mix, source, discr, model_filename)


# set output options
output_options = list(summary_save = TRUE,
                      summary_name = "mm_results/geese_ss",
                      sup_post = FALSE,
                      plot_post_save_pdf = F,
                      plot_post_name = "posterior_density",
                      sup_pairs = FALSE,
                      plot_pairs_save_pdf = F,
                      plot_pairs_name = "pairs_plot",
                      sup_xy = TRUE,
                      plot_xy_save_pdf = FALSE,
                      plot_xy_name = "xy_plot",
                      gelman = TRUE,
                      heidel = FALSE,
                      geweke = TRUE,
                      diag_save = T,
                      diag_name = "mm_results/geese_diag",
                      indiv_effect = FALSE,
                      plot_post_save_png = FALSE,
                      plot_pairs_save_png = FALSE,
                      plot_xy_save_png = FALSE,
                      diag_save_ggmcmc = FALSE)
# display outputs
output_JAGS(jags.g, mix, source, output_options)

# Example 3 bonefish (S, Nested, Concentration dependence)----
# load consumer data
mix = load_mix_data(file("data/bf.csv"),
                    iso_names=c("d13C","d15N","d34S"),
                    factors= c("Site","ID"),
                    fac_random=c(F, T),
                    fac_nested=c(F, T),
                    cont_effects=NULL)

# load source data
source = load_source_data(file("data/sFLbayMAN.csv"),
                          source_factors=NULL,
                          conc_dep=TRUE,
                          data_type="means",
                          mix)
# load TEF data
discr = load_discr_data(file("data/FLbayTEF3.5.csv"), mix)

# Make an isospace plot
plot_data(filename="isospace_plot", plot_save_pdf=FALSE, plot_save_png=FALSE, mix,source,discr)

# Write the JAGS model file
model_filename = "MixSIAR_model.txt"
resid_err = T
process_err = T
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

# run jags
jags.bf = run_model(run="test", mix, source, discr, model_filename,
                    alpha.prior = 1, resid_err, process_err)

# Process JAGS output
output_bf = list(summary_save = TRUE,
                 summary_name = "mm_results/bf_ss",
                 sup_post = FALSE,
                 plot_post_save_pdf = FALSE,
                 plot_post_name = "lower_posterior_density",
                 sup_pairs = FALSE,
                 plot_pairs_save_pdf = FALSE,
                 plot_pairs_name = "lower_pairs_plot",
                 sup_xy = TRUE,
                 plot_xy_save_pdf = FALSE,
                 plot_xy_name = "lower_xy_plot",
                 gelman = TRUE,
                 heidel = FALSE,
                 geweke = TRUE,
                 diag_save = TRUE,
                 diag_name = "mm_results/bf_diag",
                 indiv_effect = FALSE,
                 plot_post_save_png = F,
                 plot_pairs_save_png = FALSE,
                 plot_xy_save_png = FALSE)

output_JAGS(jags.bf, mix, source, output_bf)
