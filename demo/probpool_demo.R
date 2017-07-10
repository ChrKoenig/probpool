#' @include probpool_class.R
#' @include probpool_function.R
#' @include probpool_utils
NULL

# # create interaction matrix
# load("data/Ranunculaceae_occurrences.RData")
# spec_site = adply(names(occ.rst.stack), 1, .fun = function(name){values(occ.rst.stack[[name]])}, .id = NULL)
# spec_site = spec_site[, apply(spec_site, 2, function(x){!all(is.na(x))})]
# rownames(spec_site) = names(occ.rst.stack)
# interaction_matrix = cooccur(spec_site, thresh = FALSE, spp_names = TRUE, only_effects = TRUE, eff_standard = TRUE,
#                              eff_matrix = TRUE) # takes > 10 minutes
# save(interaction_matrix, file = "data/Ranunculaceae_interaction_matrix.RData")
# 
##########
# load(file = "data/Ranunculaceae_bio_prob.RData")
# load(file = "data/Ranunculaceae_env_prob.RData")
# load(file = "data/Ranunculaceae_occurrences.RData")
# load(file = "data/Ranunculaceae_interaction_matrix.RData")
test_probpool = probpool(env_pool = example1_env,
                          disp_pool = example1_disp,
                          occurrences = example1_occ,
                          interaction_matrix = example1_int, interaction_method = 1)
class(test_probpool)
print(test_probpool)
summary(test_probpool)
plot(test_probpool)
plot(test_probpool, species = "Thal_simp")
plot(test_probpool, focal.unit = c(1242,1241))
# 
# test_probpool = probpool(env_pool = suit.rst.stack,
#                           disp_pool = disp.rst.stack,
#                           occurrences = occ.rst.stack,
#                           interaction_matrix = interaction_matrix, interaction_method = 2)
# 
# summary(test_probpool)
# plot(test_probpool)
# 
# test_probpool = probpool(env_pool = suit.rst.stack,
#                           disp_pool = disp.rst.stack)
# 
# summary(test_probpool)
# plot(test_probpool)
# 
# test_probpool = probpool(env_pool = suit.rst.stack,
#                           disp_pool = disp.rst.stack)
# 
# summary(test_probpool)
# plot(test_probpool)
# 
# test_probpool = probpool(occurrences = occ.rst.stack,
#                           interaction_matrix = interaction_matrix, interaction_method = 1)
# 
# summary(test_probpool)
# plot(test_probpool)
# plot(test_probpool, species = "Thal_simp")
# plot(test_probpool, focal.unit = extent(c(7,8,49,53)))
