# ####-------------------------- Test Code -----------------------------####
# rm(list=ls())
# 
# load("data/Ranunculaceae_occurrences.RData")
# load("data/Ranunculaceae_disp_prob.RData")
# load("data/Ranunculaceae_env_prob.RData")
# load("data/Ranunculaceae_interactions.RData")
# 
# # create interaction matrix
# require(cooccur)
# require(plyr)
# spec_site = adply(names(occ.rst.stack), 1, .fun = function(name){values(occ.rst.stack[[name]])}, .id = NULL)
# spec_site = spec_site[, apply(spec_site, 2, function(x){!all(is.na(x))})]
# rownames(spec_site) = names(occ.rst.stack)
# interaction.matrix = cooccur(spec_site, thresh = FALSE, spp_names = TRUE, only_effects = TRUE, eff_standard = TRUE,
#                              eff_matrix = TRUE) # takes > 10 minutes
# save(interaction.matrix, file = "data/Ranunculaceae_interactions.RData")
# 
# ##########
# load(file = "data/Ranunculaceae_disp_prob.RData")
# load(file = "data/Ranunculaceae_bio_prob.RData")
# load(file = "data/Ranunculaceae_env_prob.RData")
# load(file = "data/Ranunculaceae_interactions.RData")
# test_probpool = prob.pool(env.pool = suit.rst.stack, 
#                           disp.pool = disp.rst.stack, 
#                           occurrences = occ.rst.stack,
#                           interaction.matrix = interaction.matrix, interaction.method = 1)
# 
# summary(test_probpool)
# plot(test_probpool)
# plot(test_probpool, species = "Thal_simp")
# plot(test_probpool, focal.unit = c(1242,1241))
# 
# test_probpool = prob.pool(env.pool = suit.rst.stack, 
#                           disp.pool = disp.rst.stack, 
#                           occurrences = occ.rst.stack,
#                           interaction.matrix = interaction.matrix, interaction.method = 2)
# 
# summary(test_probpool)
# plot(test_probpool)
# 
# test_probpool = prob.pool(env.pool = suit.rst.stack, 
#                           disp.pool = disp.rst.stack, 
#                           interaction.matrix = interaction.matrix, interaction.method = 1)
# 
# summary(test_probpool)
# plot(test_probpool)
# 
# test_probpool = prob.pool(env.pool = suit.rst.stack, 
#                           disp.pool = disp.rst.stack)
# 
# summary(test_probpool)
# plot(test_probpool)
# 
# test_probpool = prob.pool(occurrences = occ.rst.stack,
#                           interaction.matrix = interaction.matrix, interaction.method = 1)
# 
# summary(test_probpool)
# plot(test_probpool)
# plot(test_probpool, species = "Thal_simp")
# plot(test_probpool, focal.unit = extent(c(7,8,49,53)))
