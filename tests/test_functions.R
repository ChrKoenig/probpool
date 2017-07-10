# rm(list = ls())
# source("R/probpool_functions.R")
# source("R/probpool_class.R")
# 
# ###########################################################
# # Ranunculaceae
# load("data/Ranunculaceae_disp_prob.RData")
# load("data/Ranunculaceae_occurrences.RData")
# load("data/Ranunculaceae_env_prob.RData")
#  
# ###########
# par(mfrow = c(1,3))
# # Only Env
# test1 = prob.pool(env.pool = suit.rst.stack)
# summary(test1)
# plot(sum(test1@pools$env.pool), main = "env", col = rev(terrain.colors(35)))
# plot(sum(test1@pools$disp.pool), main = "disp")
# plot(sum(test1@pools$prob.pool), main = "prob", col = rev(terrain.colors(35)))
# 
# # Only Disp
# test2 = prob.pool(disp.pool = disp.rst.stack)
# summary(test2)
# plot(sum(test2@pools$env.pool), main = "env", col = rev(terrain.colors(35)),  breaks = 1:35)
# plot(sum(test2@pools$disp.pool), main = "disp")
# plot(sum(test2@pools$prob.pool), main = "prob", col = rev(terrain.colors(35)))
# 
# # Env + Disp
# test3 = prob.pool(env.pool = suit.rst.stack, disp.pool = disp.rst.stack)
# summary(test3)
# par(mfrow=c(1,3))
# plot(sum(test3@pools$env.pool), main = "env", col = rev(terrain.colors(35)),  breaks = 1:35)
# plot(sum(test3@pools$disp.pool), main = "disp")
# plot(sum(test3@pools$prob.pool), main = "prob", col = rev(terrain.colors(35)), breaks = 1:35)
# 
# ###########################################################
# # Simulated data
# load("data/Simu_occurrences.RData")
# load("data/Simu_env_prob.RData")
# load("data/Simu_InteractionMatrix.RData") # Interactions are very small and even out
# 
# ###########
# par(mfrow = c(1,3))
# # Only Env
# test4 = prob.pool(env.pool = Simu_env_prob)
# summary(test4)
# plot(sum(test4@pools$env.pool), main = "env")
# plot(sum(test4@pools$disp.pool), main = "disp")
# plot(sum(test4@pools$prob.pool), main = "prob")
# 
# # Env + Interaction (modification)
# test5 = prob.pool(env.pool = Simu_env_prob, 
#                   interaction.matrix = Simu_InteractionMatrix)
# summary(test5)
# par(mfrow=c(1,3))
# plot(sum(test5@pools$env.pool), main = "env")
# plot(sum(test5@pools$disp.pool), main = "disp")
# plot(sum(test5@pools$prob.pool), main = "prob")
# 
# test5 = prob.pool(env.pool = Simu_env_prob, 
#                   interaction.matrix = Simu_InteractionMatrix)
# 
# # Env + Interaction (multiplication)
# test6 = prob.pool(env.pool = Simu_env_prob, 
#                   interaction.matrix = Simu_InteractionMatrix,
#                   interaction.method = 2)
# summary(test6)
# par(mfrow=c(1,3))
# plot(sum(test6@pools$env.pool), main = "env")
# plot(sum(test6@pools$disp.pool), main = "disp")
# plot(sum(test6@pools$prob.pool), main = "prob")
# 
# # Interaction just from occurrence (modification)
# test7 = prob.pool(occurrences = simu_occ_rst_stack, 
#                   interaction.matrix = Simu_InteractionMatrix,
#                   interaction.method = 1)
# summary(test7)
# par(mfrow=c(1,2))
# plot(sum(test7@pools$occurrences), main = "occurrence")
# plot(sum(test7@pools$prob.pool), main = "prob") # Interaction matrix is nearly neutral, so no effect on prob.pool
# 
# 
# # Interaction just from occurrence (multiplication)
# test8 = prob.pool(occurrences = simu_occ_rst_stack, 
#                   interaction.matrix = Simu_InteractionMatrix,
#                   interaction.method = 2)
# summary(test8)
# par(mfrow=c(1,2))
# plot(sum(test8@pools$occurrences), main = "occurrence")
# plot(sum(test8@pools$prob.pool), main = "prob") # Probabilities (--> species numbers) halved due to rescaling + multiplication
# 
