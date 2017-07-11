library(probpool)

### Create a new Probpool object from Ranunculaceae dataset
example_pp = probpool(env_pool = Ranunculaceae$Ran_env,
                          disp_pool = Ranunculaceae$Ran_disp,
                          occurrences = Ranunculaceae$Ran_occ,
                          interaction_matrix = Ranunculaceae$Ran_int, interaction_method = 1)

### Inspect Probpool
class(example_pp)
print(example_pp)
show(example_pp)
summary(example_pp)

### Plot Probpool
# 1. Plot entiree probpool
plot(example_pp, main = "Ranunculaceae dataset")

# 2. Plot individual species
plot(example_pp, focal_species = "Thal_simp", main = "Thalictrum simplex")

# 3. Plot individual spatial subset
plot(example_pp, focal_unit = c(1242,1241))
