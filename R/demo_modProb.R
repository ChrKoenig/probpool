# Function that implements the formula we came up with during the workshop
modify_prob = function(p, b){
  # p is the base probability that comes from the env.pool and/or the disp.pool
  # b is the biotic interaction factor (not actually a probability)
  if(is.na(p) | is.na(b)){return(NA)}
  if(b > 0){return(p + (1-p)*b)}
  if(b <= 0){return(p + (p*b))}
}

# Test sequence of bio.pool values (-1: full competition, 0: neutral, +1: full facilitation)
myseq = seq(-1, 1, length.out = 20)
plot(myseq)

# We wanted to be able to estimate a probabilistic species pool when ONLY biotic interactions are known. In this case, our best idea
# was to rescale the biopool to a range between 0 and 1 and treat the values as probabilities that we add up.
# I found this very ugly for two reasons:
#         1. The results are not anymore interpretable as probabilities because biotic interactions can only reduce the 
#            input probability. That means, if we know that a species occurs on a given site (P_ij = 1), this method will return
#            a probability of 0.5 under neutral biotic conditions. That's kind of illogical. Users will be confused
#         2. Within the same R-function, we use two different methods for producing outputs that are superficially similar but not
#            comparable at all. We could write huge disclaimers in the help pages, but I suppose: Users will be confused.

# However, I think I found a solution to this problem. Instead of rescaling the biotic interactions, we assume a base probability 
# that we can modify using the same method as always. When assuming a base probability of 0.5 (presence/absence completely unknown), 
# this approach even produces the same results as the rescaling:

# Rescaling of the sequence
myseq_scaled = (myseq+1)/2

# Modified sequence with base probability of 0.5
myseq_modified = sapply(myseq, modify_prob, p = 0.5)

myseq_scaled == myseq_modified

# Patrick: I like this! One thought to not forget: In the original formula we forgot to account for positive interactions.
# This is why we do all these modifications here. In case a user only has negative interactions measured the original approach
# can still be used and should therefore be implemented as well (maxybe not as default though). We can the explain that to
# account for positive interactions we modified the original approach as follows...

# But assuming a base probability of 0.5 is still rather arbititrary, and in fact in most cases wrong (i.e. it will yield incorrect 
# estimates of species per site). Even when the user does not provide an env.pool and a disp.pool, he/she has to provide 
# distribution data for the species that serve as a basis for the biotic interactions. This gives us two important pieces of 
# information:
#         1. The total number of species in the dataset (N_total)
#         2. The mean number of species per site (N_mean)
# With that, we can calculate a mathematically justified base probability as N_site/N_total

##############
# AN EXAMPLE #
##############
load("data/Ranunculaceae_occurrences.RData")
load("data/Ranunculaceae_bio_prob.RData")

N_total = dim(occ.rst.stack)[3]
N_mean = mean(raster::values(sum(occ.rst.stack)), na.rm = T)
richness = sum(occ.rst.stack)

# 1. Rescaling approach
probs_rescaled = (bio.rst.stack+1)/2
richness_rescaled = sum(probs_rescaled)
N_mean_rescaled = mean(raster::values(richness_rescaled), na.rm = T)

N_mean 
N_mean_rescaled # Estimates way too high

layout(matrix(1:2, ncol = 2))
plot(richness, col = terrain.colors(30), breaks = 1:30)
plot(richness_rescaled, col = terrain.colors(30), breaks = 1:30) # I mean really, no one would consider this a valid result.

# 2. Modification approach
interaction_values = values(sum(bio.rst.stack))
interaction_values_modified = sapply(interaction_values, FUN = modify_prob, p = N_mean)
richness_modified = richness
values(richness_modified) = interaction_values_modified
N_mean_modified = mean(values(richness_modified), na.rm = T)

N_mean
N_mean_modified # not perfect, but closer

plot(richness, col = terrain.colors(30), breaks = 1:30)
plot(richness_modified, col = terrain.colors(30), breaks = 1:30)

# This is a much better estimate of the occurence probabilities, everything else (env + disp) being equal. Don't you think?
# However, be aware that adopting this would invalidate the formula published in the original paper, since the bio.pool can not be 
# treated as another factor that is simply multiplied with the rest. 

# Finally, let's see if both approaches produce similar patterns, irrespective of the total species number.
plot(richness_rescaled)
plot(richness_modified)

# The argument during the retreat was, that only the relative differences between cells are telling, but not the absolute values. The maps show that the relative differences are very similar, but the new approach returns realistic species number estimates on top of that. I therefore don't see any convincing argument for keeping the rescaling+multiplication approach at all. 
