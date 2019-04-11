## BEDASSLE analysis code from Proceedings of the Royal Society B. DOI: 10.1098/rspb.2019.0431
## Run BEDASSLE for 10 million generations using the beta binomial model with geographic distance and genetic distance normalized
## 
## Date Nov 2017

library(BEDASSLE)


###################### RUN AN MCMC model #########################

# MCMC(counts, sample_sizes, D, E, k, loci, delta, aD_stp, aE_stp, a2_stp, thetas_stp, mu_stp, ngen, printfreq, savefreq, samplefreq, directory = NULL, prefix = "", continue = FALSE, continuing.params = NULL)

#upload metadata- sample and population
m2 <- read.csv("hybrid_score.csv", header = TRUE)

#upload an OTU table, make it presence / absence, and collapse it by population such the rownames are population and colnames are OTUs
otu1 <- read.table("otu_table_no_singletons_CSS_transformed.txt", header = TRUE, sep = '\t', row.names = 1)
#make a key of otu ids and taxonomy
otu_taxa_key <- subset(otu1, select = c("taxonomy"))
#subset to samples in metadata
otu1 <- otu1[,colnames(otu1) %in% m2$sample]
#drop otus that appear in 0 or 1 samples
otu1$sumz <- rowSums(otu1 != 0)
otu1 <- subset(otu1, sumz > 1)
otu1$sumz <- NULL
#make it p/a
otu1[otu1 > 0] <- 1
#transform it
t2 <- as.data.frame(t(otu1))
# add a column for group
m3 <- subset(m2, select = c("sample", "population"))
t2 <- merge(m3, t2, by.x = "sample", by.y = "row.names")
t3 <- aggregate(t2[3:length(t2)], by=list(t2$population), "sum")
rownames(t3) <- t3$Group.1
t3$Group.1 <- NULL

#make the sample size
t4 <- plyr::count(m2$population)

t5 <- as.data.frame(matrix(ncol = length(t3), nrow = nrow(t4)))
rownames(t5)  <- t4$x
colnames(t5) <- colnames(t3)

t5[,1:length(t5)] <- t4$freq

## counts: number of individuals in the population the OTU was found
allele.counts <- as.matrix(t3)

## sample_sizes: number of individuals in the population
sample.sizes <- as.matrix(t5)

## D: Pairwise geographic distance between populations
geom <- read.csv("geographic_distance_matrix_in_km.csv", row.names = 1, header = TRUE)
geom <- as.matrix(geom)
#normalize it by sd
geom <- geom/sd(geom)

## E: Pairwise ecological distance(s). Make them into a list.
geology_bc <- read.csv("geology_bray-curtis.csv", row.names = 1, header = TRUE)
geology_bc <- as.matrix(geology_bc)

env1 <- read.csv("environmental_metadata.csv", header = TRUE, row.names = 1)
env1 <- subset(env1, select = c("ph", "total.nitrogen_k_per_kg", "exchangeable.sodium.percent"))
# run the pca
bc.pca <- labdsv::pca(env1)
#create a distance matrix from PC1
soil_pc1 <- as.matrix(dist(bc.pca$scores[,1], upper=TRUE, diag=TRUE))

genetic_fst <- read.csv("population_FST.csv", row.names = 1, header = TRUE)
genetic_fst <- as.matrix(genetic_fst)
#normalize it by sd
genetic_fst <- genetic_fst/sd(genetic_fst)

eco_list <- list(geology_bc, soil_pc1, genetic_fst)

#Call the Markov chain Monte Carlo for the beta binomial model
MCMC_BB(
  counts = allele.counts,
  sample_sizes = sample.sizes,
  D = geom,
  E = eco_list,
  k = 14,
  loci = 2711,
  delta = 0.000001,  #
  aD_stp = 0.01, 
  aE_stp = 0.1, 
  a2_stp = 0.03, 
  thetas_stp = 0.1,
  phi_stp = 30, 
  mu_stp = 0.3, 
  ngen = 10000000,  
  printfreq = 10000,
  savefreq = 10000,
  samplefreq = 100,
  directory =  "/bedassle/",
  prefix = "r1_",
  continue = FALSE,
  continuing.params = NULL
)

## run below code on local machine
plot_all_acceptance_rates('r1_MCMC_output1.Robj') #shown in viewer

#make a vector popuilation names
pop_vector <- as.vector(rownames(t5))
plot_all_trace('r1_MCMC_output1.Robj', percent.burnin = 0, thinning = 1, population.names = pop_vector)

posterior.predictive.sample(MCMC.output = "r1_MCMC_output1.Robj", posterior.predictive.sample.size = 1000, output.file = "post1000", prefix = "r1_")
plot_posterior_predictive_samples("r1_post1000.Robj", save.figure = TRUE, figure.name = "r1_post1000_1.jpg")

#use a burn-in of the first 8 million
plot_all_trace('r1_MCMC_output1.Robj', percent.burnin = 80, thinning = 1, population.names = pop_vector)

show(load("r1_MCMC_output1.Robj"))

#to correct the aE : aD ratios, multiply by sd(geom) and sd(genetic_fst)
geom <- read.csv("geographic_distance_matrix_in_km.csv", row.names = 1, header = TRUE)
geom <- as.matrix(geom)

genetic_fst <- read.csv("population_FST.csv", row.names = 1, header = TRUE)
genetic_fst <- as.matrix(genetic_fst)

#correct the aE : aD ratios
#1 is geology
geology_bray_curtis_ratio <- ((as.vector(aE[1,])* sd(geom))/(aD))
#2 is soil pc1
soil_pca_pc1_ratio <- ((as.vector(aE[2,])* sd(geom))/(aD))
#3 is genetic fst
genetic_fst_ratio <- ((as.vector(aE[3,])* sd(geom))/(aD * sd(genetic_fst)))

#calculate output after 80% burn-in
library(coda)
all.param = as.mcmc(cbind(geology_bray_curtis_ratio[-(1:80000)], soil_pca_pc1_ratio[-(1:80000)], genetic_fst_ratio[-(1:80000)]))
s1 <- summary(all.param)
s2 <- as.data.frame(s1$statistics)
s2 <- s2[1:3,1:2]
s3 <- as.data.frame(HPDinterval(all.param, prob=0.95))
s3 <- cbind(s2, s3)
rownames(s3) <- c("geology_bray_curtis_ratio", "soil_pca_pc1_ratio", "genetic_fst_ratio")
s3 <- round(s3, digits = 5)
