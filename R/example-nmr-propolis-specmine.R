library(specmine)

source("[change]/Propolis-NMR/scripts/propolis_metadata.R")

# reading data and metadata
prop.nmr.metadata.file = "[change]/Propolis-NMR/metadata/metadata_propolis.csv"
prop.nmr.data.folder = "[change]/Propolis-NMR/data/2010/"

get.metadata.agro(prop.nmr.data.folder, write.file = TRUE, file.name = prop.nmr.metadata.file)
prop.nmr.metadata = read_metadata(prop.nmr.metadata.file)
prop.nmr.metadata

peaks.lists = read_csvs_folder(prop.nmr.data.folder)

table(prop.nmr.metadata$seasons)
#au sm sp wi 
#15 16 15 13 

table(prop.nmr.metadata$agroregions)
#Highlands     Plain   Plateau 
#12        11        36 


# statistics on peaks per samples
sum(peaks_per_samples(peaks.lists))
# 25403
mean(peaks_per_samples(peaks.lists))
# 430.6

#all distinct frequencies
#all.freqs = get_overall_freq_list(peaks.lists)
#length(all.freqs)
# 943

prop.nmr.ds = group_peaks(peaks.lists, type = "nmr-peaks", 
                          metadata = prop.nmr.metadata, 
                          description = "NMR propolis",
                          label.x = "ppm", label.values = "intensity")
sum_dataset(prop.nmr.ds)

prop.nmr.ds$metadata$seasons = factor(prop.nmr.ds$metadata$seasons)
prop.nmr.ds$metadata$agroregions = factor(prop.nmr.ds$metadata$agroregions)

sum_dataset(prop.nmr.ds)
prop.nmr.ds$metadata$agroregions

# remove data peaks that will not be used in the analysis
get_data_values(prop.nmr.ds, seq(3.29, 3.31, 0.01))
get_data_values(prop.nmr.ds, seq(4.85, 5.00, 0.01))
get_data_values(prop.nmr.ds, 0)
to.remove = c(0, seq(3.29, 3.31, 0.01), seq(4.85, 5.00, 0.01))
prop.nmr.ds = remove_data_variables(prop.nmr.ds, to.remove)
sum_dataset(prop.nmr.ds)

# remove peak groups with less than 25% of values
nsamps = num_samples(prop.nmr.ds)
prop.nmr.ds = remove_variables_by_nas(prop.nmr.ds,  0.75*nsamps)
sum_dataset(prop.nmr.ds)
# 243 peaks, 59 samples

# missing values imputation - replace by low value
count_missing_values(prop.nmr.ds)
# 2756

prop.nmr.na = missingvalues_imputation(prop.nmr.ds, method="value", value = 0.00005)
sum_dataset(prop.nmr.na)
count_missing_values(prop.nmr.na)

prop.nmr.na$metadata$agroregions

# log transform
prop.nmr.log = transform_data(prop.nmr.na, method="log")

# auto scaling
prop.nmr.scal = specmine::scaling(prop.nmr.na, "auto")
prop.nmr.log.scal = specmine::scaling(prop.nmr.log, "auto")
sum_dataset(prop.nmr.log.scal)

# ANOVA analysis

# scaling does not change results of ANOVA - equal w/ ou w/o scaling
anova.prop.nmr.scal = aov_all_vars(prop.nmr.scal, "seasons")
anova.prop.nmr.scal[1:20,]
# pvalues     logs          fdr                      tukey
# 5.99 2.314264e-07 6.635587 5.623662e-05        sp-au; sp-sm; wi-sp
# 4.66 1.746928e-06 5.757725 2.122518e-04        sm-au; sp-sm; wi-sm
# 4.5  1.843384e-05 4.734384 1.493141e-03 sp-au; wi-au; sp-sm; wi-sm
# 4.45 2.617632e-05 4.582091 1.513982e-03        sp-au; sp-sm; wi-sm
# 4.55 3.115188e-05 4.506516 1.513982e-03        sm-au; sp-sm; wi-sm

anova.prop.nmr.log.scal = aov_all_vars(prop.nmr.log.scal,  "seasons")
anova.prop.nmr.log.scal[1:20,]
# pvalues      logs          fdr                             tukey
# 4.66 9.584707e-26 25.018421 2.329084e-23               sm-au; sp-sm; wi-sm
# 4.58 3.384728e-17 16.470476 4.112445e-15               sm-au; sp-sm; wi-sm
# 4.55 6.092443e-14 13.215209 4.934879e-12 sm-au; sp-au; wi-au; sp-sm; wi-sm
# 4.63 1.043963e-13 12.981315 6.342078e-12               sm-au; sp-sm; wi-sm
# 4.71 2.082569e-13 12.681401 1.012128e-11               sm-au; sp-sm; wi-sm

anova.prop.nmr.log.scal = aov_all_vars(prop.nmr.log.scal, "seasons")
anova.prop.nmr.log.scal[1:20,]

anova.prop.nmr.log.scal = aov_all_vars(prop.nmr.log.scal, "agroregions")
anova.prop.nmr.log.scal[1:20,]

# PCA analysis - used only scaling since PCA should always be done with scaled data

pca.prop.nmr.scal = pca_analysis_dataset(prop.nmr.scal)
summary(pca.prop.nmr.scal)
# Importance of components:
#                         PC1    PC2     PC3     PC4    PC5     PC6     PC7     PC8     PC9    PC10    PC11
# Standard deviation     7.1619 5.8524 4.14906 3.92434 3.4648 3.29923 3.16847 2.83972 2.57166 2.49603 2.26681
# Proportion of Variance 0.2111 0.1409 0.07084 0.06338 0.0494 0.04479 0.04131 0.03319 0.02722 0.02564 0.02115
# Cumulative Proportion  0.2111 0.3520 0.42287 0.48625 0.5356 0.58045 0.62176 0.65495 0.68216 0.70780 0.72895
pca_screeplot(pca.prop.nmr.scal)
# need 24 PCs to reach 90% of variance

pca.prop.nmr.log.scal = pca_analysis_dataset(prop.nmr.log.scal)
summary(pca.prop.nmr.log.scal)
# Importance of components:
#   PC1    PC2     PC3     PC4     PC5     PC6     PC7     PC8     PC9    PC10
# Standard deviation     6.8746 5.1918 4.83004 3.45039 3.15233 2.93118 2.63436 2.56172 2.51600 2.45092
# Proportion of Variance 0.1945 0.1109 0.09601 0.04899 0.04089 0.03536 0.02856 0.02701 0.02605 0.02472
# Cumulative Proportion  0.1945 0.3054 0.40142 0.45041 0.49130 0.52666 0.55522 0.58222 0.60828 0.63300
pca_screeplot(pca.prop.nmr.log.scal)
# need 32 PCs to reach 90% of variance

#pca.pairs.plot(autoscaling.propolis, "seasons", pca.analysis.result)
pca_scoresplot2D(prop.nmr.log.scal, pca.prop.nmr.log.scal, "seasons", ellipses = T)
pca_scoresplot2D(prop.nmr.log.scal, pca.prop.nmr.log.scal, "agroregions", ellipses = T)
pca_pairs_plot(prop.nmr.log.scal, pca.prop.nmr.log.scal, "seasons", pcas = 1:3)
pca_pairs_plot(prop.nmr.log.scal, pca.prop.nmr.log.scal, "agroregions", pcas = 1:3)

# MACHINE LEARNING

ml.prop.nmr = train_models_performance(prop.nmr.na, c("pls", "rf"), "seasons", 
                                      "repeatedcv", num.folds = 5, num.repeats = 5, 
                                       tunelength = 20)
ml.prop.nmr$performance
# Accuracy     Kappa AccuracySD   KappaSD
# pls 0.6736364 0.5619682 0.13388963 0.1790710
# rf  0.8319464 0.7750361 0.09338878 0.1248002

ml.prop.nmr.scal = train_models_performance(prop.nmr.scal, c("pls", "rf"), "seasons", 
                                       "repeatedcv", num.folds = 5, num.repeats = 10, tunelength = 10)
ml.prop.nmr.scal$performance
# Accuracy     Kappa AccuracySD   KappaSD
# pls 0.8207576 0.7596847  0.1015642 0.1368017
# rf  0.8342890 0.7781355  0.1023998 0.1369167

ml.prop.nmr.log = train_models_performance(prop.nmr.log, c("pls", "rf"), "seasons", 
                                                "repeatedcv", num.folds = 5, num.repeats = 10, tunelength = 10)
ml.prop.nmr.log$performance
# Accuracy     Kappa AccuracySD   KappaSD
# pls 0.8468182 0.7944511 0.10349007 0.1389097
# rf  0.8368298 0.7815238 0.08881368 0.1185289 

ml.prop.nmr.log.scal = train_models_performance(prop.nmr.log.scal, c("pls", "rf"), "seasons", 
                                            "repeatedcv", num.folds = 5, num.repeats = 5, 
                                            tunelength = 20)
ml.prop.nmr.log.scal$performance
# Accuracy     Kappa AccuracySD   KappaSD
# pls 0.8524242 0.8025702 0.09720762 0.1301584
# rf  0.8243124 0.7648034 0.09021085 0.1204002

# confusion matrices
ml.prop.nmr.log.scal$confusion.matrices$pls
ml.prop.nmr.log.scal$confusion.matrices$rf

# variable importance
summary_var_importance(ml.prop.nmr.log.scal, 10)

# pls plot
library(scatterplot3d)
pls.model = ml.prop.nmr.log.scal$final.models$pls
seasons = prop.nmr.log.scal$metadata$seasons
scatterplot3d(pls.model$scores[,1:3], color=as.integer(seasons), pch=17, xlab = "Component 1", ylab = "Component 2", zlab = "Component 3")
legend(-1.5, 2.5, levels(seasons), col = 1:4, pt.cex = 1.2, pch= 17)

# ml for regions
ml.prop.nmr.reg = train_models_performance(prop.nmr.log.scal, c("pls", "rf"), "agroregions", 
                                                "repeatedcv", num.folds = 5, num.repeats = 5, 
                                                tunelength = 20)
ml.prop.nmr.reg$performance
# Accuracy     Kappa AccuracySD   KappaSD
# pls 0.7297702 0.4583102 0.09627403 0.1905289
# rf  0.7079720 0.4152847 0.11712046 0.2105181

ml.prop.nmr.reg$confusion.matrices$pls
ml.prop.nmr.reg$confusion.matrices$rf

summary_var_importance(ml.prop.nmr.reg, 5)


# Correlation analysis
correl.prop.nmr.scal = correlations_dataset(prop.nmr.scal, method = "pearson")
heatmap(correl.prop.nmr.scal, col =  topo.colors(256))

# clustering
hc.prop.nmr.scal = clustering(prop.nmr.scal, method = "hc", distance = "euclidean", clustMethod="complete")
dendrogram_plot(prop.nmr.scal, hc.prop.nmr.scal)

kmeans.prop.nmr.scal = clustering(prop.nmr.scal, method = "kmeans", num.clusters = 4)
kmeans_plot(prop.nmr.scal, kmeans.prop.nmr.scal)

pca_kmeans_plot2D(prop.nmr.scal, pca.prop.nmr.scal, kmeans.result = kmeans.prop.nmr.scal, ellipses = T)

# Other ML models

ml2.prop.nmr.log.scal = train_models_performance(prop.nmr.log.scal, c("svmLinear"), "seasons", 
                                        "repeatedcv", num.folds = 5, num.repeats = 10, tunelength = 10)
ml2.prop.nmr.log.scal$performance
# Accuracy     Kappa AccuracySD   KappaSD
# svmLinear 0.8336131 0.7770978  0.1226224 0.1641851



# Feature selection

fs.prop.nmr.log.scal = feature_selection(prop.nmr.log.scal, "seasons", method="rfe", 
                                             functions = rfFuncs, validation = "repeatedcv", 
                                             repeats = 10, number = 5, subsets = 2^(1:6))
fs.prop.nmr.log.scal
# Variables Accuracy  Kappa AccuracySD KappaSD Selected
# 2   0.6254 0.4952    0.09209  0.1207         
# 4   0.6950 0.5898    0.09157  0.1212         
# 8   0.7388 0.6489    0.11397  0.1525         
# 16   0.8056 0.7403    0.13762  0.1837         
# 32   0.8265 0.7680    0.10957  0.1460         
# 64   0.8319 0.7751    0.11078  0.1478        *
# 243   0.8236 0.7641    0.10395  0.1384         

#The top 5 variables (out of 64):
#  X4.66, X4.58, X4.38, X4.5, X4.17

plot(fs.prop.nmr.log.scal, type=c("g","o"))

