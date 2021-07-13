library("MutationalPatterns")
library("devtools")
setwd('~/helmsman')
install_github("carjed/musigtools")
library("musigtools")
raw <- read.table("subtype_count_matrix.txt", header=T, stringsAsFactors=F)
msm <- format_counts(raw, 'MutationalPatterns')
msm
signatures = read.table('SBS_signatures.csv',sep=',', header = TRUE)
row.names(signatures) = paste(substring(signatures$SubType,1,1), "[", signatures$Type, "]", substring(signatures$SubType,3,3), sep = "" )
sigmat = as.matrix(signatures[,3:61], "numeric")
sigmat
plot_96_profile(msm, ymax = .03)
cos_sim_signatures = cos_sim_matrix(msm, sigmat)
cos_sim_signatures
fit_res = fit_to_signatures(msm, sigmat)
fit_res
hclust <- cluster_signatures(msm, method = 'complete')
plot(hclust)
fit_res = fit_to_signatures(msm, sigmat)
#select <- which(rowSums(fit_res$contribution) > 0)
plot_cosine_heatmap(fit_res$contribution, cluster_rows = T)

select <- which(rowSums(fit_res$contribution) > 15000)
plot_cosine_heatmap(cos_sim_signatures[,select], cluster_rows = T)
#now do de novo extraction
estimate <- nmf(msm, rank=2:5, method='brunet',nrun=10)
plot(estimate)
nmf_res <- extract_signatures(msm, rank = 3, nrun = 50)
colnames(nmf_res$signatures) <- c('Signature A','Signature B','Signature C')
rownames(nmf_res$contribution) <- c('Signature A','Signature B','Signature C')
plot_96_profile(nmf_res$signatures, condensed = T)
plot_contribution(nmf_res$contribution, nmf_res$signature, mode = 'relative')
plot_cosine_heatmap(nmf_res$contribution, cluster_rows = T)
