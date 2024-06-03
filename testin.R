library(coda)
library(Hmsc)
load("examples/basic_example/mapp_post_test/other_test/models/models_thin_100_samples_250_chains_4.Rdata")
test <- convertToCodaObject(m_list$longterm)
traceplot(test$Beta)
plot(test$Beta)
par(mfrow = c(2, 4))
plot(test$Beta[,1:3])
plot(test$Gamma[,1:3])
plot(test$V[,1:3])
plot(test$Sigma[,1:3])
plot(test$Eta[[1]])
plot(test$Lambda[[1]])
plot(test$Omega[[1]][,1:3])
plot(test$Alpha[[1]][,1:3])
plot(test$Psi[[1]][,1:3])
plot(test$Delta[[1]])

plot(test$Beta[,1])
mean(test$Beta[[1]][,1])
sd(test$Beta[[1]][,1])

#Create a sequence of 1000 x values based on population mean and standard deviation
x <- seq(-4, 4, length = 1000) * sd(test$Beta[[1]][,1])+ mean(test$Beta[[1]][,1])

#create a vector of values that shows the height of the probability distribution
#for each value in x
y <- dnorm(x, mean(test$Beta[[1]][,1]), sd(test$Beta[[1]][,1]))

plot(x,y, type = "l", lwd = 2)




autocorr(test$Beta, lags = c(0, 1, 5, 10, 50), relative=TRUE)


sp <- colnames(m_list$longterm$YData)
var <- c("Intercept",(unlist(strsplit(as.character(m_list$longterm$XFormula)[2]," \\+ "))))

library(magrittr)
library(dplyr)
library(ggplot2)
library(tidyr)
plot_tmp_beta <- NULL
plot_var_list <- NULL
for (var_tmp in var){
    plot_tmp_beta <- NULL
    for(chain in 1:4){
        tmp_beta <- test$Beta[[chain]]
        tmp_beta <- as.data.frame(tmp_beta)
        tmp_beta_var <- tmp_beta%>%
            select(contains(var_tmp))
        melted_tmp_beta_var <- reshape2::melt(tmp_beta_var)
        melted_tmp_beta_var$chain <- as.factor(chain)
        melted_tmp_beta_var_splt <- tidyr::separate(melted_tmp_beta_var,variable,c("variable","species"),sep=",")
        plot_tmp_beta <- rbind(plot_tmp_beta,melted_tmp_beta_var_splt)
    }
    plot_var_list[[var_tmp]] <- ggplot(plot_tmp_beta)+
        geom_jitter(aes(x= value, y=species, color=chain),width=0, height = 0.12, alpha=0.4)+
        labs(title = paste0("Beta distribution - ",var_tmp))+
        theme_light()
}
plot_var_list$Intercept


library(jsonify)
library(Hmsc)

load("examples/basic_example/mapp_post_test/other_test/models/unfitted_models_noran.RData")
m_list_noran <- m_list
load("examples/basic_example/mapp_post_test/other_test/models/unfitted_models_hurdle.RData")
load("examples/basic_example/mapp_post_test/other_test/models/unfitted_models.RData")


test_hpc <-  from_json(readRDS(file = "examples/basic_example/mapp_post_test/other_test/models/models_thin_100_samples_250_chains_4_longterm_post_file.rds")[[1]])
test_hpc_noran <- from_json(readRDS(file = 'examples/basic_example/mapp_post_test/other_test/models/models_thin_100_samples_250_chains_4_longterm_post_file_noran.rds')[[1]])
test_hpc_mabuscaleX <- from_json(readRDS(file = 'examples/basic_example/mapp_post_test/other_test/models/models_thin_100_samples_250_chains_4_mabuscaleX_post_file_hurdle.rds')[[1]])
test_hpc_mabuscaleXY <- from_json(readRDS(file = 'examples/basic_example/mapp_post_test/other_test/models/models_thin_100_samples_250_chains_4_mabuscaleXY_post_file_hurdle.rds')[[1]])
test_hpc_mabu <- from_json(readRDS(file = 'examples/basic_example/mapp_post_test/other_test/models/models_thin_100_samples_250_chains_4_mabu_post_file_hurdle.rds')[[1]])
test_hpc_mpa <- from_json(readRDS(file = 'examples/basic_example/mapp_post_test/other_test/models/models_thin_100_samples_250_chains_4_mpa_post_file_hurdle.rds')[[1]])


postList <- test_hpc[1:4]
postList_noran <- test_hpc_noran[1:4]
postList_mpa <- test_hpc_mpa[1:4]
postList_mabu <- test_hpc_mabu[1:4]
postList_mabuscaleX <- test_hpc_mabuscaleX[1:4]
postList_mabuscaleXY <- test_hpc_mabuscaleXY[1:4]


cat(sprintf("fitting time %.1f sec\n", test_hpc[[4+1]]))
cat(sprintf("fitting time %.1f sec\n", test_hpc_noran[[4+1]]))
cat(sprintf("fitting time %.1f sec\n", test_hpc_mpa[[4+1]]))
cat(sprintf("fitting time %.1f sec\n", test_hpc_mabu[[4+1]]))
cat(sprintf("fitting time %.1f sec\n", test_hpc_mabuscaleX[[4+1]]))
cat(sprintf("fitting time %.1f sec\n", test_hpc_mabuscaleXY[[4+1]]))

fitTF <- importPosteriorFromHPC(m_list$longterm, postList, 250, 100, ceiling(0.5*5*1))
fitTF_noran <- importPosteriorFromHPC(m_list_noran$longterm, postList_noran, 250, 100, ceiling(0.5*5*1))
fitF_mpa <- importPosteriorFromHPC(m_list_hurdle$mpa, postList_mpa, 250, 100, ceiling(0.5*5*1))
fitF_mabu <- importPosteriorFromHPC(m_list_hurdle$mabu, postList_mabu, 250, 100, ceiling(0.5*5*1))
fitF_mabuscaleX <- importPosteriorFromHPC(m_list_hurdle$mabuscaleX, postList_mabuscaleX, 250, 100, ceiling(0.5*5*1))
fitF_mabuscaleXY <- importPosteriorFromHPC(m_list_hurdle$mabuscaleXY, postList_mabuscaleXY, 250, 100, ceiling(0.5*5*1))

plotVariancePartitioning(fitTF, computeVariancePartitioning(fitTF), args.legend=list(x="bottomright")) #longterm HPC
plotVariancePartitioning(m_list$longterm, computeVariancePartitioning(m_list$longterm), args.legend=list(x="bottomright")) #longterm classic R
plotVariancePartitioning(fitTF_noran, computeVariancePartitioning(fitTF_noran), args.legend=list(x="bottomright")) #longterm HPC
plotVariancePartitioning(fitF_mpa, computeVariancePartitioning(fitF_mpa), args.legend=list(x="bottomright"))# longterm PA
plotVariancePartitioning(fitF_mabu, computeVariancePartitioning(fitF_mabu), args.legend=list(x="bottomright"))# longterm ABU
plotVariancePartitioning(fitF_mabuscaleX, computeVariancePartitioning(fitF_mabuscaleX), args.legend=list(x="bottomright"))# longterm ABU Y scaled
plotVariancePartitioning(fitF_mabuscaleXY, computeVariancePartitioning(fitF_mabuscaleXY), args.legend=list(x="bottomright"))# longterm ABU Y scaled


preds <- computePredictedValues(m_list$longterm)
test1 <- evaluateModelFit(hM=m_list$longterm, predY=preds)
test1$R2
