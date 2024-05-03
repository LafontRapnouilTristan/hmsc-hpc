library(jsonify)
library(Hmsc)

importFromHPC <- from_json(readRDS(file = "examples/basic_example/mapp_post_test/other_test/models/models_thin_1_samples_5_chains_4_longterm_post_file.rds")[[1]])
postList <- importFromHPC[1:nChains]
cat(sprintf("fitting time %.1f sec\n", importFromHPC[[nChains+1]]))

nSamples <- 5
transient <- ceiling(0.5*samples*thin)
fitTF = importPosteriorFromHPC(m_list$longterm, postList, nSamples, thin, transient)



plotVariancePartitioning(fitTF, computeVariancePartitioning(fitTF), args.legend=list(x="bottomright"))

?from_json
??fom_json

readRDS(file = "examples/basic_example/mapp_post_test/models_thin_10_samples_250_chains_4_longterm_post_file.rds")

readRDS(file = "examples/basic_example/mapp_post_test/models_thin_10_samples_250_chains_4_null_post_file.rds")

library(Hmsc)
readRDS(file = "examples/basic_example/post_file.rds")


nChains <- 4
importFromHPC <- from_json(readRDS(file = "examples/basic_example/post_file.rds")[[1]])
postList <- importFromHPC[1:nChains]
cat(sprintf("fitting time %.1f sec\n", importFromHPC[[nChains+1]]))
