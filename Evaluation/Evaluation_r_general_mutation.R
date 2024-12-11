#!/usr/bin/env Rscript
#
# Supplementary Figures and statistics
#

library(ggplot2)
library(cowplot)
library(ggrepel)
library(stringr)
library(RColorBrewer)
library(Biostrings)
library(reshape2)

################################################################################
# Create Output Scaffolding
################################################################################

base_dir <- '/Users/utente/Desktop/Varcovid/GenSeq/Model/Generattion/Sequence_Generated/SpikeGTP'

s <- ifelse(!dir.exists(file.path(base_dir, 'fig')), dir.create(file.path(base_dir, 'fig')), FALSE)
s <- ifelse(!dir.exists(file.path(base_dir, 'fig/supp')), dir.create(file.path(base_dir, 'fig/supp')), FALSE)

s <- ifelse(!dir.exists(file.path(base_dir, 'calc')), dir.create(file.path(base_dir, 'calc')), FALSE)
s <- ifelse(!dir.exists(file.path(base_dir, 'calc/supp')), dir.create(file.path(base_dir, 'calc/supp')), FALSE)

################################################################################
# Supplmentary Figure -- Pairwise Frequency Correlations (Train; Test; Predict)
################################################################################

df <- read.csv(file.path(base_dir, 'calc/merged_substitution_info.csv'))
df$location <- sapply(df$substitution, function(x) { 
  as.numeric(substr(x, 2, str_length(x)-1)) })
df$group <- factor(df$group, levels=c('Train', 'Test', 'Predict'), ordered=T)

for (l in list(c('Train', 'Test'), c('Train', 'Predict'), c('Test', 'Predict'))) {
  f1 <- l[[1]]
  f2 <- l[[2]]
  
  common.subs.df <- as.data.frame.matrix(table(df$substitution, df$group))
  common.subs <- rownames(common.subs.df)[common.subs.df[,f1] == 1
                                          & common.subs.df[,f2] == 1]
  f1.sub <- subset(df, group == f1)
  f2.sub <- subset(df, group == f2)
  
  merged.df <- data.frame(
    substitution = common.subs,
    f1 = f1.sub$frequency[match(common.subs, f1.sub$substitution)],
    f2 = f2.sub$frequency[match(common.subs, f2.sub$substitution)]
  )
  merged.df$f1.scaled <- merged.df$f1 / sum(merged.df$f1)
  merged.df$f2.scaled <- merged.df$f2 / sum(merged.df$f2)
  
  ggplot(merged.df, aes(x=log10(f1), y=log10(f2))) +
    geom_point() +
    geom_smooth(method='lm', formula= y~x) +
    theme_cowplot() +
    background_grid()
  ggsave(file.path(base_dir, paste0('fig/supp', f1, '_', f2, '_corr.png')), width=8, height=7.5)
  
  print(paste0('=====', f1, ' - ', f2, ' ====='))
  print(cor.test(merged.df$f1, merged.df$f2, method = 'spearman'))
  print('==========')
}

################################################################################
# Supplmentary Calculation -- Comparison of Test Accuracy Against Other Methods
################################################################################

####
# LLM Prediction
####
set.seed(42)

model.0.predictions <- subset(df, group=='Predict')$substitution
model.0.predictions <- model.0.predictions[
  ! model.0.predictions %in% subset(df, group=='Train')$substitution]

sum(model.0.predictions %in% df[df$group == 'Test',]$substitution)

####
# Comparator Models Setup
####
p_start = 14
p_end = 700
valid.aa <- 'ARNDCQEGHILKMFPSTWYV'
valid.aa.vec <- strsplit(valid.aa, '')[[1]]
locs.to.substitute <- seq(from=p_start, to=p_end, by=1)

primary.rbm.seq <- str_sub(
  df[match(seq(from=p_start, to=p_end, by=1), df$location),]$substitution,
  start=1, end=1)
names(primary.rbm.seq) <- seq(from=p_start, to=p_end, by=1)

####
# Model 1: Uniform Sampling and Mutation
####
set.seed(42)

n.samples <- length(model.0.predictions)
model.1.predictions <- c()
while(length(model.1.predictions) < n.samples) {
  # choose location
  tmp.loc <- sample(locs.to.substitute, 1)
  
  # choose amino acid change
  tmp.loc.chr <- as.character(tmp.loc)
  tmp.sub <- sample(valid.aa.vec[
    ! valid.aa.vec %in% primary.rbm.seq[tmp.loc.chr]], 1)
  tmp.variant <- paste0(primary.rbm.seq[tmp.loc.chr], tmp.loc.chr, tmp.sub)
  
  if (tmp.variant %in% subset(df, group == 'Train')$substitution) {
    next;
  }
  model.1.predictions <- c(
    model.1.predictions,
    paste0(primary.rbm.seq[tmp.loc.chr], tmp.loc.chr, tmp.sub))
}

sum(model.1.predictions %in% df[df$group == 'Test',]$substitution)

####
# Model 2: Weighted Sample from Existing Variant Positions
####
set.seed(42)

train.locs <- unique(subset(df, group=='Train')$location)
train.locs <- train.locs[train.locs %in% as.character(seq(from=p_start, to=p_end, by=1))]
weight.vec <- c()
for (l in train.locs) {
  weight.vec <- c(
    weight.vec,
    sum(subset(df, group=='Train')$location == l) / nrow(subset(df, group=='Train'))
  )
}

n.samples <- length(model.0.predictions)
model.2.predictions <- c()
while(length(model.2.predictions) < n.samples) {
  # choose location
  tmp.loc <- sample(train.locs, 1, prob = weight.vec)
  
  # choose amino acid change
  tmp.loc.chr <- as.character(tmp.loc)
  tmp.sub <- sample(valid.aa.vec[
    ! valid.aa.vec %in% primary.rbm.seq[tmp.loc.chr]], 1)
  tmp.variant <- paste0(primary.rbm.seq[tmp.loc.chr], tmp.loc.chr, tmp.sub)
  
  if (tmp.variant %in% subset(df, group == 'Train')$substitution) {
    next;
  }
  model.2.predictions <- c(
    model.2.predictions,
    paste0(primary.rbm.seq[tmp.loc.chr], tmp.loc.chr, tmp.sub))
}


sum(model.2.predictions %in% df[df$group == 'Test',]$substitution)
sum(model.2.predictions %in% df[df$group == 'Test',]$substitution
    & !model.2.predictions %in% df[df$group == 'Train',]$substitution)

################################################################################
# Supplmentary Calculation -- Number needed to Simulate (NNS)
################################################################################

n.nns.trials <- 15000
set.seed(42)

test.only.substitutions <- subset(df, group == 'Test')$substitution
test.only.substitutions <- test.only.substitutions[
  !test.only.substitutions %in% subset(df, group == 'Train')$substitution
]

# Model 0: LLM, bootstrapped sampling

run.model.0 <- function() {
  n.trials <- 1
  while(T) {
    tmp <- sample(model.0.predictions, 1)
    if (tmp %in% test.only.substitutions) { return (n.trials) }
    n.trials <- n.trials + 1
  }
}

model.0.nns.trials <- c()
for (i in 1:n.nns.trials) {
  model.0.nns.trials <- c(model.0.nns.trials, run.model.0())
}

# Model 1: random sampling

run.model.1 <- function() {
  n.trials <- 1
  while(T) {
    # choose location
    tmp.loc <- sample(locs.to.substitute, 1)
    
    # choose amino acid change
    tmp.loc.chr <- as.character(tmp.loc)
    tmp.sub <- sample(valid.aa.vec[
      ! valid.aa.vec %in% primary.rbm.seq[tmp.loc.chr]], 1)
    tmp.variant <- paste0(primary.rbm.seq[tmp.loc.chr], tmp.loc.chr, tmp.sub)
    
    if (tmp.variant %in% test.only.substitutions) { return (n.trials) }
    n.trials <- n.trials + 1
  }
}

model.1.nns.trials <- c()
for (i in 1:n.nns.trials) {
  model.1.nns.trials <- c(model.1.nns.trials, run.model.1())
}

# Model 2: weighted sampling
run.model.2 <- function() {
  n.trials <- 1
  while(T) {
    # choose location
    tmp.loc <- sample(train.locs, 1, prob = weight.vec)
    
    # choose amino acid change
    tmp.loc.chr <- as.character(tmp.loc)
    tmp.sub <- sample(valid.aa.vec[
      ! valid.aa.vec %in% primary.rbm.seq[tmp.loc.chr]], 1)
    tmp.variant <- paste0(primary.rbm.seq[tmp.loc.chr], tmp.loc.chr, tmp.sub)
    
    if (tmp.variant %in% test.only.substitutions) { return (n.trials) }
    n.trials <- n.trials + 1
  }
}

model.2.nns.trials <- c()
for (i in 1:n.nns.trials) {
  model.2.nns.trials <- c(model.2.nns.trials, run.model.2())
}

# Manually create the BLOSUM62 matrix
PAM30 <- matrix(c(
  6,  -7,  -4,  -5,  -7,  -3,  -4,  -6,  -5,  -7,  -7,  -7,  -6,  -9,  -2,   1,  -1, -16,  -9,  -2,
  -7,   8,  -6, -13, -10,   0,  -6,  -8,  -3,  -9,  -9,   1,  -7, -11,  -6,  -1,  -4,  -9,  -9,  -9,
  -4,  -6,   8,   2, -11,  -2,  -2,  -6,   0,  -8, -10,  -2,  -9, -11,  -7,   1,  -2, -17,  -8,  -9,
  -5, -13,   2,   9, -14,  -4,   3,  -4,  -2, -11, -13,  -6, -12, -16,  -9,  -2,  -5, -18, -13, -10,
  -7, -10, -11, -14,   9, -14, -15,  -9, -11,  -4,  -5, -14,  -7,  -5, -10,  -7,  -9, -15,  -4,  -5,
  -3,   0,  -2,  -4, -14,   8,  -1,  -7,   1, -10,  -9,   0,  -4, -13,  -5,  -2,  -4,  -9, -10,  -9,
  -4,  -6,  -2,   3, -15,  -1,   8,  -7,  -5,  -9, -10,  -3,  -8, -15,  -6,  -4,  -5, -16, -11,  -9,
  -6,  -8,  -6,  -4,  -9,  -7,  -7,   9,  -8,  -9, -10,  -8, -10, -10,  -7,  -4,  -7, -17, -11,  -7,
  -5,  -3,   0,  -2, -11,   1,  -5,  -8,  13,  -8,  -9,  -4,  -4,  -8,  -5,  -2,  -3,  -8,   2,  -8,
  -7,  -9,  -8, -11,  -4, -10,  -9,  -9,  -8,   8,   2,  -9,   2,   0,  -9,  -7,  -5, -10,  -4,   5,
  -7,  -9, -10, -13,  -5,  -9, -10, -10,  -9,   2,   9,  -9,   6,   4,  -9,  -8,  -5,  -5,  -6,   3,
  -7,   1,  -2,  -6, -14,   0,  -3,  -8,  -4,  -9,  -9,   8,  -6, -14,  -6,  -2,  -5, -12, -11,  -9,
  -6,  -7,  -9, -12,  -7,  -4,  -8, -10,  -4,   2,   6,  -6,  10,   0,  -9,  -7,  -5,  -6,  -5,   0,
  -9, -11, -11, -16,  -5, -13, -15, -10,  -8,   0,   4, -14,   0,  11, -12,  -9,  -8,  -4,   2,  -4,
  -2,  -6,  -7,  -9, -10,  -5,  -6,  -7,  -5,  -9,  -9,  -6,  -9, -12,  10,  -2,  -4, -15, -11,  -8,
  1,  -1,   1,  -2,  -7,  -2,  -4,  -4,  -2,  -7,  -8,  -2,  -7,  -9,  -2,   6,   1, -10,  -4,  -6,
  -1,  -4,  -2,  -5,  -9,  -4,  -5,  -7,  -3,  -5,  -5,  -5,  -5,  -8,  -4,   1,   7, -12,  -5,  -4,
  -16,  -9, -17, -18, -15,  -9, -16, -17,  -8, -10,  -5, -12,  -6,  -4, -15, -10, -12,  12,  -5, -17,
  -9,  -9,  -8, -13,  -4, -10, -11, -11,   2,  -4,  -6, -11,  -5,   2, -11,  -4,  -5,  -5,   9,  -5,
  -2,  -9,  -9, -10,  -5,  -9,  -9,  -7,  -8,   5,   3,  -9,   0,  -4,  -8,  -6,  -4, -17,  -5,   8
), nrow = 20, byrow = TRUE)

rownames(PAM30) <- colnames(PAM30) <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")

run.model.3 <- function() {
  n.trials <- 1
  while(TRUE) {
    # choose location
    tmp.loc <- sample(train.locs, 1, prob = weight.vec)
    
    # choose amino acid change
    tmp.loc.chr <- as.character(tmp.loc)
    tmp.loc.aa <- primary.rbm.seq[tmp.loc.chr]
    affinity.subset <- !names(BLOSUM62[tmp.loc.aa, ]) %in% tmp.loc.aa
    tmp.prob.vec <- BLOSUM62[tmp.loc.aa, ][affinity.subset] -
      (min(BLOSUM62[tmp.loc.aa, ][affinity.subset]) - 1)
    tmp.prob.vec <- tmp.prob.vec / sum(tmp.prob.vec)
    
    tmp.sub <- sample(names(BLOSUM62[tmp.loc.aa, ])[affinity.subset],
                      prob = tmp.prob.vec, 1)
    tmp.variant <- paste0(tmp.loc.aa, tmp.loc.chr, tmp.sub)
    
    if (tmp.variant %in% test.only.substitutions) { return (n.trials) }
    n.trials <- n.trials + 1
  }
}

model.3.nns.trials <- c()
for (i in 1:n.nns.trials) {
  model.3.nns.trials <- c(model.3.nns.trials, run.model.3())
}

# Print the results
model.3.nns.trials



# plot and run statistics
plot.df <- rbind(
  data.frame(
    'n'=model.0.nns.trials,
    'model'='LLM'
  ),
  data.frame(
    'n'=model.1.nns.trials,
    'model'='Model 1'
  ),
  data.frame(
    'n'=model.2.nns.trials,
    'model'='Model 2'
  ),
  data.frame(
    'n'=model.3.nns.trials,
    'model'='Model 3'
  )
)
plot.df$model <- factor(plot.df$model,
                        levels=c('Model 1', 'Model 2', 'Model 3', 'LLM'),
                        ordered = T)

ggplot(plot.df, aes(y=n, x=model, fill=model)) +
  geom_violin() +
  geom_boxplot(width=0.1, fill='#DDDDDD', alpha=0.5, outlier.shape = NA) +
  scale_fill_manual(values=rev(brewer.pal(4, 'Set2'))) +
  theme_cowplot()
ggsave(file.path(base_dir, 'fig/model_comparison_violins.png'), width=8.5, height=8)

ggplot(plot.df, aes(y=log10(n), x=model, fill=model)) +
  geom_violin() +
  geom_boxplot(width=0.1, fill='#DDDDDD', alpha=0.5, outlier.shape = NA) +
  scale_fill_manual(values=rev(brewer.pal(4, 'Set2'))) +
  theme_cowplot()
ggsave(file.path(base_dir, 'fig/model_comparison_violins_logscale.png'), width=8.5, height=6)

# test differences
model.names <- c('Model 1', 'Model 2', 'Model 3', 'LLM')
p.val.matrix <- matrix(0, nrow=length(model.names), ncol=length(model.names))
stat.matrix <- matrix(0, nrow=length(model.names), ncol=length(model.names))
for (i in 1:length(model.names)) {
  for (j in 1:length(model.names)) {
    tst <- t.test(
      log10(subset(plot.df, model==model.names[i])$n),
      log10(subset(plot.df, model==model.names[j])$n))
    p.val.matrix[i, j] <- tst$p.value
    stat.matrix[i, j] <- tst$statistic
  }
}
rownames(p.val.matrix) <- model.names
colnames(p.val.matrix) <- model.names

rownames(stat.matrix) <- model.names
colnames(stat.matrix) <- model.names

print('All done!')