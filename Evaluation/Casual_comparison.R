#!/usr/bin/env Rscript
#
# Supplementary Figures and Statistics with Multiple LLMs using PAM30
#

library(ggplot2)
library(cowplot)
library(ggrepel)
library(stringr)
library(RColorBrewer)
library(reshape2)

################################################################################
# Define LLMs and Their Corresponding Directories
################################################################################

llm_dirs <- list(
  'SARITA S' = '/Users/utente/Desktop/Varcovid/GenSeq/Model/Generattion/Sequence_Generated/run_clm_256_s',
  'SARITA M' = '/Users/utente/Desktop/Varcovid/GenSeq/Model/Generattion/Sequence_Generated/run_clm_256_m',
  'SARITA L' = '/Users/utente/Desktop/Varcovid/GenSeq/Model/Generattion/Sequence_Generated/run_clm_256_l',
  'SARITA XL' = '/Users/utente/Desktop/Varcovid/GenSeq/Model/Generattion/Sequence_Generated/run_clm_256_xl'
)

################################################################################
# Initialize Variables
################################################################################

# Base LLM names (keys of the llm_dirs list)
llm_names <- names(llm_dirs)

# Initialize lists to store data and predictions
df_llm_list <- list()
llm_predictions_list <- list()

# Set seed for reproducibility
set.seed(42)

################################################################################
# Process Each LLM
################################################################################

for (llm_name in llm_names) {
  base_dir <- llm_dirs[[llm_name]]
  
  # Create Output Scaffolding if necessary
  s <- ifelse(!dir.exists(file.path(base_dir, 'fig')), dir.create(file.path(base_dir, 'fig')), FALSE)
  s <- ifelse(!dir.exists(file.path(base_dir, 'fig/supp')), dir.create(file.path(base_dir, 'fig/supp')), FALSE)
  
  s <- ifelse(!dir.exists(file.path(base_dir, 'calc')), dir.create(file.path(base_dir, 'calc')), FALSE)
  s <- ifelse(!dir.exists(file.path(base_dir, 'calc/supp')), dir.create(file.path(base_dir, 'calc/supp')), FALSE)
  
  # Read the merged substitution info for this LLM
  df_llm <- read.csv(file.path(base_dir, 'calc/merged_substitution_info.csv'))
  df_llm$location <- sapply(df_llm$substitution, function(x) { 
    as.numeric(substr(x, 2, str_length(x)-1)) })
  df_llm$group <- factor(df_llm$group, levels=c('Train', 'Test', 'Predict'), ordered=T)
  
  # Store the dataframe for later use
  df_llm_list[[llm_name]] <- df_llm
  
  # Get LLM predictions
  model_0_predictions <- subset(df_llm, group=='Predict')$substitution
  model_0_predictions <- model_0_predictions[
    ! model_0_predictions %in% subset(df_llm, group=='Train')$substitution]
  
  # Store the predictions
  llm_predictions_list[[llm_name]] <- model_0_predictions
}

################################################################################
# Set Up Common Variables for Models
################################################################################

# Using the first LLM's data as the main data (assuming they share the same train/test splits)
df_main <- df_llm_list[[llm_names[1]]]

p_start = 14
p_end = 700
valid.aa <- 'ARNDCQEGHILKMFPSTWYV'
valid.aa.vec <- strsplit(valid.aa, '')[[1]]
locs.to.substitute <- seq(from=p_start, to=p_end, by=1)

primary.rbm.seq <- str_sub(
  df_main[match(seq(from=p_start, to=p_end, by=1), df_main$location),]$substitution,
  start=1, end=1)
names(primary.rbm.seq) <- seq(from=p_start, to=p_end, by=1)

# Define train.locs and weight.vec using df_main
train.locs <- unique(subset(df_main, group=='Train')$location)
train.locs <- train.locs[train.locs %in% as.character(seq(from=p_start, to=p_end, by=1))]
weight.vec <- c()
for (l in train.locs) {
  weight.vec <- c(
    weight.vec,
    sum(subset(df_main, group=='Train')$location == l) / nrow(subset(df_main, group=='Train'))
  )
}

################################################################################
# Comparator Models
################################################################################

#### Model 3: PAM30-Based Sampling (Rand PAM30)

# Manually create the PAM30 matrix
aa_order <- valid.aa.vec

PAM30_values <- c(
  6, -7, -4, -5, -7, -3, -4, -6, -5, -7, -7, -7, -6, -9, -2,  1, -1,-16, -9, -2,
  -7,  8, -6,-13,-10,  0, -6, -8, -3, -9, -9,  1, -7,-11, -6, -1, -4, -9, -9, -9,
  -4, -6,  8,  2,-11, -2, -2, -6,  0, -8,-10, -2, -9,-11, -7,  1, -2,-17, -8, -9,
  -5,-13,  2,  9,-14, -4,  3, -4, -2,-11,-13, -6,-12,-16, -9, -2, -5,-18,-13,-10,
  -7,-10,-11,-14,  9,-14,-15, -9,-11, -4, -5,-14, -7, -5,-10, -7, -9,-15, -4, -5,
  -3,  0, -2, -4,-14,  8, -1, -7,  1,-10, -9,  0, -4,-13, -5, -2, -4, -9,-10, -9,
  -4, -6, -2,  3,-15, -1,  8, -7, -5, -9,-10, -3, -8,-15, -6, -4, -5,-16,-11, -9,
  -6, -8, -6, -4, -9, -7, -7,  9, -8, -9,-10, -8,-10,-10, -7, -4, -7,-17,-11, -7,
  -5, -3,  0, -2,-11,  1, -5, -8, 13, -8, -9, -4, -4, -8, -5, -2, -3, -8,  2, -8,
  -7, -9, -8,-11, -4,-10, -9, -9, -8,  8,  2, -9,  2,  0, -9, -7, -5,-10, -4,  5,
  -7, -9,-10,-13, -5, -9,-10,-10, -9,  2,  9, -9,  6,  4, -9, -8, -5, -5, -6,  3,
  -7,  1, -2, -6,-14,  0, -3, -8, -4, -9, -9,  8, -6,-14, -6, -2, -5,-12,-11, -9,
  -6, -7, -9,-12, -7, -4, -8,-10, -4,  2,  6, -6, 10,  0, -9, -7, -5, -6, -5,  0,
  -9,-11,-11,-16, -5,-13,-15,-10, -8,  0,  4,-14,  0, 11,-12, -9, -8, -4,  2, -4,
  -2, -6, -7, -9,-10, -5, -6, -7, -5, -9, -9, -6, -9,-12, 10, -2, -4,-15,-11, -8,
  1, -1,  1, -2, -7, -2, -4, -4, -2, -7, -8, -2, -7, -9, -2,  6,  1,-10, -4, -6,
  -1, -4, -2, -5, -9, -4, -5, -7, -3, -5, -5, -5, -5, -8, -4,  1,  7,-12, -5, -4,
  -16, -9,-17,-18,-15, -9,-16,-17, -8,-10, -5,-12, -6, -4,-15,-10,-12, 12, -5,-17,
  -9, -9, -8,-13, -4,-10,-11,-11,  2, -4, -6,-11, -5,  2,-11, -4, -5, -5,  9, -5,
  -2, -9, -9,-10, -5, -9, -9, -7, -8,  5,  3, -9,  0, -4, -8, -6, -4,-17, -5,  8
)

PAM30 <- matrix(PAM30_values, nrow = 20, ncol = 20, byrow = TRUE)
rownames(PAM30) <- aa_order
colnames(PAM30) <- aa_order

# Model 3 Predictions (Rand PAM30)
n.samples <- length(subset(df_main, group=='Predict')$substitution)
model_3_predictions <- c()
while(length(model_3_predictions) < n.samples) {
  # choose location
  tmp.loc <- sample(train.locs, 1, prob = weight.vec)
  
  # choose amino acid change based on PAM30
  tmp.loc.chr <- as.character(tmp.loc)
  tmp.loc.aa <- primary.rbm.seq[tmp.loc.chr]
  affinity.subset <- !names(PAM30[tmp.loc.aa, ]) %in% tmp.loc.aa
  tmp.prob.vec <- PAM30[tmp.loc.aa, ][affinity.subset] -
    (min(PAM30[tmp.loc.aa, ][affinity.subset]) - 1)
  tmp.prob.vec <- tmp.prob.vec / sum(tmp.prob.vec)
  
  tmp.sub <- sample(names(PAM30[tmp.loc.aa, ])[affinity.subset],
                    prob = tmp.prob.vec, 1)
  tmp.variant <- paste0(tmp.loc.aa, tmp.loc.chr, tmp.sub)
  
  if (tmp.variant %in% subset(df_main, group == 'Train')$substitution) {
    next;
  }
  model_3_predictions <- c(
    model_3_predictions,
    tmp.variant)
}

################################################################################
# Number Needed to Simulate (NNS) Calculations
################################################################################

n.nns.trials <- 15000

test.only.substitutions <- subset(df_main, group == 'Test')$substitution
test.only.substitutions <- test.only.substitutions[
  !test.only.substitutions %in% subset(df_main, group == 'Train')$substitution
]

#### Model 3 NNS (Rand PAM30)
run_model_3 <- function() {
  n.trials <- 1
  while(TRUE) {
    # choose location
    tmp.loc <- sample(train.locs, 1, prob = weight.vec)
    
    # choose amino acid change based on PAM30
    tmp.loc.chr <- as.character(tmp.loc)
    tmp.loc.aa <- primary.rbm.seq[tmp.loc.chr]
    affinity.subset <- !names(PAM30[tmp.loc.aa, ]) %in% tmp.loc.aa
    tmp.prob.vec <- PAM30[tmp.loc.aa, ][affinity.subset] -
      (min(PAM30[tmp.loc.aa, ][affinity.subset]) - 1)
    tmp.prob.vec <- tmp.prob.vec / sum(tmp.prob.vec)
    
    tmp.sub <- sample(names(PAM30[tmp.loc.aa, ])[affinity.subset],
                      prob = tmp.prob.vec, 1)
    tmp.variant <- paste0(tmp.loc.aa, tmp.loc.chr, tmp.sub)
    
    if (tmp.variant %in% test.only.substitutions) { return (n.trials) }
    n.trials <- n.trials + 1
  }
}

model_3_nns_trials <- replicate(n.nns.trials, run_model_3())

#### LLM NNS Calculations

llm_nns_trials_list <- list()

for (llm_name in llm_names) {
  model_0_predictions <- llm_predictions_list[[llm_name]]
  
  run_model_0 <- function() {
    n.trials <- 1
    while(TRUE) {
      tmp <- sample(model_0_predictions, 1)
      if (tmp %in% test.only.substitutions) { return (n.trials) }
      n.trials <- n.trials + 1
    }
  }
  
  model_0_nns_trials <- replicate(n.nns.trials, run_model_0())
  
  # Store the results
  llm_nns_trials_list[[llm_name]] <- model_0_nns_trials
}

################################################################################
# Combine All Models for Plotting
################################################################################

# Prepare data frames for LLMs
llm_dfs <- lapply(llm_names, function(llm_name) {
  data.frame(
    'n' = llm_nns_trials_list[[llm_name]],
    'model' = llm_name
  )
})

# Models without data (add a star)
models_without_data <- c('Rand', 'Rand Fr', 'RITA S', 'RITA M', 'RITA L', 'RITA XL', 'SpikeGPT2')
models_without_data_df <- data.frame(
  'n' = NA,
  'model' = models_without_data
)

# Model with data (Rand PAM30)
model_3_df <- data.frame(
  'n' = model_3_nns_trials,
  'model' = 'Rand PAM30'
)

# Combine all data frames
plot.df <- rbind(
  models_without_data_df,
  model_3_df,
  do.call(rbind, llm_dfs)
)

# Update model names and factor levels
model.names <- c('Rand', 'Rand Fr', 'Rand PAM30', 'RITA S', 'RITA M', 'RITA L', 'RITA XL',
                 'SpikeGPT2', 'SARITA S', 'SARITA M', 'SARITA L', 'SARITA XL')
plot.df$model <- factor(plot.df$model,
                        levels = model.names,
                        ordered = TRUE)


# Add a column to indicate whether the model has data
plot.df$has_data <- !is.na(plot.df$n)

################################################################################
# Plotting
################################################################################

# Set color palette (you can adjust colors as needed)
color_palette <- c(
  'Rand' = '#1f77b4',
  'Rand Fr' = '#ff7f0e',
  'Rand PAM30' = '#2ca02c',
  'RITA S' = '#d62728',
  'RITA M' = '#9467bd',
  'RITA L' = '#8c564b',
  'RITA XL' = '#e377c2',
  'SpikeGPT2' = '#7f7f7f',
  'SARITA S' = '#bcbd22',
  'SARITA M' = '#17becf',
  'SARITA L' = '#aec7e8',
  'SARITA XL' = '#ff9896'
)

# Define the y-value for the star symbol (e.g., median of NNS values)
star_y_value <- median(model_3_nns_trials, na.rm = TRUE)

# Violin plot without log scale
p <- ggplot() +
  geom_violin(data = subset(plot.df, has_data),
              aes(y = n, x = model, fill = model)) +
  geom_boxplot(data = subset(plot.df, has_data),
               aes(y = n, x = model),
               width = 0.1, fill = '#DDDDDD', alpha = 0.5, outlier.shape = NA) +
  geom_point(data = subset(plot.df, !has_data),
             aes(x = model, y = star_y_value),
             shape = 8, size = 5, color = 'black') +
  scale_fill_manual(values = color_palette) +
  scale_x_discrete(limits = model.names) + 
  theme_cowplot() +
  ylab('Number Needed to Simulate (NNS)') +
  xlab('Model') +
  ggtitle('Comparison of NNS Across Models') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_blank(),
        legend.position = "none")  # Remove legend

# Save the plot
ggsave(file.path('/Users/utente/Desktop/Varcovid/GenSeq/Model/Generattion/Sequence_Generated/combined_figures', 'model_comparison_violins.png'), p, width = 12, height = 8)

# Violin plot with log scale
p_log <- ggplot() +
  geom_violin(data = subset(plot.df, has_data),
              aes(y = log10(n), x = model, fill = model)) +
  geom_boxplot(data = subset(plot.df, has_data),
               aes(y = log10(n), x = model),
               width = 0.1, fill = '#DDDDDD', alpha = 0.5, outlier.shape = NA) +
  geom_point(data = subset(plot.df, !has_data),
             aes(x = model, y = log10(star_y_value)),
             shape = 8, size = 5, color = 'black') +
  scale_fill_manual(values = color_palette) +
  scale_x_discrete(limits = model.names) + 
  theme_cowplot() +
  ylab('Log10 of NNS') +
  xlab('Model') +
  ggtitle('Comparison of NNS Across Models (Log Scale)') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_blank(),
        legend.position = "none")  # Remove legend

# Save the plot
ggsave(file.path('/Users/utente/Desktop/Varcovid/GenSeq/Model/Generattion/Sequence_Generated/combined_figures', 'model_comparison_violins_logscale.png'), p_log, width = 12, height = 8)

################################################################################
# Statistical Tests Between Models
################################################################################

# Only perform statistical tests between models with data
models_with_data <- c('Rand PAM30', llm_names)

p.val.matrix <- matrix(NA, nrow = length(model.names), ncol = length(model.names))
stat.matrix <- matrix(NA, nrow = length(model.names), ncol = length(model.names))
for (i in seq_along(model.names)) {
  for (j in seq_along(model.names)) {
    if (model.names[i] %in% models_with_data && model.names[j] %in% models_with_data) {
      tst <- t.test(
        log10(subset(plot.df, model == model.names[i])$n),
        log10(subset(plot.df, model == model.names[j])$n))
      p.val.matrix[i, j] <- tst$p.value
      stat.matrix[i, j] <- tst$statistic
    } else {
      p.val.matrix[i, j] <- NA
      stat.matrix[i, j] <- NA
    }
  }
}
rownames(p.val.matrix) <- model.names
colnames(p.val.matrix) <- model.names

rownames(stat.matrix) <- model.names
colnames(stat.matrix) <- model.names

# Print p-value matrix
print("P-Value Matrix:")
print(p.val.matrix)

# Print test statistic matrix
print("Test Statistic Matrix:")
print(stat.matrix)

print('All done!')


