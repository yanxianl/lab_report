## Import data and tidy --------------------------------------------------------
library(tidyverse) # tidy, manipulate and plot data
df <- read_csv("data/AqFl2_qPCR_ref_PI.csv", col_names = T, na = "") 
str(df)

## Tidy data
df <- df %>% 
  arrange(Ref_gene, desc(Diet)) %>% 
  mutate(Quantity = `^`(PE, -Cq)) %>% # calculate mRNA quantity
  group_by(Ref_gene, Diet)  %>%
  mutate(Mean = mean(Quantity), SD = sd(Quantity)) %>% # calculate average mRNA quantity and sd
  ungroup() %>%
  group_by(Ref_gene) %>%
  mutate(Quantity_std = scale(Quantity), # standardize mRNA quantity of each reference gene
         Quantity_nm = Quantity / Quantity[1]) # normalize mRNA quantity to the first sample

# Convert Sample_ID and Net_pen to character variables
df <- within(df, {
  Sample_ID <- as.character(Sample_ID)
  Net_pen <- as.character(Net_pen)
})

# Convert Diet and Sample_ID to factors, specifying the desired orders for plotting
df$Diet <- factor(df$Diet, levels = c("REF", "IM"))
df$Sample_ID <- factor(df$Sample_ID, levels = unique(df$Sample_ID)) 

## Exploratory data analysis (EDA) --------------------------------------------- 
# Descriptive statistics --------------------------------------------
# Convert dataframe to a wider form
library(data.table)
dt <- as.data.table(df)
df_wider <- dcast(dt, Sample_ID + Diet + Net_pen ~ Ref_gene, value.var = colnames(dt)[5:13])

# Descriptive statistics with *summarytools* --------------
library(summarytools) 
view(dfSummary(df_wider))

# Descriptive statistics with *skimr* ---------------------
library(skimr)
skim(df_wider) 

# EDA by visualization ----------------------------------------------
library(cowplot) #  a ggplot2 add-on to provide a publication-ready theme 

# Dot chart -----------------------------------------------
# Make a dataframe containing the grand mean of mRNA quantity 
gm <- df %>% 
  group_by(Ref_gene) %>%
  summarize(gmean = mean(Quantity))

# Make dot charts showing individual expression profile
p1 <- df %>% 
  ggplot(aes(x = Sample_ID, y = Quantity, group = 1)) +
    geom_point(aes(color = Diet)) +
    geom_line() +
    facet_wrap(~ Ref_gene, ncol = 1, scales = "free_y") +
    labs(y = "mRNA quantity", tag = "A") +
    scale_y_continuous(labels = scales::scientific) + 
    geom_hline(aes(yintercept = gmean), gm, color = "blue") + # The blue line shows the grand mean.
    theme_bw() 
p1

# Make dot charts showing normalized expression levels of all the reference genes in one plot
p2 <- df %>% 
  ggplot(aes(x = Sample_ID, y = Quantity_nm, color = Ref_gene)) +
    geom_point() +
    geom_line(aes(group = Ref_gene)) +
    labs(y = "normalized mRNA quantity", tag = "B") +
    theme_bw() +
    scale_fill_brewer(palette = "Dark2") 

p2

# Combine the plots
plot_grid(p1, p2, ncol = 1, align = "v", rel_heights = c(3, 2))

# Boxplot overlaying violin plot --------------------------
# Convert "Sample_ID"" to character, otherwise ggrepel won't work properly
df$Sample_ID <- as.character(df$Sample_ID)

# Define a function for identifying outliers
is_outlier <- function(x) {
  x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x)
}

# Make plots using standardized mRNA quantity
set.seed(1910)

df %>% 
  group_by(Ref_gene) %>%
  mutate(outlier = is_outlier(Quantity_std)) %>%
  ggplot(aes(x = Ref_gene, y = Quantity_std, 
             label = ifelse(outlier, Sample_ID, NA))) +
  geom_violin(trim = T) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_point(aes(fill = Diet, colour = Diet), size = 3, shape = 21, 
             position = position_jitterdodge(0.2)) +
  ggrepel::geom_text_repel(aes(colour = Diet), 
                           position = position_jitterdodge(0.2)) + # label outliers with Sample_ID
  ylab("standardized mRNA quantity") +
  theme_cowplot() +
  guides(colour = F)

# Bar plot ------------------------------------------------
df %>%
  ggplot(aes(x = Diet, y =  Mean, fill = Diet)) +
  geom_bar(stat = "identity", position = position_dodge(), colour = "black") +
  geom_errorbar(aes(ymin = Mean, ymax = Mean + SD), size = 0.3, width = 0.2) + # add error bar (sd)
  facet_wrap(~ Ref_gene, nrow = 1, scales = "free_y") +
  scale_y_continuous(limits = c(0, NA), 
                     expand = expand_scale(mult = c(0, 0.1))) + 
  ylab("mRNA quantity") +
  theme_cowplot()

# Correlation ---------------------------------------------
# Make a subset of data for correlation analysis
df1 <- df %>%
  select(Sample_ID, Ref_gene, Quantity) %>%
  spread(key = "Ref_gene", value = "Quantity")

# Correlation by chart.Correlation() ------------
library("PerformanceAnalytics")
chart.Correlation(df1[2:ncol(df1)], histogram = TRUE, pch = 19)

# Correlation by corrplot() ---------------------
# Compute a correlation matrix
corr <- cor(df1[2:ncol(df1)], method = "pearson", use = "complete.obs")

# Make the correlation plot
library(corrplot)
corrplot.mixed(corr, order = "hclust", tl.col = "black")

# heatmap -------------------------------------------------
# Make a numeric metrix for plotting heatmap
mat <- as.matrix(df1[2:ncol(df1)]) %>% 
  `row.names<-`(df1[[1]]) # Sample_ID as rownames

# Define color scheme
col <- colorRampPalette(c("blue", "white", "red"))(25)

# Static heatmap --------------------------------
library(gplots)
heatmap.2(mat,                     
          Rowv = TRUE,
          Colv = TRUE,
          distfun = function(x) dist(x,method = 'euclidean'),
          hclustfun = function(x) hclust(x,method = 'ward.D2'),
          ylab = "Sample ID",
          dendrogram = "column",
          density.info = "none",
          scale = "column",
          trace = "none",
          srtCol = 45,
          col = col)                       

# Interactive heatmap ---------------------------
library(d3heatmap)
d3heatmap(mat,                     
          Rowv = TRUE,
          Colv = TRUE,
          distfun = function(x) dist(x,method = 'euclidean'),
          hclustfun = function(x) hclust(x,method = 'ward.D2'),
          dendrogram = "column",
          density.info = "none",
          scale = "column",
          trace = "none",
          col = col)

## Evaluation of reference gene stability --------------------------------------
# Ranks by coefficient of variance (cv) -----------------------------
cv <- df %>% 
  group_by(Ref_gene) %>%
  summarize(Cq_max = max(Cq), 
            Cq_min = min(Cq), 
            Cq_range = max(Cq) - min(Cq), 
            Quantity_mean = mean(Quantity), 
            Quantity_SD = sd(Quantity),
            Quantity_CV = 100 * Quantity_SD / Quantity_mean) %>%
  mutate(CV_rank = dense_rank(Quantity_CV)) %>% # rank reference genes by CV
  arrange(Ref_gene)

# Ranks by F-statistic ----------------------------------------------
library(broom) # tidy statistical outputs
F_stat <- df %>% 
  group_by(Ref_gene) %>% 
  do(tidy(aov(Quantity ~ Diet, data = .))) %>% # run ANOVA and tidy statistical outputs
  na.omit() %>%
  select(Ref_gene, statistic) %>%
  rename(Quantity_F = statistic) %>%
  as.data.frame() %>%
  mutate(F_rank = dense_rank(Quantity_F)) %>% # rank reference genes by F statistic
  arrange(Ref_gene) 

## Summary of ranks -------------------------------------------------
# Merge tables
smr <- full_join(cv, F_stat, by = "Ref_gene")

# Format digits
smr$Quantity_mean <- formatC(smr$Quantity_mean, format = "e", digits = 1)
smr$Quantity_SD <- formatC(smr$Quantity_SD, format = "e", digits = 2)

library(kableExtra)
kable(smr, format = "html", digits = 2, align = "l") %>%
      kable_styling(bootstrap_options = c("striped", "hover"), 
                    font_size = 14, 
                    position = "left")

# Session info -----------------------------------------------------------------
sessionInfo()
