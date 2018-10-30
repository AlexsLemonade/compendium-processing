`%>%` <- dplyr::`%>%`
library(ggplot2)

#### samples -------------------------------------------------------------------

samples_df <- data.frame(
  num_samples = c(seq(100, 1000, 100), 5000, 10000, 20000, 40000),
  num_genes = rep(500, 14),
  run_time = c(
    0.048 + (0.002 * 100),
    0.141 + (0.002 * 200),
    0.295 + (0.002 * 300),
    0.535 + (0.002 * 400),
    0.825 + (0.002 * 500),
    1.284 + (0.002 * 600),
    1.587 + (0.002 * 700),
    2.158 + (0.002 * 800),
    2.717 + (0.002 * 900),
    3.371 + (0.002 * 1000),
    117.030 + (0.004 * 5000),
    472.211 + (0.005 * 10000),
    1742.242 + (0.011 * 20000),
    6735.273 + (0.012 * 40000)
  )
)

samples_df %>%
  ggplot(aes(x = num_samples, y = run_time)) +
  geom_point() +
  geom_smooth(method = "loess", formula = y ~ x) +
  theme_classic() +
  labs(x = "number of samples",
       y = "run time (s)",
       title = "KNN Impute (K=10)",
       subtitle = "n genes = 500")

#### genes ---------------------------------------------------------------------

genes_df <- data.frame(
  num_samples = rep(500, 6),
  num_genes = c(3000, 6000, 9000, 12000, 15000, 18000),
  run_time = c(
    6.962 + (0.011 * 500),
    13.917 + (0.021 * 500),
    20.366 + (0.031 * 500),
    27.177 + (0.043 * 500),
    34.242 + (0.054 * 500),
    40.276 + (0.061 * 500)
  )
)

genes_df %>%
  ggplot(aes(x = num_genes, y = run_time)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic() +
  labs(x = "number of genes",
       y = "run time (s)",
       title = "KNN Impute (K=10)",
       subtitle = "n samples = 500")

#### "sweep" -------------------------------------------------------------------

sweep_df <- data.frame(
  num_samples = rep(c(500, 1000, 10000, 25000), 2),
  num_genes = c(rep(3000, 4), rep(15000, 4)),
  run_time = c(
    6.311 + (0.010 * 500), 
    24.531 + (0.012 * 1000),
    2455.273 + (0.029 * 10000),
    15232.715 + (0.061 * 25000),
    30.947 + (0.053 * 500),
    122.254 + (0.059 * 1000),
    12271.147 + (0.138 * 10000),
    75158.243 + (0.3 * 25000)
  )
)

df <- rbind(samples_df, genes_df, sweep_df)

df %>%
  dplyr::mutate(num_genes = as.factor(num_genes)) %>%
  ggplot(aes(x = num_samples, y = run_time, colour = num_genes)) +
  geom_point() +
  geom_smooth(method = "loess", formula = y ~ x) +
  theme_classic() +
  labs(x = "number of samples",
       y = "run time (s)",
       title = "KNN Impute (K=10)",
       colour = "number of genes")

