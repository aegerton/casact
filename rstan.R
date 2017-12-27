library(tidyverse)
library(actuar)
library(rstan)
library(scales)
library(bayesplot)
library(cowplot)
### File was downloaded from http://www.casact.org/research/reserve_data/ppauto_pos.csv
options(width = 80L
        ,warn  = 1
        ,mc.cores = parallel::detectCores()
)
rstan_options(auto_write = TRUE)
theme_set(theme_cowplot())
set.seed(42)
stan_seed <- 42
source("custom_functions.R")
data_files <- dir(pattern = "ppauto_pos.csv", full.names = TRUE)

data_cols <- cols(GRCODE = col_character())

rawdata_tbl <- data_files %>%
  map(read_claim_datafile, col_type = data_cols) %>%
  bind_rows
rawdata_tbl <- as.tibble(read.csv("ppauto_pos.csv"))
glimpse(rawdata_tbl)
claimdata_tbl <- rawdata_tbl %>%
  mutate(acc_year   = as.character(AccidentYear)
         ,dev_year   = DevelopmentYear
         ,dev_lag    = DevelopmentLag
         ,premium    = EarnedPremDIR_B
         ,cum_loss   = CumPaidLoss_B
         ,loss_ratio = cum_loss / premium) %>%
  select(GRCODE, GRNAME, acc_year, dev_year, dev_lag, premium, cum_loss, loss_ratio)
use_grcode <- c(43,353,388,620)
carrier_snapshot_tbl <- claimdata_tbl  %>%
  filter(GRCODE %in% use_grcode
         ,dev_year < 1998)
ggplot(carrier_snapshot_tbl) +
  geom_line(aes(x = dev_lag, y = loss_ratio, colour = as.character(acc_year))
            ,size = 0.3) +
  expand_limits(y = c(0,1)) +
  facet_wrap(~GRCODE) +
  xlab('Development Time') +
  ylab('Loss Ratio') +
  ggtitle('Snapshot of Loss Curves for 10 Years of Loss Development'
          ,subtitle = 'Private Passenger Auto Insurance for Single Organisation') +
  guides(colour = guide_legend(title = 'Cohort Year'))
snapshot_tbl <- carrier_snapshot_tbl %>%
  filter(GRCODE %in% use_grcode[1])
snapshot_tbl %>%
  select(acc_year, dev_lag, premium, cum_loss) %>%
  spread(dev_lag, cum_loss) %>%
  print
snapshot_tbl %>%
  select(acc_year, dev_lag, premium, loss_ratio) %>%
  spread(dev_lag, loss_ratio) %>%
  print.data.frame(digits = 2)
ggplot(snapshot_tbl) +
  geom_line(aes(x = dev_lag, y = loss_ratio, colour = acc_year)
            ,size = 0.3) +
  expand_limits(y = 0) +
  xlab('Development Time') +
  ylab('Loss Ratio') +
  ggtitle("Loss Ratio Curves by Development Time") +
  guides(colour = guide_legend(title = 'Cohort Year'))
t_seq <- seq(0, 15, by = 0.01)
loglogistic_func <- function(t, om, th) 1 - exp(-(t/th)^om)
weibull_func     <- function(t, om, th) t^om / (t^om + th^om)
weibull_tbl <- tibble(
  label = 'Weibull'
  ,t = t_seq
  ,value = weibull_func(t_seq, 1.5, 2.2)
)
loglogistic_tbl <- tibble(
  label = 'Log-logistic'
  ,t = t_seq
  ,value = loglogistic_func(t_seq, 1.5, 2.2)
)
plot_tbl <- bind_rows(weibull_tbl, loglogistic_tbl)
ggplot(plot_tbl) +
  geom_line(aes(x = t, y = value, colour = label)) +
  xlab(expression(t)) +
  ylab(expression("Growth Factor for (" * omega * "=1.5, " * theta * "=2.2)")) +
  ggtitle("Sample Curves for Log-Logistic and Weibull Forms")
ll_vals <- loglogistic_tbl$value
new_param_func <- function(x) {
  omega <- x[1]
  theta <- x[2]
  
  new_vals <- weibull_func(t_seq, omega, theta)
  
  tot_ss <- sum((new_vals - ll_vals)^2)
  
  return(tot_ss)
}
optim_params <- optim(c(1, 1), new_param_func)
fittedweibull_tbl <- tibble(
  label = 'Weibull (fitted)'
  ,t = t_seq
  ,value = weibull_func(t_seq
                        ,optim_params$par[1]
                        ,optim_params$par[2])
)
plot_tbl <- bind_rows(weibull_tbl
                      ,loglogistic_tbl
                      ,fittedweibull_tbl)
ggplot(plot_tbl) +
  geom_line(aes(x = t, y = value, colour = label)) +
  xlab(expression(t)) +
  ylab(expression("Functional Forms for Growth/Development Factors")) +
  ggtitle("Comparison Plot for Weibull and Log-Logistic Curves")   
modeldata_tbl <- claimdata_tbl %>%
  filter(GRCODE == use_grcode[1])

usedata_tbl <- modeldata_tbl %>%
  filter(dev_year < 1998)

cohort_maxtime <- usedata_tbl %>%
  group_by(acc_year) %>%
  summarise(maxtime = max(dev_lag)) %>%
  arrange(acc_year) %>%
  pull(maxtime)

cohort_premium <- usedata_tbl %>%
  group_by(acc_year) %>%
  summarise(premium = unique(premium)) %>%
  pull(premium)

t_values <- usedata_tbl %>%
  select(dev_lag) %>%
  arrange(dev_lag) %>%
  unique %>%
  pull(dev_lag)

standata_lst <- list(
  growthmodel_id = 1   # Use weibull rather than loglogistic
  ,n_data         = usedata_tbl %>% nrow
  ,n_time         = usedata_tbl %>% select(dev_lag)  %>% unique %>% nrow
  ,n_cohort       = usedata_tbl %>% select(acc_year) %>% unique %>% nrow
  ,cohort_id      = get_character_index(usedata_tbl$acc_year)
  ,cohort_maxtime = cohort_maxtime
  ,t_value        = t_values
  ,t_idx          = get_character_index(usedata_tbl$dev_lag)
  ,premium        = cohort_premium
  ,loss           = usedata_tbl$cum_loss
)
stan_file <- "losscurves_sislob.stan"
cat(read_lines(stan_file), sep = "\n")
model_sislob_stanmodel <- stan_model(stan_file)
model_sislob_stanfit <- sampling(
  object = model_sislob_stanmodel
  ,data   = standata_lst
  ,iter   = 500
  ,chains = 8
  ,seed   = stan_seed
)