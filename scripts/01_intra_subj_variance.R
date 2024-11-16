# Intra-subject variability

suppressPackageStartupMessages({
  library("tidyverse")
  library("here")
  library("rio")
  library("multilevelTools")
  library("lmerTest")
  library("JWileymisc") # testDistribution()
  library("extraoperators") # %!in%
  library("sjPlot") # plot_model()
  library("brms")
  library("cmdstanr")
  library("MplusAutomation")
  library("gt")
  library("glue")
  library("kableExtra")
  library("misty")
  library("missRanger")
})

d <- readRDS(here::here("data", "prep", "ema", "data_for_manuscript.RDS")) |> 
  dplyr::rename(
    sc = psc,
    usc = nsc,
    zsc = zpsc,
    zusc = znsc
  )

cor(d$sc, d$usc)


# Compute intra-subject standard deviation for each variable across all measurement occasions
intra_subject_sd <- d |> 
  group_by(user_id) |> 
  summarize(
    sd_sc = sd(sc, na.rm = TRUE),
    sd_usc = sd(usc, na.rm = TRUE),
    sd_neg_aff = sd(neg_aff, na.rm = TRUE),
    sd_context = sd(context, na.rm = TRUE)
  )

# View the result
print(intra_subject_sd)

plot(intra_subject_sd$sd_sc, intra_subject_sd$usc)

plot(intra_subject_sd$sd_neg_aff, intra_subject_sd$sd_sc)
plot(intra_subject_sd$sd_neg_aff, intra_subject_sd$sd_usc)

plot(intra_subject_sd$sd_context, intra_subject_sd$sd_sc)
plot(intra_subject_sd$sd_context, intra_subject_sd$sd_usc)

ders <- rio::import(
  here::here("data", "prep", "quest", "ders_scores.csv")
) |> 
  dplyr::select(user_id, ders_ts)


tot_df <- left_join(intra_subject_sd, ders, by = "user_id")

d1 <- left_join(d, ders, by = "user_id")


fm1 <- lm(scale(sd_sc) ~ scale(ders_ts), data = tot_df)
summary(fm1)

fm2 <- lm(scale(sd_usc) ~ scale(ders_ts), data = tot_df)
summary(fm2)

# Increased variability in SC is associated with high DERS scores (difficulti in 
# emotion regulation)
# There is no relation between USC variability and DERS scores.


# Location Scale Model ----------------------------------------------------

# options(mc.cores = parallel::detectCores())
# 
# d1$psc <- scale(d1$psc) |> as.numeric()
# d1$nsc <- scale(d1$nsc) |> as.numeric()
# d1$ders <- scale(d1$ders_ts) |> as.numeric()
# 
# 
# mod1_ls <- brm(
#   bf(
#     psc ~ 1 + ders + (1 + ders | user_id / bysubj_day),  
#     sigma ~ 1 + ders + (1 | user_id)
#   ),
#   family = gaussian(),  # Use function for clarity
#   data = d1,
#   seed = 1234,
# 
#   backend = "cmdstanr",
#   silent = 2, refresh = 0
#   # algorithm = "meanfield"
# )
# pp_check(mod1_ls)
# summary(mod1_ls)
# bayes_R2(mod1_ls)


# Save data for Mplus -----------------------------------------------------

# Glimpse the original data
glimpse(d1)

# Prepare the data for Mplus
ema_data <- d1 %>%
  # Select relevant variables
  dplyr::select(
    user_id, 
    bysubj_day, 
    time_window, 
    zsc,       # SC (Self-Compassion)
    zusc,      # USC (Uncompassionate Self-Compassion)
    neg_aff,   # Negative Affect (momentary)
    ders       # DERS (criterion variable, measured once per participant)
  ) %>%
  
  # Rename columns for Mplus compatibility
  dplyr::rename(
    SC_or_USC = zsc,       # Use SC for this analysis (change to zusc for USC analysis)
    NAF = neg_aff,          # Negative Affect
    DERS = ders            # Difficulties in Emotion Regulation Scale
  ) %>%
  
  # Ensure no duplicated rows (if applicable)
  distinct()

# Check for within-cluster variation in DERS
ema_data %>%
  group_by(user_id) %>%
  summarise(DERS_var = var(DERS)) %>%
  dplyr::filter(DERS_var > 0)  # Should return no rows, as DERS is constant within participants

# Aggregate DERS to be constant for each user_id
ema_data <- ema_data %>%
  group_by(user_id) %>%
  mutate(DERS = mean(DERS, na.rm = TRUE)) %>%  # Ensures DERS is constant
  ungroup()

# Verify that DERS is constant within clusters
ema_data %>%
  group_by(user_id) %>%
  summarise(DERS_var = var(DERS, na.rm = TRUE)) %>%
  dplyr::filter(DERS_var > 0)  # Should produce no rows

# Summarize and structure data for Mplus
ema_data <- ema_data %>%
  group_by(user_id) %>%
  summarise(
    bysubj_day = first(bysubj_day),    # Retain representative value
    time_window = first(time_window), # Retain representative value
    SC_or_USC = list(SC_or_USC),      # Retain within-level variable
    NAF = list(NAF),                    # Retain within-level Negative Affect
    DERS = mean(DERS, na.rm = TRUE)   # Ensure DERS is constant
  ) %>%
  unnest(cols = c(SC_or_USC, NAF))  # Restore within-level structure

# Convert user_id to numeric for compatibility with Mplus
ema_data <- ema_data |> 
  mutate(user_id_numeric = as.numeric(factor(user_id))) %>%
  dplyr::select(user_id_numeric, bysubj_day, time_window, SC_or_USC, NAF, DERS) |> 
  dplyr::rename(
    user_num = user_id_numeric,
    day_win = bysubj_day,
    time_win = time_window,
    CMP = SC_or_USC,
    NEGAFF = NAF
  )

# Save data to a .dat file for Mplus
write.table(
  ema_data,
  here::here("scripts", "mplus", "ema_data.dat"),
  row.names = FALSE,
  col.names = FALSE,
  sep = " "
)

# Verify saved .dat file by reloading
ema_data_check <- read.table(
  here::here("scripts", "mplus", "ema_data.dat"),
  header = FALSE, 
  sep = " "
)

# Assign column names for clarity
colnames(ema_data_check) <- 
  c("user_num", "day_win", "time_win", "CMP", "NEGAFF", "DERS")

# Check for consistency in the reloaded data
ema_data_check %>%
  group_by(user_num) %>%
  summarise(DERS_var = var(DERS)) %>%
  filter(DERS_var > 0)  # Should return no rows


# Set the_dir for further use when specifying the location of Mplus models.
the_dir <- here::here("scripts",  "mplus")


# Model 1: Unconditional ------------------------------------------------

runModels(paste0(the_dir, "/m1_unconditional.inp"), showOutput = TRUE)
m1_unc <- MplusAutomation::readModels(paste0(the_dir, "/m1_unconditional.out"))
summary(m1_unc)


# Model 1: Conditional ------------------------------------------------

runModels(paste0(the_dir, "/mod1_cond.inp"), showOutput = TRUE)
m1_unc <- MplusAutomation::readModels(paste0(the_dir, "/m1_unconditional.out"))
summary(m1_unc)


runModels(paste0(the_dir, "/m3_conditional.inp"), showOutput = TRUE)
m3_unc <- MplusAutomation::readModels(paste0(the_dir, "/m3_unconditional.out"))
summary(m1_unc)

#### UCS -----------------------

# Prepare the data for Mplus
ema_data <- d1_imp %>%
  # Select relevant variables
  dplyr::select(user_id, bysubj_day, time_window, psc, nsc, ders) %>%
  
  # Rename columns for Mplus compatibility
  rename(
    CurAve_d = nsc,       # Use NSC (negative SC) as the main variable
    PSC = psc,            # Rename PSC for potential future use
    FSMean = ders         # DERS as the Level-2 variable (renamed for clarity)
  ) %>%
  
  # Ensure no duplicated rows (if applicable)
  distinct()

# Check if FSMean varies within user_id clusters
ema_data %>%
  group_by(user_id) %>%
  summarise(FSMean_var = var(FSMean)) %>%
  filter(FSMean_var > 0)

# Aggregate FSMean to be constant for each user_id
ema_data <- ema_data %>%
  group_by(user_id) %>%
  mutate(FSMean = mean(FSMean, na.rm = TRUE)) %>%
  ungroup()

# Verify that FSMean is constant within clusters
ema_data %>%
  group_by(user_id) %>%
  summarise(FSMean_var = var(FSMean, na.rm = TRUE)) %>%
  filter(FSMean_var > 0)  # Should return no rows

# Summarize and structure data for Mplus
ema_data <- ema_data %>%
  group_by(user_id) %>%
  summarise(
    bysubj_day = first(bysubj_day),       # Retain representative value
    time_window = first(time_window),    # Retain representative value
    CurAve_d = list(CurAve_d),           # Retain within-level variable
    FSMean = mean(FSMean, na.rm = TRUE)  # Ensure FSMean is constant
  ) %>%
  unnest(cols = c(CurAve_d))  # Restore within-level structure

# Convert user_id to numeric for compatibility with Mplus
ema_data <- ema_data %>%
  mutate(user_id_numeric = as.numeric(factor(user_id))) %>%
  dplyr::select(user_id_numeric, bysubj_day, time_window, CurAve_d, FSMean)

# Save data to a .dat file for Mplus
write.table(
  ema_data,
  here::here("workflows", "scripts", "ema", "mplus_models", "ema_data.dat"),
  row.names = FALSE,
  col.names = FALSE,
  sep = " "
)


runModels(paste0(the_dir, "/m3_conditional.inp"), showOutput = TRUE)
m3_unc <- MplusAutomation::readModels(paste0(the_dir, "/m3_unconditional.out"))
summary(m1_unc)
