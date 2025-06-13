targets::tar_load(c(model_arm_lp_data,model_thigh_lp_data, ranef_cors))

variances_arm <- VarCorr(model_arm_lp_data)
variances_thigh <- VarCorr(model_thigh_lp_data)

b0 = 0
b_time = 0.05
b_cond_time = 0.025
sd_site_arm = variances_arm$site_id[1]
sd_participant_arm = variances_arm$participant_id[1]
sigma_arm = sigma(model_arm_lp_data)
sd_site_thigh = variances_thigh$site_id[1]
sd_participant_thigh = variances_thigh$participant_id[1]
sigma_thigh = sigma(model_thigh_lp_data)
u_participant_cor = ranef_cors$cor[1]
u_site_cor = ranef_cors$cor[2]


# b0 = 0, b_time = 0.05, b_cond_time = 0.025
# sd_site = variances_model$site_id[1]
# sd_participant = variances_model$participant_id[1]
# sigma = sigma(model)

data <- add_random(site_id = 24) |>
  add_random(participant_id = 100, .nested_in = "site_id") |>
  add_between("participant_id", cond = c("lower_volume", "higher_volume"), .shuffle = TRUE) |>
  add_recode("cond", "cond_dummy", lower_volume = 0, higher_volume = 1) |>
  add_within("participant_id", 
             time = c(0,1), 
             measurement = c(1,2,3),
             outcome = c("arm", "thigh")) |>
  add_recode("outcome", "outcome_dummy", arm = -0.5, thigh = 0.5) |>
  add_ranef("site_id", u_site_arm = sd_site_arm, u_site_thigh = sd_site_thigh, .cors = u_site_cor) |>
  add_ranef("participant_id", u_participant_arm = sd_participant_arm, u_participant_thigh = sd_participant_thigh, .cors = u_participant_cor) |>
  add_ranef(error_arm = sigma_arm,
            error_thigh = sigma_thigh) |>
  # rescale the standardised effects to the current variance estimates
  mutate(b_time_arm = b_time * (sd_participant_arm^2 + sd_site_arm^2 + error_arm^2), 
         b_cond_time_arm = b_cond_time * (sd_participant_arm^2 + sd_site_arm^2 + error_arm^2),
         b_time_thigh = b_time * (sd_participant_thigh^2 + sd_site_thigh^2 + error_thigh^2), 
         b_cond_time_thigh = b_cond_time * (sd_participant_thigh^2 + sd_site_thigh^2 + error_thigh^2)) |>
  # calculate outcomes based on assumed model
  mutate(ma_z = case_when(
    outcome == "thigh" ~ (b0 + u_site_thigh + u_participant_thigh) + (b_time_thigh * time) + (b_cond_time_thigh * cond_dummy * time) + error_thigh,
    outcome == "arm" ~ (b0 + u_site_arm + u_participant_arm) + (b_time_arm * time) + (b_cond_time_arm * cond_dummy * time) + error_arm,
  ))|>
  
  select(site_id, participant_id, cond_dummy, outcome_dummy, time, measurement, ma_z)  

# simulate loss to follow
follow_up_participants <- data |>
  filter(time == 1) %>%
  distinct(participant_id) %>%
  slice_sample(prop = 1-0.15)

# Filter the original data
data <- data %>%
  filter(
    time == 0 | 
      (time == 1 & participant_id %in% follow_up_participants$participant_id)
  )

data_arm <- data |>
  filter(
    outcome_dummy == -0.5
  )

data_thigh <- data |>
  filter(
    outcome_dummy == 0.5
  )


model_all <- lmer(ma_z ~ time + time:cond_dummy + (1|site_id) + (1 | participant_id),
                  data = data, REML = TRUE)

model_all_outcome_fixed <- lmer(ma_z ~ outcome_dummy + time + time:outcome_dummy + time:cond_dummy + time:cond_dummy:outcome_dummy + 
                                  (1|site_id) + (1 | participant_id),
                                data = data, REML = TRUE)

model_all_outcome_random <- lmer(ma_z ~ outcome_dummy + time + time:outcome_dummy + time:cond_dummy + time:cond_dummy:outcome_dummy + 
                                   (1|site_id/outcome_dummy) + (1 | participant_id/outcome_dummy),
                                 data = data, REML = TRUE)

model_all_outcome_disp <- glmmTMB(ma_z ~ time + time:cond_dummy + (1|site_id) + (1 | participant_id),
                                  dispformula = ~ outcome_dummy,
                                  data = data, REML = TRUE)

model_upper <- lmer(ma_z ~ time + time:cond_dummy + (1|site_id) + (1 | participant_id),
                    data = data_arm, REML = TRUE)

model_lower <- lmer(ma_z ~ time + time:cond_dummy + (1|site_id) + (1 | participant_id),
                    data = data_thigh, REML = TRUE)

coef_upper <- tidy(model_upper, effects = "fixed") %>%
  filter(term == "time:cond_dummy") %>%
  select(estimate, std.error)

coef_lower <- tidy(model_lower, effects = "fixed") %>%
  filter(term == "time:cond_dummy") %>%
  select(estimate, std.error)

meta_input <- bind_rows(
  mutate(coef_upper, study = "upper"),
  mutate(coef_lower, study = "lower")
) %>%
  rename(
    yi = estimate,   # Effect size
    sei = std.error  # Standard error
  )

meta_result_re <- rma(yi = yi, sei = sei, data = meta_input, method = "REML")
meta_result_fe <- rma(yi = yi, sei = sei, data = meta_input, method = "FE")

meta_tests <- bind_rows(
  tibble(
    estimate = coef(meta_result_re),
    se = meta_result_re$se,
    p.value = meta_result_re$pval,
    model = "meta_re"
  ),
  tibble(
    estimate = coef(meta_result_fe),
    se = meta_result_fe$se,
    p.value = meta_result_fe$pval,
    model = "meta_fe"
  )
) |>
  group_by(model) |>
  mutate(
    z_eq_lo = (estimate - (-0.1)) / se,
    z_eq_hi = (estimate - 0.1) / se,
    p_lo = 1 - pnorm(z_eq_lo),
    p_hi = pnorm(z_eq_hi),
    p.value.equiv = max(c(p_lo,p_hi)),
    term = "time:cond_dummy"
  ) |>
  select(model, term, p.value, p.value.equiv)


tests <- bind_rows(
  hypotheses(model_all, equivalence = c(-0.1,0.1), df = insight::get_df(model_all)) |> mutate(model = "pooled"),
  hypotheses(model_all_outcome_fixed, equivalence = c(-0.1,0.1), df = insight::get_df(model_all_outcome_fixed)) |> mutate(model = "pooled_fixed"),
  hypotheses(model_all_outcome_random, equivalence = c(-0.1,0.1), df = insight::get_df(model_all_outcome_random)) |> mutate(model = "pooled_random"),
  hypotheses(model_all_outcome_disp, equivalence = c(-0.1,0.1), df = insight::get_df(model_all_outcome_disp)) |> mutate(model = "pooled_disp"),
  hypotheses(model_upper, equivalence = c(-0.1,0.1), df = insight::get_df(model_upper)) |> mutate(model = "upper"),
  hypotheses(model_lower, equivalence = c(-0.1,0.1), df = insight::get_df(model_lower)) |> mutate(model = "lower")
) |>
  mutate(term = case_when(
    term == "conditional_time:cond_dummy" ~ "time:cond_dummy",
    .default = term
  )) |>
  filter(term == "time:cond_dummy") |>
  select(model, term, p.value, p.value.equiv)

# combine models p values
all_tests <- bind_rows(meta_tests, tests)

# dat <- add_random(site_id = 15) |>
#   add_random(participant_id = 30, .nested_in = "site_id") |>
#   add_between("participant_id", order = c("lower_upper", "upper_lower")) |>
#   add_within("participant_id", study = c(1,2)) |>
#   add_between("participant_id", cond_1 = c("lower_volume", "higher_volume"), .shuffle = TRUE) |>
#   add_between("participant_id", cond_2 = c("lower_volume", "higher_volume"), .shuffle = TRUE) |>
#   mutate(
#     cond = case_when(
#       study == 1 ~ cond_1,
#       study == 2 ~ cond_2
#     )
#   ) |>
#   select(-cond_1, -cond_2) |>
#   add_recode("cond", "cond_dummy", lower_volume = 0, higher_volume = 1) |>
#   add_within("participant_id", time = c(0,1), 
#              measurement = c(1,2,3)) |>
#   add_ranef("site_id", u_site_arm = sd_site_arm, u_site_thigh = sd_site_thigh, .cors = u_site_cor) |>
#   add_ranef("participant_id", u_participant_arm = sd_participant_arm, u_participant_thigh = sd_participant_thigh, .cors = u_participant_cor) |>
#   add_ranef(error_arm = sigma_arm,
#             error_thigh = sigma_thigh) |>
#   # rescale the standardised effects to the current variance estimates
#   mutate(b_time_arm = b_time * sqrt(sd_participant_arm^2 + sd_site_arm^2 + sigma_arm^2), 
#          b_cond_time_arm = b_cond_time * sqrt(sd_participant_arm^2 + sd_site_arm^2 + sigma_arm^2),
#          b_time_thigh = b_time * sqrt(sd_participant_thigh^2 + sd_site_thigh^2 + sigma_thigh^2), 
#          b_cond_time_thigh = b_cond_time * sqrt(sd_participant_thigh^2 + sd_site_thigh^2 + sigma_thigh^2)) |>
#   # calculate outcomes based on assumed model
#   mutate(ma_z = case_when(
#     order == "lower_upper" & study == 1 ~ (b0 + u_site_thigh + u_participant_thigh) + (b_time_thigh * time) + (b_cond_time_thigh * cond_dummy * time) + error_thigh,
#     order == "lower_upper" & study == 2 ~ (b0 + u_site_arm + u_participant_arm) + (b_time_arm * time) + (b_cond_time_arm * cond_dummy * time) + error_arm,
#     order == "upper_lower" & study == 2 ~ (b0 + u_site_thigh + u_participant_thigh) + (b_time_thigh * time) + (b_cond_time_thigh * cond_dummy * time) + error_thigh,
#     order == "upper_lower" & study == 1 ~ (b0 + u_site_arm + u_participant_arm) + (b_time_arm * time) + (b_cond_time_arm * cond_dummy * time) + error_arm
#     
#     ))|>
#   
#   # add outcome code
#   mutate(outcome_dummy = case_when(
#     (order == "lower_upper" & study == 2) | (order == "upper_lower" & study == 1) ~ -0.5, # arm
#     (order == "lower_upper" & study == 1) | (order == "upper_lower" & study == 2) ~ 0.5 # thigh
#   )) |>
#   
#   select(site_id, participant_id, order, study, cond_dummy, outcome_dummy, time, measurement, ma_z)
#   
#   # simulate loss to follow
#   study_2_participants <- dat |>
#     filter(study == 2) %>%
#     distinct(participant_id) %>%
#     slice_sample(prop = 1-0.15)
#   
#   # Filter the original data
#   dat <- dat %>%
#     filter(
#       study == 1 | 
#         (study == 2 & participant_id %in% study_2_participants$participant_id)
#     )
# 
#   
#   
#   dat_arm <- dat |>
#     filter(
#       (order == "lower_upper" & study == 2) | (order == "upper_lower" & study == 1)
#     )
#   
#   dat_thigh <- dat |>
#     filter(
#       (order == "lower_upper" & study == 1) | (order == "upper_lower" & study == 2)
#     )
#   
#   
#   model_all_outcome_fixed <- lmer(ma_z ~ outcome_dummy + time + time:outcome_dummy + time:cond_dummy + time:cond_dummy:outcome_dummy + 
#                                     (1|site_id) + (1 | participant_id),
#                                   data = dat, REML = TRUE)
#   
#   model_all_outcome_random <- lmer(ma_z ~ outcome_dummy + time + time:outcome_dummy + time:cond_dummy + time:cond_dummy:outcome_dummy + 
#                                      (1|site_id/outcome_dummy) + (1 | participant_id/outcome_dummy),
#                                    data = dat, REML = TRUE)
#   
#   model_all_outcome_disp <- glmmTMB(ma_z ~ time + time:cond_dummy + (1|site_id) + (1 | participant_id),
#                                     dispformula = ~ outcome_dummy,
#                                     data = dat, REML = TRUE)
#   
#   model_upper <- lmer(ma_z ~ time + time:cond_dummy + (1|site_id) + (1 | participant_id),
#                       data = dat_arm, REML = TRUE)
#   
#   model_lower <- lmer(ma_z ~ time + time:cond_dummy + (1|site_id) + (1 | participant_id),
#                       data = dat_thigh, REML = TRUE)
#   
#   
#   
#   
#   t <- hypotheses(model_all_outcome_disp, equivalence = c(-0.1,0.1), df = insight::get_df(model_all_outcome_disp)) |> mutate(model = "pooled_disp")
#   
#   
#   
#   # lp study
#   
#   # set up data structure
#   dat <- add_random(site = 15) |>
#     add_random(participant = 20, .nested_in = "site") |>
#     # add_between("participant", cond = c("full_ROM", "lengthened_partial")) |>
#     # add_recode("cond", "cond_dummy", full_ROM = 0, lengthened_partial = 1) |>
#     
#     # add study replicate
#     add_between("participant", order = c("lower_upper", "upper_lower")) |>
#     add_within("participant", study = c(1,2)) |>
#     add_between("participant", cond_1 = c("lower_volume", "higher_volume"), .shuffle = TRUE) |>
#     add_between("participant", cond_2 = c("lower_volume", "higher_volume"), .shuffle = TRUE) |>
#     mutate(
#       cond = case_when(
#         study == 1 ~ cond_1,
#         study == 2 ~ cond_2
#       )
#     ) |>
#     select(-cond_1, -cond_2) |>
#     add_recode("cond", "cond_dummy", lower_volume = 0, higher_volume = 1) |>
#     
#     add_within("participant", time = seq(0,1, by = 1), 
#                measurement = seq(1:3)) |>
#     add_ranef("site", u_site = 0.205) |>
#     add_ranef("participant", u_participant_arm = 0.629) |>
#     add_ranef("participant", u_participant_thigh = 0.597) |>
#     add_ranef(arm_ma_error = 0.166) |>
#     add_ranef(thigh_ma_error = 0.198) |>
#     mutate(arm_ma = case_when(
#       order == "lower_upper" & study == 2 ~ (b0 + u_site + u_participant_arm) + (b_time * time) + 0 + (b_cond_time * cond_dummy * time) + arm_ma_error,
#       order == "upper_lower" & study == 1 ~ (b0 + u_site + u_participant_arm) + (b_time * time) + 0 + (b_cond_time * cond_dummy * time) + arm_ma_error,
#       .default = NA
#     ),
#     thigh_ma = case_when(
#       order == "lower_upper" & study == 1 ~ (b0 + u_site + u_participant_thigh) + (b_time * time) + 0 + (b_cond_time * cond_dummy * time) + thigh_ma_error,
#       order == "upper_lower" & study == 2 ~ (b0 + u_site + u_participant_thigh) + (b_time * time) + 0 + (b_cond_time * cond_dummy * time) + thigh_ma_error,
#       .default = NA
#       )
#     )
#   
#   dat_long <- dat |>
#     select(participant, site, time, cond_dummy, arm_ma, thigh_ma) |>
#     pivot_longer(5:6,
#                  names_to = "muscle",
#                  values_to = "ma_z") |>
#     filter(!is.na(ma_z))
#   
#   model_all <- lmer(ma_z ~ time + time:cond_dummy + (1|site) + (1 | participant),
#                     data = dat_long, REML = TRUE)
#   
#   
#   # run mixed effect model and return relevant values
#   model_arm <- lmer(arm_ma ~ time + time:cond_dummy + (1|site) + (1 | participant),
#                     data = dat, REML = TRUE)
#   
#   model_thigh <- lmer(thigh_ma ~ time + time:cond_dummy + (1|site) + (1 | participant),
#                       data = dat, REML = TRUE)
#   
#   
#   # if (time_n == 1 && measurements_n == 1) {
#   #   # run mixed effect model and return relevant values
#   #   model_arm <- lmer(arm_ma ~ time*cond_dummy + (1|site) + (1 | participant),
#   #                     data = dat, REML = TRUE)
#   #   
#   #   model_thigh <- lmer(thigh_ma ~ time*cond_dummy + (1|site) + (1 | participant),
#   #                       data = dat, REML = TRUE)
#   # } else {
#   #   # run mixed effect model and return relevant values
#   #   model_arm <- lmer(arm_ma ~ time*cond_dummy + (1|site) + (time | participant),
#   #                     data = dat, REML = TRUE)
#   #   
#   #   model_thigh <- lmer(thigh_ma ~ time*cond_dummy + (1|site) + (time | participant),
#   #                       data = dat, REML = TRUE)
#   # }
#   
#   hypothesis_tests_arm <- hypotheses(model_arm, equivalence = c(NA,0.1), 
#                                      # vcov = "satterthwaite", 
#                                      df = insight::get_df(model_arm))[3,] |>
#     mutate(muscle_site = "arm")
#   
#   
#   hypothesis_tests_thigh <- hypotheses(model_thigh, equivalence = c(NA,0.1), 
#                                        # vcov = "satterthwaite", 
#                                        df = insight::get_df(model_thigh))[3,] |>
#     mutate(muscle_site = "thigh")
#   
  