##### Functions for pre-registration and sample size estimation ----
prep_lp_data <- function(data) {
  data <- data |>
    filter(time != "1") |>
    mutate(
      thigh_ma_1 =(4.68*thigh_circum_1)-(2.09*thigh_skinfold_1)-80.99,
      thigh_ma_2 =(4.68*thigh_circum_2)-(2.09*thigh_skinfold_2)-80.99,
      thigh_ma_3 =(4.68*thigh_circum_3)-(2.09*thigh_skinfold_3)-80.99,
      arm_ma_1 = if_else(sex == "M" ,((arm_circum_1 - (pi * (arm_skinfold_1/10)))^2 / (4*pi))-10, ((arm_circum_1 - (pi * (arm_skinfold_1/10)))^2 / (4*pi))-6.5),
      arm_ma_2 = if_else(sex == "M" ,((arm_circum_2 - (pi * (arm_skinfold_2/10)))^2 / (4*pi))-10, ((arm_circum_2 - (pi * (arm_skinfold_2/10)))^2 / (4*pi))-6.5),
      arm_ma_3 = if_else(sex == "M" ,((arm_circum_3 - (pi * (arm_skinfold_3/10)))^2 / (4*pi))-10, ((arm_circum_3 - (pi * (arm_skinfold_3/10)))^2 / (4*pi))-6.5)
    ) |>
    select(site_id, participant_id, condition, time, contains("ma")) |>
    pivot_longer(5:10,
                 names_to = "muscle",
                 values_to = "estimated_ma") |>
    mutate(
      measurement = case_when(
        str_detect(muscle, "1") == TRUE ~ 1,
        str_detect(muscle, "2") == TRUE ~ 2,
        str_detect(muscle, "3") == TRUE ~ 3
      ),
      muscle = case_when(
        str_detect(muscle, "thigh") == TRUE ~ "thigh",
        str_detect(muscle, "arm") == TRUE ~ "arm"
      )
    ) |>
    mutate(cond_dummy = case_when(
      condition == "fROM" ~ -0.5,
      condition == "lpROM" ~ 0.5
    ),
    timepoint = factor(case_when(
      time == 2 ~ "Pre",
      time == 3 ~ "Post"
    ), levels = c("Pre", "Post"))
    ) |>
    select(site_id, participant_id, condition, cond_dummy, time, timepoint, muscle, estimated_ma)
  
  # calculate z-score numerator
  std_lm_arm <- lm(estimated_ma ~ cond_dummy,
                   data = data |>
                     filter(muscle == "arm"))
  
  std_arm <- summary(std_lm_arm)$sigma
  
  
  std_lm_thigh <- lm(estimated_ma ~ cond_dummy,
                     data = data |>
                       filter(muscle == "thigh"))
  
  std_thigh <- summary(std_lm_thigh)$sigma
  
  
  
  data <- data |>
    group_by(muscle) |>
    mutate(
      estimated_ma_z = case_when(
        muscle == "arm" ~ (estimated_ma - mean(estimated_ma, na.rm=TRUE)) / std_arm,
        muscle == "thigh" ~ (estimated_ma - mean(estimated_ma, na.rm=TRUE)) / std_thigh
        
      ) 
    )
  
  
}

fit_model_lp_data <- function(data) {
  lmer(estimated_ma_z ~ time + time:cond_dummy + (1|site_id) + (1 | participant_id),
       data = data, REML = TRUE)
}

# Simulating from a DGP of parameters from the analysis of the LP study data pooled across arm/thigh muscles
generate_simulation_grid_pooled <- function(model) {
  
  variances_model <- VarCorr(model)

  crossing(
    rep = 1:100,
    participant_n = seq(5, 100, by = 5),
    site_n = 15,
    measurements_n = 3,
    dropout_prop = c(0,0.15)
  ) %>%
    mutate(
      b0 = 0, b_time = 0.05, b_cond_time = 0.025,
      sd_site = variances_model$site_id[1],
      sd_participant = variances_model$participant_id[1],
      sigma = sigma(model)
    )
}

sim_pooled <- function(params) {
  
  with(params, {
    # set up data structure
    data <- add_random(site_id = site_n) |>
      add_random(participant_id = participant_n, .nested_in = "site_id") |>
      add_between("participant_id", order = c("lower_upper", "upper_lower")) |>
      add_within("participant_id", study = c(-0.5,0.5)) |>
      add_between("participant_id", cond_1 = c("lower_volume", "higher_volume"), .shuffle = TRUE) |>
      add_between("participant_id", cond_2 = c("lower_volume", "higher_volume"), .shuffle = TRUE) |>
      mutate(
        cond = case_when(
          study == -0.5 ~ cond_1,
          study == 0.5 ~ cond_2
        )
      ) |>
      select(-cond_1, -cond_2) |>
      add_recode("cond", "cond_dummy", lower_volume = 0, higher_volume = 1) |>
      add_within("participant_id", time = c(0,1), 
                 measurement = c(1,2,3)) |>
      add_ranef("site_id", u_site = sd_site) |>
      add_ranef("participant_id", u_participant = sd_participant) |>
      add_ranef(error = sigma) |>
      # rescale the standardised effects to the current variance estimates
      mutate(b_time_rescale = b_time * (sd_participant^2 + sd_site^2 + error^2), 
             b_cond_time_rescale = b_cond_time * (sd_participant^2 + sd_site^2 + error^2)) |>
      # calculate outcomes based on assumed model
      mutate(ma_z = (b0 + u_site + u_participant) + (b_time_rescale * time) + (b_cond_time_rescale * cond_dummy * time) + error)|>
      
      # add outcome code
      mutate(outcome_dummy = case_when(
        (order == "lower_upper" & study == 0.5) | (order == "upper_lower" & study == -0.5) ~ -0.5, # arm
        (order == "lower_upper" & study == -0.5) | (order == "upper_lower" & study == 0.5) ~ 0.5 # thigh
      )) |>
      
      select(site_id, participant_id, order, study, cond_dummy, outcome_dummy, time, measurement, ma_z) 
    
    # simulate loss to follow
    study_2_participants <- data |>
      filter(study == 0.5) %>%
      distinct(participant_id) %>%
      slice_sample(prop = 1-dropout_prop)
    
    # Filter the original data
    data <- data %>%
      filter(
        study == -0.5 | 
          (study == 0.5 & participant_id %in% study_2_participants$participant_id)
      )
    
    data_arm <- data |>
      filter(
        (order == "lower_upper" & study == 0.5) | (order == "upper_lower" & study == -0.5)
      )
    
    data_thigh <- data |>
      filter(
        (order == "lower_upper" & study == -0.5) | (order == "upper_lower" & study == 0.5)
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
    
    # Add parameter info to each row
    bind_cols(params, all_tests)
  })
  
  
}

# Simulating from a DGP of parameters from separate analysis of the LP study data by each muscle and the correlation of random effects

check_ranef_cors <- function(data, model_arm, model_thigh) {
  
  ranef_arm_participant <- tibble(
    id = rownames(ranef(model_arm)$participant_id),
    arm_ranef = ranef(model_arm)$participant_id
  )
  
  ranef_thigh_participant <- tibble(
    id = rownames(ranef(model_thigh)$participant_id),
    thigh_ranef = ranef(model_thigh)$participant_id
  )
  
  ranefs_participant <- left_join(ranef_arm_participant, ranef_thigh_participant, by = "id") |>
    mutate(param = "Random Intercept: Participant")
  
  ranef_cor_participant <- cor.test(ranefs_participant$arm_ranef$`(Intercept)`, ranefs_participant$thigh_ranef$`(Intercept)`)
  
  ranef_arm_site <- tibble(
    id = rownames(ranef(model_arm)$site_id),
    arm_ranef = ranef(model_arm)$site_id
  )
  
  ranef_thigh_site <- tibble(
    id = rownames(ranef(model_thigh)$site_id),
    thigh_ranef = ranef(model_thigh)$site_id
  )
  
  ranefs_site <- left_join(ranef_arm_site, ranef_thigh_site, by = "id") |>
    mutate(param = "Random Intercept: Site")
  
  ranef_cor_site <- cor.test(ranefs_site$arm_ranef$`(Intercept)`, ranefs_site$thigh_ranef$`(Intercept)`)
  
  
  ranefs <- bind_rows(ranefs_participant, ranefs_site)
  
  ranef_cors <- bind_rows(
    tibble(param = "Random Intercept: Participant",
           cor = ranef_cor_participant$estimate,
           conf.low = ranef_cor_participant$conf.int[1],
           conf.high = ranef_cor_participant$conf.int[2],
           y = max(ranefs_participant$arm_ranef$`(Intercept)`)),
    tibble(param = "Random Intercept: Site",
           cor = ranef_cor_site$estimate,
           conf.low = ranef_cor_site$conf.int[1],
           conf.high = ranef_cor_site$conf.int[2],
           y = max(ranefs_site$arm_ranef$`(Intercept)`)),
  )
  
  return(ranef_cors)
}

plot_check_ranef_cors <- function(data, model_arm, model_thigh) {
  
  ranef_arm_participant <- tibble(
    id = rownames(ranef(model_arm)$participant_id),
    arm_ranef = ranef(model_arm)$participant_id
  )
  
  ranef_thigh_participant <- tibble(
    id = rownames(ranef(model_thigh)$participant_id),
    thigh_ranef = ranef(model_thigh)$participant_id
  )
  
  ranefs_participant <- left_join(ranef_arm_participant, ranef_thigh_participant, by = "id") |>
    mutate(param = "Random Intercept: Participant")
  
  ranef_cor_participant <- cor.test(ranefs_participant$arm_ranef$`(Intercept)`, ranefs_participant$thigh_ranef$`(Intercept)`)
  
  ranef_arm_site <- tibble(
    id = rownames(ranef(model_arm)$site_id),
    arm_ranef = ranef(model_arm)$site_id
  )
  
  ranef_thigh_site <- tibble(
    id = rownames(ranef(model_thigh)$site_id),
    thigh_ranef = ranef(model_thigh)$site_id
  )
  
  ranefs_site <- left_join(ranef_arm_site, ranef_thigh_site, by = "id") |>
    mutate(param = "Random Intercept: Site")
  
  ranef_cor_site <- cor.test(ranefs_site$arm_ranef$`(Intercept)`, ranefs_site$thigh_ranef$`(Intercept)`)
  
  
  ranefs <- bind_rows(ranefs_participant, ranefs_site)
  
  ranef_cors <- bind_rows(
    tibble(param = "Random Intercept: Participant",
           cor = ranef_cor_participant$estimate,
           conf.low = ranef_cor_participant$conf.int[1],
           conf.high = ranef_cor_participant$conf.int[2],
           y = max(ranefs_participant$arm_ranef$`(Intercept)`)),
    tibble(param = "Random Intercept: Site",
           cor = ranef_cor_site$estimate,
           conf.low = ranef_cor_site$conf.int[1],
           conf.high = ranef_cor_site$conf.int[2],
           y = max(ranefs_site$arm_ranef$`(Intercept)`)),
  )
  
  
  ranefs |>
    ggplot(aes(x=thigh_ranef$`(Intercept)`, y=arm_ranef$`(Intercept)`)) +
    geom_point() +
    geom_smooth(method = "lm") +
    geom_text(aes(x=0, y = y*1.2, 
                  label = glue::glue("Correlation: {round(cor, 2)} [95% CI: {round(conf.low, 2)},  {round(conf.high, 2)}]")),
              size = 3,
              data = ranef_cors) +
    labs(
      x = "Random effect for thigh muscle area model",
      y = "Random effect for arm muscle area model",
      caption = "Using data from DOI: 10.51224/SRXIV.485\nFitting model estimated_ma_z ~ time + time:cond_dummy + (1|site_id) + (1 | participant_id) to each muscle separately"
    ) +
    ggh4x::facet_wrap2("param", scales = "free")
}

generate_simulation_grid_separate <- function(model_arm, model_thigh, ranef_cors) {
  
  variances_arm <- VarCorr(model_arm)
  variances_thigh <- VarCorr(model_thigh)
  
  
  crossing(
    rep = 1:100,
    participant_n = seq(5, 100, by = 5),
    site_n = 15,
    measurements_n = 3,
    dropout_prop = c(0,0.15)
  ) %>%
    mutate(
      b0 = 0, b_time = 0.05, b_cond_time = 0.025,
      sd_site_arm = variances_arm$site_id[1],
      sd_participant_arm = variances_arm$participant_id[1],
      sigma_arm = sigma(model_arm),
      sd_site_thigh = variances_thigh$site_id[1],
      sd_participant_thigh = variances_thigh$participant_id[1],
      sigma_thigh = sigma(model_thigh),
      u_participant_cor = ranef_cors$cor[1],
      u_site_cor = ranef_cors$cor[2]
      
    )
}

sim_separate <- function(params) {
  
  with(params, {
    # set up data structure
    data <- add_random(site_id = site_n) |>
      add_random(participant_id = participant_n, .nested_in = "site_id") |>
      add_between("participant_id", order = c("lower_upper", "upper_lower")) |>
      add_within("participant_id", study = c(-0.5,0.5)) |>
      add_between("participant_id", cond_1 = c("lower_volume", "higher_volume"), .shuffle = TRUE) |>
      add_between("participant_id", cond_2 = c("lower_volume", "higher_volume"), .shuffle = TRUE) |>
      mutate(
        cond = case_when(
          study == -0.5 ~ cond_1,
          study == 0.5 ~ cond_2
        )
      ) |>
      select(-cond_1, -cond_2) |>
      add_recode("cond", "cond_dummy", lower_volume = 0, higher_volume = 1) |>
      add_within("participant_id", time = c(0,1), 
                 measurement = c(1,2,3)) |>
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
        order == "lower_upper" & study == -0.5 ~ (b0 + u_site_thigh + u_participant_thigh) + (b_time_thigh * time) + (b_cond_time_thigh * cond_dummy * time) + error_thigh,
        order == "lower_upper" & study == 0.5 ~ (b0 + u_site_arm + u_participant_arm) + (b_time_arm * time) + (b_cond_time_arm * cond_dummy * time) + error_arm,
        order == "upper_lower" & study == -0.5 ~ (b0 + u_site_thigh + u_participant_thigh) + (b_time_thigh * time) + (b_cond_time_thigh * cond_dummy * time) + error_thigh,
        order == "upper_lower" & study == 0.5 ~ (b0 + u_site_arm + u_participant_arm) + (b_time_arm * time) + (b_cond_time_arm * cond_dummy * time) + error_arm
        
      ))|>
      
      # add outcome code
      mutate(outcome_dummy = case_when(
        (order == "lower_upper" & study == 0.5) | (order == "upper_lower" & study == -0.5) ~ -0.5, # arm
        (order == "lower_upper" & study == -0.5) | (order == "upper_lower" & study == 0.5) ~ 0.5 # thigh
      )) |>
      
      select(site_id, participant_id, order, study, cond_dummy, outcome_dummy, time, measurement, ma_z)  
    
    # simulate loss to follow
    study_2_participants <- data |>
      filter(study == 0.5) %>%
      distinct(participant_id) %>%
      slice_sample(prop = 1-dropout_prop)
    
    # Filter the original data
    data <- data %>%
      filter(
        study == -0.5 | 
          (study == 0.5 & participant_id %in% study_2_participants$participant_id)
      )
    
    data_arm <- data |>
      filter(
        (order == "lower_upper" & study == 0.5) | (order == "upper_lower" & study == -0.5)
      )
    
    data_thigh <- data |>
      filter(
        (order == "lower_upper" & study == -0.5) | (order == "upper_lower" & study == 0.5)
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

    # Add parameter info to each row
    bind_cols(params, all_tests)
  })
  
  
}

sim_separate_uncor <- function(params) {
  
  with(params, {
    # set up data structure
    data <- add_random(site_id = site_n) |>
      add_random(participant_id = participant_n, .nested_in = "site_id") |>
      add_between("participant_id", order = c("lower_upper", "upper_lower")) |>
      add_within("participant_id", study = c(-0.5,0.5)) |>
      add_between("participant_id", cond_1 = c("lower_volume", "higher_volume"), .shuffle = TRUE) |>
      add_between("participant_id", cond_2 = c("lower_volume", "higher_volume"), .shuffle = TRUE) |>
      mutate(
        cond = case_when(
          study == -0.5 ~ cond_1,
          study == 0.5 ~ cond_2
        )
      ) |>
      select(-cond_1, -cond_2) |>
      add_recode("cond", "cond_dummy", lower_volume = 0, higher_volume = 1) |>
      add_within("participant_id", time = c(0,1), 
                 measurement = c(1,2,3)) |>
      add_ranef("site_id", u_site_arm = sd_site_arm, u_site_thigh = sd_site_thigh) |>
      add_ranef("participant_id", u_participant_arm = sd_participant_arm, u_participant_thigh = sd_participant_thigh) |>
      add_ranef(error_arm = sigma_arm,
                error_thigh = sigma_thigh) |>
      # rescale the standardised effects to the current variance estimates
      mutate(b_time_arm = b_time * (sd_participant_arm^2 + sd_site_arm^2 + error_arm^2), 
             b_cond_time_arm = b_cond_time * (sd_participant_arm^2 + sd_site_arm^2 + error_arm^2),
             b_time_thigh = b_time * (sd_participant_thigh^2 + sd_site_thigh^2 + error_thigh^2), 
             b_cond_time_thigh = b_cond_time * (sd_participant_thigh^2 + sd_site_thigh^2 + error_thigh^2)) |>
      # calculate outcomes based on assumed model
      mutate(ma_z = case_when(
        order == "lower_upper" & study == -0.5 ~ (b0 + u_site_thigh + u_participant_thigh) + (b_time_thigh * time) + (b_cond_time_thigh * cond_dummy * time) + error_thigh,
        order == "lower_upper" & study == 0.5 ~ (b0 + u_site_arm + u_participant_arm) + (b_time_arm * time) + (b_cond_time_arm * cond_dummy * time) + error_arm,
        order == "upper_lower" & study == -0.5 ~ (b0 + u_site_thigh + u_participant_thigh) + (b_time_thigh * time) + (b_cond_time_thigh * cond_dummy * time) + error_thigh,
        order == "upper_lower" & study == 0.5 ~ (b0 + u_site_arm + u_participant_arm) + (b_time_arm * time) + (b_cond_time_arm * cond_dummy * time) + error_arm
        
      ))|>
      
      # add outcome code
      mutate(outcome_dummy = case_when(
        (order == "lower_upper" & study == 0.5) | (order == "upper_lower" & study == -0.5) ~ -0.5, # arm
        (order == "lower_upper" & study == -0.5) | (order == "upper_lower" & study == 0.5) ~ 0.5 # thigh
      )) |>
      
      select(site_id, participant_id, order, study, cond_dummy, outcome_dummy, time, measurement, ma_z)  
    
    # simulate loss to follow
    study_2_participants <- data |>
      filter(study == 0.5) %>%
      distinct(participant_id) %>%
      slice_sample(prop = 1-dropout_prop)
    
    # Filter the original data
    data <- data %>%
      filter(
        study == -0.5 | 
          (study == 0.5 & participant_id %in% study_2_participants$participant_id)
      )
    
    data_arm <- data |>
      filter(
        (order == "lower_upper" & study == 0.5) | (order == "upper_lower" & study == -0.5)
      )
    
    data_thigh <- data |>
      filter(
        (order == "lower_upper" & study == -0.5) | (order == "upper_lower" & study == 0.5)
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
    
    # Add parameter info to each row
    bind_cols(params, all_tests)
  })
  
  
}


generate_simulation_grid_lp <- function() {
  
  crossing(rep = 1:1000, # number of replicates
           participant_n = seq(5,100, by = 5), # range of participant N
           site_n = 15, # fixed site N
           time_n = c(0.5,1), # to examine inclusion of midpoint test or not
           measurements_n = c(1:3), # range of measurements N
  ) |>
    mutate(
      b0 = 0, b_time = 0.05, b_cond = 0, b_cond_time = 0.025,         # fixed effects 
      u_site = 0.205,  # random intercept site
      u_participant_arm = 0.629,   # random intercept participant arm
      u_participant_thigh = 0.597,  # random intercept participant thigh
      arm_ma_error = 0.166,           # measurement error arm
      thigh_ma_error = 0.198         # measurement error thigh  
    ) 
}

generate_simulation_grid_separate_old_lp_vars <- function() {
  
   crossing(
    rep = 1:100,
    participant_n = seq(5, 100, by = 5),
    site_n = 15,
    measurements_n = 3,
    dropout_prop = c(0,0.15)
  ) %>%
    mutate(
      b0 = 0, b_time = 0.05, b_cond_time = 0.025,
      sd_site_arm = 0.205,
      sd_participant_arm = 0.629,
      sigma_arm = 0.166,
      sd_site_thigh = 0.205,
      sd_participant_thigh = 0.597,
      sigma_thigh = 0.198,
      u_participant_cor = 0,
      u_site_cor = 0
      
    )
}

# plot_power <- function(sim_result) {
#   
#   sim_result_power <- sim_result |> 
#     bind_rows() |>
#     filter(term == "time:cond_dummy") |>
#     pivot_longer(c("p.value", "p.value.equiv"), 
#                  names_to = "test", 
#                  values_to = "p.value") |> 
#     mutate(test = case_when(test == "p.value" ~ "Difference",
#                             test == "p.value.equiv" ~ "Equivalence")) |>
#     group_by(participant_n, dropout_prop, test) |>
#     summarise(total = sum(p.value < .01),
#               power = mean(p.value < .01),
#               ci.lower = prop.test(total, 1000)$conf.int[1],
#               ci.upper = prop.test(total, 1000)$conf.int[2],
#               .groups = "drop")
#   
#   sim_result_power |>
#     mutate(dropout_label = "Study 2 Dropout Proportion") |>
#     ggplot(aes(x = participant_n, y = power)) +
#     geom_hline(yintercept = c(0.95,0.8), linetype = "dashed") +
#     geom_pointrange(aes(ymin = ci.lower, ymax = ci.upper), size = 0.1) +
#     geom_line() +
#     ggh4x::facet_nested(dropout_label + dropout_prop ~ test) +
#     labs(
#       x = "Total Number of Participants per Site (N)",
#       y = "Power (1-\u03b2)",
#       color = "Number of Measurements per Timepoint (N)",
#       title = "Sample Size Estimation (N per Site)",
#       subtitle = expression(~italic(y)[ijtk]~" = (\u03b2"[0]~"+ "~u[k]~" + "~u[i]~") + \u03b2"[1]~"time"[t]~" + \u03b2"[2]~"condition:time"[jt]~" + \u03f5"[ijtk]),
#       caption = "Simulated to detect a standardised effect, or equivalence, for condition:time of 0.025 with a SESOI of 0.1\nSite number k fixed at 15, \u03b1 = 0.01, 1000 simulations per sample size"
#     ) +
#     theme_bw() +
#     theme(legend.position = "bottom")
#   
# }




### testing old lp sims

generate_simulation_grid_lp <- function() {
  
  crossing(rep = 1:1000, # number of replicates
  participant_n = seq(5,100, by = 5), # range of participant N
  site_n = 15, # fixed site N
  time_n = c(0.5,1), # to examine inclusion of midpoint test or not
  measurements_n = c(1:3), # range of measurements N
  ) |>
  mutate(
    b0 = 0, b_time = 0.05, b_cond = 0, b_cond_time = 0.025,         # fixed effects 
    u_site = 0.205,  # random intercept site
    u_participant_arm = 0.629,   # random intercept participant arm
    u_participant_thigh = 0.597,  # random intercept participant thigh
    arm_ma_error = 0.166,           # measurement error arm
    thigh_ma_error = 0.198         # measurement error thigh  
  ) 
}


sim_lp_change <- function(params) {

  with(params, {
    # set up data structure
    dat <- add_random(site = site_n) |>
      add_random(participant = participant_n, .nested_in = "site") |>
      add_between("participant", cond = c("full_ROM", "lengthened_partial")) |>
      add_recode("cond", "cond_dummy", full_ROM = 0, lengthened_partial = 1) |>
      add_within("participant", time = seq(0,1, by = time_n), 
                 measurement = seq(1:measurements_n)) |>
      add_ranef("site", u_site = u_site) |>
      add_ranef("participant", u_participant_arm = u_participant_arm) |>
      add_ranef("participant", u_participant_thigh = u_participant_thigh) |>
      add_ranef(arm_ma_error = arm_ma_error) |>
      add_ranef(thigh_ma_error = thigh_ma_error) |>
      mutate(b_time_arm = b_time * (u_participant_arm^2 + u_site^2 + arm_ma_error^2), 
             b_cond_time_arm = b_cond_time * (u_participant_arm^2 + u_site^2 + arm_ma_error^2),
             b_time_thigh = b_time * (u_participant_thigh^2 + u_site^2 + thigh_ma_error^2), 
             b_cond_time_thigh = b_cond_time * (u_participant_thigh^2 + u_site^2 + thigh_ma_error^2)) |>
      mutate(arm_ma = (b0 + u_site + u_participant_arm) + (b_time_arm * time) + (b_cond * cond_dummy) + (b_cond_time_arm * cond_dummy * time) + arm_ma_error,
             thigh_ma = (b0 + u_site + u_participant_thigh) + (b_time_thigh * time) + (b_cond * cond_dummy) + (b_cond_time_thigh * cond_dummy * time) + thigh_ma_error) 
    

    # run mixed effect model and return relevant values
    model_arm <- lmer(arm_ma ~ time + cond_dummy + time:cond_dummy + (1|site) + (1 | participant),
                      data = dat, REML = TRUE)

    model_thigh <- lmer(thigh_ma ~ time + cond_dummy + time:cond_dummy + (1|site) + (1 | participant),
                        data = dat, REML = TRUE)

    
    # if (time_n == 1 && measurements_n == 1) {
    #   # run mixed effect model and return relevant values
    #   model_arm <- lmer(arm_ma ~ time*cond_dummy + (1|site) + (1 | participant),
    #                     data = dat, REML = TRUE)
    #   
    #   model_thigh <- lmer(thigh_ma ~ time*cond_dummy + (1|site) + (1 | participant),
    #                       data = dat, REML = TRUE)
    # } else {
    #   # run mixed effect model and return relevant values
    #   model_arm <- lmer(arm_ma ~ time*cond_dummy + (1|site) + (time | participant),
    #                     data = dat, REML = TRUE)
    #   
    #   model_thigh <- lmer(thigh_ma ~ time*cond_dummy + (1|site) + (time | participant),
    #                       data = dat, REML = TRUE)
    # }
    
    hypothesis_tests_arm <- hypotheses(model_arm, equivalence = c(NA,0.1), 
                                       # vcov = "satterthwaite", 
                                       df = insight::get_df(model_arm))[4,] |>
      mutate(muscle_site = "arm")
    
    
    hypothesis_tests_thigh <- hypotheses(model_thigh, equivalence = c(NA,0.1), 
                                         # vcov = "satterthwaite", 
                                         df = insight::get_df(model_thigh))[4,] |>
      mutate(muscle_site = "thigh")
    
    tests <- bind_rows(hypothesis_tests_arm, hypothesis_tests_thigh)
    
    bind_cols(params, tests)
  })
}

sim_lp_ancova <- function(params) {
  
  with(params, {
    # set up data structure
    dat <- add_random(site = site_n) |>
      add_random(participant = participant_n, .nested_in = "site") |>
      add_between("participant", cond = c("full_ROM", "lengthened_partial")) |>
      add_recode("cond", "cond_dummy", full_ROM = 0, lengthened_partial = 1) |>
      add_within("participant", time = seq(0,1, by = time_n), 
                 measurement = seq(1:measurements_n)) |>
      add_ranef("site", u_site = u_site) |>
      add_ranef("participant", u_participant_arm = u_participant_arm) |>
      add_ranef("participant", u_participant_thigh = u_participant_thigh) |>
      add_ranef(arm_ma_error = arm_ma_error) |>
      add_ranef(thigh_ma_error = thigh_ma_error) |>
      mutate(b_time_arm = b_time * (u_participant_arm^2 + u_site^2 + arm_ma_error^2), 
             b_cond_time_arm = b_cond_time * (u_participant_arm^2 + u_site^2 + arm_ma_error^2),
             b_time_thigh = b_time * (u_participant_thigh^2 + u_site^2 + thigh_ma_error^2), 
             b_cond_time_thigh = b_cond_time * (u_participant_thigh^2 + u_site^2 + thigh_ma_error^2)) |>
      mutate(arm_ma = (b0 + u_site + u_participant_arm) + (b_time_arm * time) + (b_cond * cond_dummy) + (b_cond_time_arm * cond_dummy * time) + arm_ma_error,
             thigh_ma = (b0 + u_site + u_participant_thigh) + (b_time_thigh * time) + (b_cond * cond_dummy) + (b_cond_time_thigh * cond_dummy * time) + thigh_ma_error) 
    
    
    # run mixed effect model and return relevant values
    model_arm <- lmer(arm_ma ~ time + time:cond_dummy + (1|site) + (1 | participant),
                      data = dat, REML = TRUE)
    
    model_thigh <- lmer(thigh_ma ~ time + time:cond_dummy + (1|site) + (1 | participant),
                        data = dat, REML = TRUE)
    
    
    # if (time_n == 1 && measurements_n == 1) {
    #   # run mixed effect model and return relevant values
    #   model_arm <- lmer(arm_ma ~ time*cond_dummy + (1|site) + (1 | participant),
    #                     data = dat, REML = TRUE)
    #   
    #   model_thigh <- lmer(thigh_ma ~ time*cond_dummy + (1|site) + (1 | participant),
    #                       data = dat, REML = TRUE)
    # } else {
    #   # run mixed effect model and return relevant values
    #   model_arm <- lmer(arm_ma ~ time*cond_dummy + (1|site) + (time | participant),
    #                     data = dat, REML = TRUE)
    #   
    #   model_thigh <- lmer(thigh_ma ~ time*cond_dummy + (1|site) + (time | participant),
    #                       data = dat, REML = TRUE)
    # }
    
    hypothesis_tests_arm <- hypotheses(model_arm, equivalence = c(NA,0.1), 
                                       # vcov = "satterthwaite", 
                                       df = insight::get_df(model_arm))[3,] |>
      mutate(muscle_site = "arm")
    
    
    hypothesis_tests_thigh <- hypotheses(model_thigh, equivalence = c(NA,0.1), 
                                         # vcov = "satterthwaite", 
                                         df = insight::get_df(model_thigh))[3,] |>
      mutate(muscle_site = "thigh")
    
    tests <- bind_rows(hypothesis_tests_arm, hypothesis_tests_thigh)
    
    bind_cols(params, tests)
  })
}

sim_lp_all_ancova <- function(params) {
  
  with(params, {
    # set up data structure
    dat <- add_random(site = site_n) |>
      add_random(participant = participant_n, .nested_in = "site") |>
      # add study replicate
      add_between("participant", order = c("lower_upper", "upper_lower")) |>
      add_within("participant", study = c(1,2)) |>
      add_between("participant", cond_1 = c("lower_volume", "higher_volume"), .shuffle = TRUE) |>
      add_between("participant", cond_2 = c("lower_volume", "higher_volume"), .shuffle = TRUE) |>
      mutate(
        cond = case_when(
          study == 1 ~ cond_1,
          study == 2 ~ cond_2
        )
      ) |>
      select(-cond_1, -cond_2) |>
      add_recode("cond", "cond_dummy", lower_volume = 0, higher_volume = 1) |>
      add_within("participant", time = seq(0,1, by = time_n), 
                 measurement = seq(1:measurements_n)) |>
      add_ranef("site", u_site = u_site) |>
      add_ranef("participant", u_participant_arm = u_participant_arm) |>
      add_ranef("participant", u_participant_thigh = u_participant_thigh) |>
      add_ranef(arm_ma_error = arm_ma_error) |>
      add_ranef(thigh_ma_error = thigh_ma_error) |>
      mutate(b_time_arm = b_time * (u_participant_arm^2 + u_site^2 + arm_ma_error^2), 
             b_cond_time_arm = b_cond_time * (u_participant_arm^2 + u_site^2 + arm_ma_error^2),
             b_time_thigh = b_time * (u_participant_thigh^2 + u_site^2 + thigh_ma_error^2), 
             b_cond_time_thigh = b_cond_time * (u_participant_thigh^2 + u_site^2 + thigh_ma_error^2)) |>
      mutate(arm_ma = case_when(
        order == "lower_upper" & study == 2 ~ (b0 + u_site + u_participant_arm) + (b_time_arm * time) + 0 + (b_cond_time_arm * cond_dummy * time) + arm_ma_error,
        order == "upper_lower" & study == 1 ~ (b0 + u_site + u_participant_arm) + (b_time_arm * time) + 0 + (b_cond_time_arm * cond_dummy * time) + arm_ma_error,
        .default = NA
      ),
      thigh_ma = case_when(
        order == "lower_upper" & study == 1 ~ (b0 + u_site + u_participant_thigh) + (b_time_thigh * time) + 0 + (b_cond_time_thigh * cond_dummy * time) + thigh_ma_error,
        order == "upper_lower" & study == 2 ~ (b0 + u_site + u_participant_thigh) + (b_time_thigh * time) + 0 + (b_cond_time_thigh * cond_dummy * time) + thigh_ma_error,
        .default = NA
      )
      )
    
    # run mixed effect model and return relevant values
    dat_long <- dat |>
      select(participant, site, time, cond_dummy, arm_ma, thigh_ma) |>
      pivot_longer(5:6,
                   names_to = "muscle",
                   values_to = "ma_z") |>
      filter(!is.na(ma_z))
    
    model_all <- lmer(ma_z ~ time + time:cond_dummy + (1|site) + (1 | participant),
                      data = dat_long, REML = TRUE)
    
    
    
    # if (time_n == 1 && measurements_n == 1) {
    #   # run mixed effect model and return relevant values
    #   model_arm <- lmer(arm_ma ~ time*cond_dummy + (1|site) + (1 | participant),
    #                     data = dat, REML = TRUE)
    #   
    #   model_thigh <- lmer(thigh_ma ~ time*cond_dummy + (1|site) + (1 | participant),
    #                       data = dat, REML = TRUE)
    # } else {
    #   # run mixed effect model and return relevant values
    #   model_arm <- lmer(arm_ma ~ time*cond_dummy + (1|site) + (time | participant),
    #                     data = dat, REML = TRUE)
    #   
    #   model_thigh <- lmer(thigh_ma ~ time*cond_dummy + (1|site) + (time | participant),
    #                       data = dat, REML = TRUE)
    # }
    
    tests <- hypotheses(model_all, equivalence = c(NA,0.1), 
                                       # vcov = "satterthwaite", 
                                       df = insight::get_df(model_all))[3,]

    
    bind_cols(params, tests)
  })
}
