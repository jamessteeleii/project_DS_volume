

#### 24 sites

# 12 weeks
# rep range 


library(tidyverse)

targets::tar_load(c(sim_result_pooled,
                    sim_result_separate,
                    sim_result_separate_uncor))

sims <- sim_result_separate

sim_result_power_meta_re <- sims |> 
  bind_rows() |>
  filter(term == "time:cond_dummy" & model == "meta_re") |>
  pivot_longer(c("p.value", "p.value.equiv"), 
               names_to = "test", 
               values_to = "p.value") |> 
  mutate(test = case_when(test == "p.value" ~ "Difference",
                          test == "p.value.equiv" ~ "Equivalence")) |>
  group_by(participant_n, dropout_prop, test) |>
  summarise(total = sum(p.value < .01, na.rm=TRUE),
            power = mean(p.value < .01, na.rm=TRUE),
            ci.lower = prop.test(total, 100)$conf.int[1],
            ci.upper = prop.test(total, 100)$conf.int[2],
            .groups = "drop") |>
  mutate(model = "meta_re")

sim_result_power_meta_fe <- sims |> 
  bind_rows() |>
  filter(term == "time:cond_dummy" & model == "meta_fe") |>
  pivot_longer(c("p.value", "p.value.equiv"), 
               names_to = "test", 
               values_to = "p.value") |> 
  mutate(test = case_when(test == "p.value" ~ "Difference",
                          test == "p.value.equiv" ~ "Equivalence")) |>
  group_by(participant_n, dropout_prop, test) |>
  summarise(total = sum(p.value < .01, na.rm=TRUE),
            power = mean(p.value < .01, na.rm=TRUE),
            ci.lower = prop.test(total, 100)$conf.int[1],
            ci.upper = prop.test(total, 100)$conf.int[2],
            .groups = "drop") |>
  mutate(model = "meta_fe")

sim_result_power_all <- sims |> 
  bind_rows() |>
  filter(term == "time:cond_dummy" & model == "pooled") |>
  pivot_longer(c("p.value", "p.value.equiv"), 
               names_to = "test", 
               values_to = "p.value") |> 
  mutate(test = case_when(test == "p.value" ~ "Difference",
                          test == "p.value.equiv" ~ "Equivalence")) |>
  group_by(participant_n, dropout_prop, test) |>
  summarise(total = sum(p.value < .01, na.rm=TRUE),
            power = mean(p.value < .01, na.rm=TRUE),
            ci.lower = prop.test(total, 100)$conf.int[1],
            ci.upper = prop.test(total, 100)$conf.int[2],
            .groups = "drop") |>
  mutate(model = "pooled")

sim_result_power_all_fixed <- sims |> 
  bind_rows() |>
  filter(term == "time:cond_dummy" & model == "pooled_fixed") |>
  pivot_longer(c("p.value", "p.value.equiv"), 
               names_to = "test", 
               values_to = "p.value") |> 
  mutate(test = case_when(test == "p.value" ~ "Difference",
                          test == "p.value.equiv" ~ "Equivalence")) |>
  group_by(participant_n, dropout_prop, test) |>
  summarise(total = sum(p.value < .01, na.rm=TRUE),
            power = mean(p.value < .01, na.rm=TRUE),
            ci.lower = prop.test(total, 100)$conf.int[1],
            ci.upper = prop.test(total, 100)$conf.int[2],
            .groups = "drop") |>
  mutate(model = "pooled_fixed")

sim_result_power_all_random <- sims |> 
  bind_rows() |>
  filter(term == "time:cond_dummy" & model == "pooled_fixed") |>
  pivot_longer(c("p.value", "p.value.equiv"), 
               names_to = "test", 
               values_to = "p.value") |> 
  mutate(test = case_when(test == "p.value" ~ "Difference",
                          test == "p.value.equiv" ~ "Equivalence")) |>
  group_by(participant_n, dropout_prop, test) |>
  summarise(total = sum(p.value < .01, na.rm=TRUE),
            power = mean(p.value < .01, na.rm=TRUE),
            ci.lower = prop.test(total, 100)$conf.int[1],
            ci.upper = prop.test(total, 100)$conf.int[2],
            .groups = "drop") |>
  mutate(model = "pooled_random")

sim_result_power_all_disp <- sims |> 
  bind_rows() |>
  filter(term == "time:cond_dummy" & model == "pooled_disp") |>
  pivot_longer(c("p.value", "p.value.equiv"), 
               names_to = "test", 
               values_to = "p.value") |> 
  mutate(test = case_when(test == "p.value" ~ "Difference",
                          test == "p.value.equiv" ~ "Equivalence")) |>
  group_by(participant_n, dropout_prop, test) |>
  summarise(total = sum(p.value < .01, na.rm=TRUE),
            power = mean(p.value < .01, na.rm=TRUE),
            ci.lower = prop.test(total, 100)$conf.int[1],
            ci.upper = prop.test(total, 100)$conf.int[2],
            .groups = "drop") |>
  mutate(model = "pooled_disp")

sim_result_power_singly <- sims |> 
  bind_rows() |>
  mutate(
    model = case_when(
      model == "upper" | model == "lower" ~ "single_studies",
      .default = model
    )
  ) |>
  filter(term == "time:cond_dummy" & model == "single_studies") |>
  pivot_longer(c("p.value", "p.value.equiv"), 
               names_to = "test", 
               values_to = "p.value") |> 
  mutate(test = case_when(test == "p.value" ~ "Difference",
                          test == "p.value.equiv" ~ "Equivalence")) |>
  group_by(participant_n, dropout_prop, test) |>
  summarise(total = sum(p.value < .005, na.rm=TRUE),
            power = mean(p.value < .005, na.rm=TRUE),
            ci.lower = prop.test(total, 200)$conf.int[1],
            ci.upper = prop.test(total, 200)$conf.int[2],
            .groups = "drop") |>
  mutate(model = "single_studies")


sim_result_power <- bind_rows(sim_result_power_meta_fe, sim_result_power_meta_re,
                              sim_result_power_all,
                              sim_result_power_all_random, sim_result_power_all_fixed,
                              sim_result_power_all_disp,
                              sim_result_power_singly) |>
  mutate(
    model = factor(model, levels = c("single_studies", "pooled", "pooled_fixed", "pooled_random", "pooled_disp", "meta_fe", "meta_re"))
  )

sim_result_power |>
  mutate(dropout_label = "Study 2 Dropout Proportion") |>
  ggplot(aes(x = participant_n, y = power, color = model)) +
  geom_hline(yintercept = c(0.95,0.8), linetype = "dashed") +
  # geom_pointrange(aes(ymin = ci.lower, ymax = ci.upper), size = 0.1, position = position_dodge(w=2)) +
  geom_point() +
  geom_line() +
  ggh4x::facet_nested(dropout_label + dropout_prop ~ test) +
  labs(
    x = "Total Number of Participants per Site (N)",
    y = "Power (1-\u03b2)",
    color = "Studies analysed singly with alpha adjustment,\npooled with fixed/random effect of study,\nor as fixed/random effect two-stage meta-analysis",
    title = "Sample Size Estimation (N per Site)",
    subtitle = expression(~italic(y)[ijtk]~" = (\u03b2"[0]~"+ "~u[k]~" + "~u[i]~") + \u03b2"[1]~"time"[t]~" + \u03b2"[2]~"condition:time"[jt]~" + \u03f5"[ijtk]),
    caption = "Simulated to detect a standardised effect, or equivalence, for condition:time of 0.025 with a SESOI of 0.1\nSite number k fixed at 15, \u03b1 = 0.01, 100 simulations per sample size"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")














#### lp plot
targets::tar_load(c(sim_result_lp_change,
                    sim_result_lp_ancova,
                    sim_result_lp_all_ancova))

sim_result_lp_all_ancova |>
  bind_rows() |>
  filter(term == "time:cond_dummy") |>
  mutate(time_n = case_when(time_n == 0.5 ~ "Include Midpoint",
                            time_n == 1 ~ "No Midpoint",)) |>
  mutate(time_n = factor(time_n, levels = c("No Midpoint", "Include Midpoint"))) |>
  pivot_longer(c("p.value", "p.value.nonsup"), 
               names_to = "test", 
               values_to = "p.value") |> 
  mutate(test = case_when(test == "p.value" ~ "Difference",
                          test == "p.value.nonsup" ~ "Non-superiority")) |>
  group_by(measurements_n, participant_n, time_n, test) |>
  summarise(total = sum(p.value < .005, na.rm=TRUE),
            power = mean(p.value < .005, na.rm=TRUE),
            ci.lower = prop.test(total, 2000)$conf.int[1],
            ci.upper = prop.test(total, 2000)$conf.int[2],
            .groups = "drop") |>
  ggplot(aes(x = participant_n, y = power, color = factor(measurements_n))) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  geom_point() +
  geom_line() +
  scale_color_manual(values = c("#000000", "#E69F00", "#56B4E9")) +
  facet_grid(time_n ~ test) +
  labs(
    x = "Total Number of Participants per Site (N)",
    y = "Power (1-\u03b2)",
    color = "Number of Measurements per Timepoint (N)",
    title = "Sample Size Estimation (N per Site)",
    # subtitle = expression(~italic(y)[ijtk]~" = (\u03b2"[0]~"+ "~u[k]~" + "~u[i]~") + \u03b2"[1]~"condition"[j]~" + \u03b2"[2]~"time"[t]~" + \u03b2"[3]~"condition:time"[jt]~" + \u03f5"[ijtk]),
    caption = "Simulated to detect a standardised effect, or non-superiority, for condition:time with a SESOI of 0.1\nSite number k fixed at 15, \u03b1 = 0.01 (corrected to 0.005 for two primary outcomes i.e., arm and thigh muscle area)\n1000 simulations per combination of participants, number of measurements, and timepoints"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")
