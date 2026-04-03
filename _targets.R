# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)
library(tarchetypes)
library(crew)

# Set target options:
tar_option_set(
  packages = c("here",
               "tidyverse",
               "faux",
               "lme4",
               "marginaleffects",
               "metafor",
               "glmmTMB",
               "broom.mixed",
               "ggdist",
               "patchwork",
               # "furrr",
               "glue",
               "progressr"),
  seed = 1988,  # <-- GLOBAL reproducible seed
  memory = "transient",
  format = "qs",
  garbage_collection = TRUE,
  storage = "worker",
  retrieval = "worker",
  controller = crew_controller_local(workers = 10) 
)

# Run the R scripts in the R/ folder with your custom functions:
tar_source("R/functions.R") # Source other scripts as needed.

# Replace the target list below with your own:
list(
  
  ##### Pre-registration simulations ----
  
  # First read, prepare, and check LP study data
  tar_target(
    lp_data,
    read_csv(url("https://raw.githubusercontent.com/jamessteeleii/lenthened_partial_trial/refs/heads/main/data/PD_data_hypertrophy.csv"))
  ),
  
  tar_target(
    lp_data_prep,
    prep_lp_data(lp_data)
  ),
  
  # Simulating from a DGP of parameters from separate analysis of the LP study data by each muscle and the correlation of random effects

  tar_target(
    model_arm_lp_data,
    fit_model_lp_data(lp_data_prep |> filter(muscle == "arm"))
  ),

  tar_target(
    model_thigh_lp_data,
    fit_model_lp_data(lp_data_prep |> filter(muscle == "thigh"))
  ),

  # checking what random intercept correlation is between arm/thigh
  tar_target(
    ranef_cors,
    check_ranef_cors(lp_data_prep, model_arm_lp_data, model_thigh_lp_data)
  ),

  tar_target(
    plot_ranef_cors,
    plot_check_ranef_cors(lp_data_prep, model_arm_lp_data, model_thigh_lp_data)
  ),
  
  tar_target(
    plot_ranef_cors_tiff,
    ggsave(
      plot = plot_ranef_cors,
      filename = "plots/plot_ranef_cors.tiff",
      device = "tiff",
      dpi = 300,
      w = 10, 
      h = 5
    )
  ),

  tar_target(
    sim_grid_separate,
    generate_simulation_grid_separate(model_arm_lp_data, model_thigh_lp_data, ranef_cors)
  ),

  tar_target(
    sim_result_separate,
    sim_separate(
      sim_grid_separate
    ),
    pattern = map(sim_grid_separate),
    iteration = "list"
  ),
  
  tar_target(
    sim_plot,
    plot_sim(sim_result_separate)
  ),
  
  tar_target(
    sim_plot_tiff,
    ggsave(
      plot = sim_plot,
      filename = "plots/sim_plot.tiff",
      device = "tiff",
      dpi = 300,
      w = 10, 
      h = 8
    )
  ),
  
  #### Update to sample estimate for pre-reg
  tar_target(
    sim_plot_update,
    plot_sim_update(sim_result_separate)
  ),
  
  tar_target(
    sim_plot_update_tiff,
    ggsave(
      plot = sim_plot_update,
      filename = "plots/sim_plot_update.tiff",
      device = "tiff",
      dpi = 300,
      w = 10, 
      h = 8
    )
  ),
  
  #### Final data and analysis ----
  
  # We use the lengthened partial study dataset to generate the denominators for standardising
  # This is because we recruited a much smaller sample here and so the LP data is a better estimate of sigma
  
  tar_target(
    denominators,
    get_z_denominators(lp_data_prep)
  ),
  
  tar_target(
    volume_data_file,
    here("data", "pd_volume_data.csv"),
    format = "file"
  ),
  
  tar_target(
    volume_data,
    read_csv(volume_data_file)
  ),
  
  tar_target(
    volume_data_prep,
    prep_volume_data_hypertrophy(volume_data, denominators)
  ),
  
  tar_target(
    results_hypertrophy,
    fit_model_get_results_hypertrophy(volume_data_prep)
  ),
  
  tar_target(
    plots_hypertrophy,
    plot_hypertrophy_results(volume_data_prep, results_hypertrophy)
  ),
  
  tar_target(
    plots_hypertrophy_tiff,
    ggsave(
      plot = plots_hypertrophy,
      filename = "plots/hypertrophy_plot.tiff",
      device = "tiff",
      dpi = 300,
      width = 6,
      height = 9
    )
  )
  
  
)
