# ============================================================
# HELPER FUNCTIONS
# ============================================================

# Apply method label recoding using centralized lookup
apply_method_labels <- function(data, label_map = method_label_map_main) {
  data %>% mutate(method = recode(method, !!!label_map))
}

# Merge benchmark results with execution times cleanly
merge_exec_times <- function(df_results, exec_times) {
  df_results %>%
    mutate(temp_key = paste0(dataset, method)) %>%
    left_join(
      exec_times %>%
        mutate(temp_key = paste0(dataset, method)) %>%
        select(temp_key, time_secs),
      by = "temp_key"
    ) %>%
    select(-temp_key)
}

# For ggplot x-axis label bolding (and optionally color)
highlight <- function(x, pat, color = "black", family = "") {
  ifelse(
    grepl(pat, x),
    glue("<b style='font-family:{family}; color:{color}'>{x}</b>"),
    x
  )
}
