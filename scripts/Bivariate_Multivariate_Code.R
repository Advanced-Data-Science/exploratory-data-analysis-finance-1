library(tidyverse)
library(Hmisc)
library(ggcorrplot)
library(GGally)
library(effectsize)
library(broom)

df <- readr::read_csv("1 Month Clean.csv")

num_cols <- df |> select(where(is.numeric)) |> names()
cat_cols <- df |> select(where(\(x) is.character(x) || is.factor(x) || is.logical(x))) |> names()

num_mat <- df |> select(all_of(num_cols)) |> as.matrix()
rc <- rcorr(num_mat, type = "pearson")

ggcorrplot(
  rc$r, type = "lower", lab = TRUE, lab_size = 3,
  p.mat = rc$P, sig.level = 0.05, insig = "blank",
  title = "Correlation heatmap (Pearson r)"
)

corr_table <-
  as_tibble(rc$r, rownames = "var_x") |>
  pivot_longer(-var_x, names_to = "var_y", values_to = "pearson_r") |>
  filter(var_x < var_y) |>
  left_join(
    as_tibble(rc$P, rownames = "var_x") |>
      pivot_longer(-var_x, names_to = "var_y", values_to = "p_value"),
    by = c("var_x","var_y")
  ) |>
  left_join(
    as_tibble(rc$n, rownames = "var_x") |>
      pivot_longer(-var_x, names_to = "var_y", values_to = "n"),
    by = c("var_x","var_y")
  ) |>
  mutate(abs_r = abs(pearson_r)) |>
  arrange(desc(abs_r))

print(corr_table |> slice_head(n = 15))

top_pairs <- corr_table |> slice_head(n = 6) |> distinct(var_x, var_y) |> slice_head(n = 3)

make_scatter <- function(x, y) {
  ggplot(df, aes_string(x = x, y = y)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE) +
    labs(title = paste0("Scatter: ", x, " vs ", y, " (with linear trend)")) +
    theme_minimal()
}

plots <- purrr::pmap(top_pairs, ~ make_scatter(..1, ..2))
plots

low_card_cats <- cat_cols[map_int(df[cat_cols], ~ n_distinct(.x, na.rm = TRUE)) <= 6]
chosen_cats <- head(low_card_cats, 2)
chosen_nums <- head(num_cols, 2)

cat_num_results <- list()
for (cvar in chosen_cats) {
  for (nvar in chosen_nums) {
    print(
      ggplot(df, aes_string(x = cvar, y = nvar)) +
        geom_boxplot(outlier.alpha = 0.5) +
        stat_summary(fun = mean, geom = "point", shape = 21, size = 2) +
        labs(title = paste(nvar, "by", cvar), x = cvar, y = nvar) +
        theme_minimal()
    )
    d2 <- df |> select(all_of(c(cvar, nvar))) |> drop_na()
    if (n_distinct(d2[[cvar]]) >= 2) {
      fit <- aov(reformulate(cvar, nvar), data = d2)
      anova_tbl <- broom::tidy(fit)
      eta <- effectsize::eta_squared(fit, partial = FALSE)
      cat_num_results[[paste(cvar, nvar, sep = "___")]] <- list(anova = anova_tbl, eta = eta)
      print(anova_tbl)
      print(eta)
    }
  }
}

missing_pct <- df |> dplyr::summarise(across(all_of(num_cols), ~ mean(is.na(.)))) |> pivot_longer(everything())
pair_cols <- missing_pct |> arrange(value) |> slice_head(n = min(6, length(num_cols))) |> pull(name)

if (length(pair_cols) >= 3) {
  df |> select(all_of(pair_cols)) |> drop_na() |>
    GGally::ggpairs(title = "Scatterplot matrix (numeric features)")
}

if (length(num_cols) >= 3) {
  col_strength <- colSums(abs(rc$r), na.rm = TRUE)
  ycol <- names(sort(col_strength, decreasing = TRUE))[1]
  preds <- names(sort(abs(rc$r[, ycol]), decreasing = TRUE))
  preds <- preds[preds != ycol]
  if (length(preds) >= 2) {
    x1 <- preds[1]; x2 <- preds[2]
    sub <- df |> select(all_of(c(ycol, x1, x2))) |> drop_na()
    fit_int <- lm(reformulate(c(x1, x2, paste0(x1, ":", x2)), ycol), data = sub)
    summ <- broom::glance(fit_int)
    coef_tbl <- broom::tidy(fit_int)
    print(summ[, c("r.squared","adj.r.squared","nobs")])
    print(coef_tbl)
    sub <- sub |> mutate(x2_tertile = Hmisc::cut2(.data[[x2]], g = 3))
    ggplot(sub, aes_string(x = x1, y = ycol, color = "x2_tertile")) +
      geom_point(alpha = 0.7) +
      geom_smooth(method = "lm", se = FALSE) +
      labs(title = paste0(ycol, " ~ ", x1, " * ", x2, " (stratified by ", x2, " tertiles)")) +
      theme_minimal()
  }
}
