# ------------------------------------------------------------------------------
# 00_check_master2018_vs_org.R
# Module: pipeline/01_data_prep/01_cleaning
# Author: Agostina Torres
# ------------------------------------------------------------------------------
# Provenance check of the resurvey data: compares MORA_Cascade_Legacy_Data_
# Master2018.xlsx (field data entry, 2014-2017) with the _Org.csv exports used
# as pipeline input, cell by cell, for the four resurvey tables.
#
# Questions:
#   1. Does each CSV only ADD records (2018-2020) to the xlsx, or does it also
#      CHANGE records (values or survey years) present in the xlsx?
#   2. Which cells differ, in which column, and how often?
#
# Method (per table):
#   - All values are compared as text after normalization: whitespace trimmed,
#     numbers rounded to 8 significant digits, dates as yyyy-mm-dd, empty
#     strings and "NA" as missing; in the date column, Excel serial numbers and
#     m/d/yy text are converted to yyyy-mm-dd. Column names are matched after lowercasing
#     and removing non-alphanumeric characters.
#   - Rows are compared as complete records on the columns shared by both
#     files (multiset comparison: duplicated rows count).
#   - Records of the xlsx with no identical CSV record are then compared with
#     each CSV record leaving one column out at a time; a match identifies the
#     single cell that differs (e.g. Year changed).
#
# Inputs : data/data_raw/**/MORA_Cascade_Legacy_Data_Master2018.xlsx
#          data/data_raw/01_raw_surveys/MORA_*_Resurvey_Org.csv
# Outputs (tmp_files/master2018_vs_org/):
#   summary.csv              per table: rows, years, identical / changed / missing
#   years.csv                rows per year in xlsx and CSV
#   column_match.csv         columns present in one file only
#   single_cell_changes.csv  records differing in exactly one cell (old/new value)
#   multi_cell_changes.csv   records changed in more than one cell, paired by
#                            identifier columns, one row per differing cell
#   unmatched_xlsx_rows.csv  xlsx records with no identical or one-cell match
# ------------------------------------------------------------------------------

library(tidyverse)
library(here)
library(readxl)

cleaning_dir <- function(...) here("pipeline", "01_data_prep", "01_cleaning", ...)
org_dir <- cleaning_dir("data", "data_raw", "01_raw_surveys")
out_dir <- cleaning_dir("tmp_files", "master2018_vs_org")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

xlsx_path <- list.files(cleaning_dir("data"), pattern = "Master2018\\.xlsx$",
                        recursive = TRUE, full.names = TRUE)
if (length(xlsx_path) != 1) stop("Expected one Master2018.xlsx under 01_cleaning/data, found ", length(xlsx_path))

# xlsx sheet -> CSV file
tables <- tribble(
  ~table,      ~sheet,       ~csv,
  "plots",     "Plots",      "MORA_Plot_Resurvey_Org.csv",
  "seedlings", "Seedlings",  "MORA_Seedlings_Resurvey_Org.csv",
  "veg",       "Veg Data",   "MORA_VegData_Resurvey_Org.csv",
  "tree",      "Tree Data",  "MORA_TreeData_Resurvey_Org.csv"
)

# ------------------------------------------------------------------------------
# Normalization
# ------------------------------------------------------------------------------
clean_names <- function(x) str_to_lower(str_remove_all(x, "[^A-Za-z0-9]"))

norm_value <- function(x) {
  if (inherits(x, c("Date", "POSIXt"))) return(format(as.Date(x), "%Y-%m-%d"))
  x <- str_squish(as.character(x))
  x[x %in% c("", "NA", "N/A", "na")] <- NA
  # dates written as m/d/yyyy or yyyy-mm-dd
  is_date <- !is.na(x) & str_detect(x, "^\\d{1,2}/\\d{1,2}/\\d{4}$|^\\d{4}-\\d{2}-\\d{2}")
  if (any(is_date)) {
    d <- suppressWarnings(as.Date(lubridate::parse_date_time(x[is_date], c("mdy", "ymd", "ymd HMS"))))
    x[is_date] <- ifelse(is.na(d), x[is_date], format(d, "%Y-%m-%d"))
  }
  num <- suppressWarnings(as.numeric(x))
  is_num <- !is.na(num) & !is_date
  x[is_num] <- format(signif(num[is_num], 8), scientific = FALSE, trim = TRUE, drop0trailing = TRUE)
  x
}

# Date columns: Excel serial numbers (days since 1899-12-30) and m/d/yy text
# are both converted to yyyy-mm-dd
norm_date <- function(x) {
  serial <- suppressWarnings(as.numeric(x))
  out <- x
  is_serial <- !is.na(serial) & serial > 20000 & serial < 60000
  out[is_serial] <- format(as.Date(serial[is_serial], origin = "1899-12-30"), "%Y-%m-%d")
  is_text <- !is.na(x) & !is_serial & str_detect(x, "^\\d{1,2}/\\d{1,2}/\\d{2,4}$")
  d <- suppressWarnings(as.Date(lubridate::parse_date_time(x[is_text], c("mdy", "mdY"))))
  out[is_text] <- ifelse(is.na(d), x[is_text], format(d, "%Y-%m-%d"))
  out
}

normalize_table <- function(df) {
  names(df) <- make.unique(clean_names(names(df)), sep = "_")
  df %>%
    mutate(across(everything(), norm_value)) %>%
    mutate(across(any_of("date"), norm_date))
}

row_key <- function(df) do.call(paste, c(df, sep = "\u001f"))

# Multiset match: for each row of a, whether an unused identical row exists in b
multiset_match <- function(key_a, key_b) {
  a <- tibble(key = key_a, i = seq_along(key_a)) %>% group_by(key) %>% mutate(k = row_number()) %>% ungroup()
  b <- tibble(key = key_b) %>% count(key, name = "n_b")
  a %>% left_join(b, by = "key") %>% arrange(i) %>% mutate(matched = !is.na(n_b) & k <= n_b) %>% pull(matched)
}

# ------------------------------------------------------------------------------
# Comparison per table
# ------------------------------------------------------------------------------
compare_table <- function(table, sheet, csv) {
  x <- read_excel(xlsx_path, sheet = sheet, guess_max = 1e5) %>% normalize_table()
  y <- read_csv(file.path(org_dir, csv), col_types = cols(.default = col_character())) %>% normalize_table()
  
  cols_common <- intersect(names(x), names(y))
  col_info <- bind_rows(tibble(table, column = setdiff(names(x), names(y)), only_in = "xlsx"),
                        tibble(table, column = setdiff(names(y), names(x)), only_in = "csv"))
  year_col <- cols_common[cols_common %in% c("year", "surveyyear")][1]
  
  xs <- x[cols_common]; ys <- y[cols_common]
  identical_row <- multiset_match(row_key(xs), row_key(ys))
  
  # One-cell differences: leave one column out and match again
  remaining <- which(!identical_row)
  one_cell <- map_dfr(cols_common, function(col) {
    if (length(remaining) == 0) return(NULL)
    others <- setdiff(cols_common, col)
    kx <- row_key(xs[remaining, others]); ky <- row_key(ys[others])
    hit <- match(kx, ky)
    tibble(table, xlsx_row = remaining[!is.na(hit)], column = col,
           xlsx_value = xs[[col]][remaining[!is.na(hit)]],
           csv_value  = ys[[col]][hit[!is.na(hit)]]) %>%
      filter(!(is.na(xlsx_value) & is.na(csv_value)), xlsx_value != csv_value | is.na(xlsx_value) != is.na(csv_value))
  }) %>%
    bind_rows(tibble(table = character(), xlsx_row = integer(), column = character(),
                     xlsx_value = character(), csv_value = character())) %>%
    distinct(table, xlsx_row, .keep_all = TRUE)
  
  unmatched <- setdiff(remaining, one_cell$xlsx_row)
  
  # Records changed in more than one cell: pair xlsx and CSV records by
  # identifier columns (year, order of entry, plot, species, subplot where
  # present) and list every differing cell
  id_cols <- intersect(c("year", "orderentered", "plotid", "species", "replicate"), cols_common)
  multi_cell <- if (length(unmatched) > 0 && length(id_cols) > 0) {
    xu <- xs[unmatched, ] %>% mutate(xlsx_row = unmatched)
    yk <- ys %>% mutate(csv_row = row_number()) %>% distinct(across(all_of(id_cols)), .keep_all = TRUE)
    xu %>%
      inner_join(yk, by = id_cols, suffix = c(".xlsx", ".csv")) %>%
      pivot_longer(-c(all_of(id_cols), xlsx_row, csv_row),
                   names_to = c("column", ".value"), names_pattern = "(.*)\\.(xlsx|csv)$") %>%
      filter(xlsx != csv | is.na(xlsx) != is.na(csv)) %>%
      transmute(table, xlsx_row, csv_row, across(all_of(id_cols)), column,
                xlsx_value = xlsx, csv_value = csv)
  } else tibble()
  
  years <- if (!is.na(year_col)) {
    full_join(count(x, year = .data[[year_col]], name = "rows_xlsx"),
              count(y, year = .data[[year_col]], name = "rows_csv"), by = "year") %>%
      mutate(table, .before = 1) %>% arrange(year)
  } else tibble()
  
  summary <- tibble(
    table, sheet, csv,
    rows_xlsx = nrow(x), rows_csv = nrow(y),
    years_xlsx = if (!is.na(year_col)) paste(sort(unique(x[[year_col]])), collapse = ", ") else NA,
    years_csv  = if (!is.na(year_col)) paste(sort(unique(y[[year_col]])), collapse = ", ") else NA,
    columns_compared = length(cols_common),
    xlsx_rows_identical_in_csv = sum(identical_row),
    xlsx_rows_one_cell_changed = nrow(one_cell),
    xlsx_rows_one_cell_changed_year = if (!is.na(year_col)) sum(one_cell$column == year_col) else NA,
    xlsx_rows_not_found        = length(unmatched),
    csv_rows_from_years_not_in_xlsx = if (!is.na(year_col)) sum(!y[[year_col]] %in% x[[year_col]]) else NA
  )
  
  list(summary = summary, years = years, cols = col_info, one_cell = one_cell, multi_cell = multi_cell,
       unmatched = xs[unmatched, ] %>% mutate(table, xlsx_row = unmatched, .before = 1))
}

res <- pmap(tables, compare_table)

summary_tbl <- map_dfr(res, "summary")
write_csv(summary_tbl,                   file.path(out_dir, "summary.csv"))
write_csv(map_dfr(res, "years"),         file.path(out_dir, "years.csv"))
write_csv(map_dfr(res, "cols"),          file.path(out_dir, "column_match.csv"))
write_csv(map_dfr(res, "one_cell"),      file.path(out_dir, "single_cell_changes.csv"))
write_csv(map_dfr(res, ~ mutate(.x$multi_cell, across(everything(), as.character))),
          file.path(out_dir, "multi_cell_changes.csv"))
write_csv(map_dfr(res, ~ mutate(.x$unmatched, across(everything(), as.character))),
          file.path(out_dir, "unmatched_xlsx_rows.csv"))

print(summary_tbl %>% select(-sheet, -csv), width = Inf)
print(map_dfr(res, "years"), n = Inf)
cat("\nOne-cell changes by column:\n")
print(map_dfr(res, "one_cell") %>% count(table, column, sort = TRUE), n = Inf)
cat("\nRecords changed in more than one cell: differing columns\n")
print(map_dfr(res, ~ mutate(.x$multi_cell, across(everything(), as.character))) %>%
        count(table, column, sort = TRUE), n = Inf)
cat("\nColumns present in one file only:\n")
print(map_dfr(res, "cols"), n = Inf)
