# =============================================================================
# R/prep_data.R  |  WRAIR Africa Dashboard: data preparation
#
# Run once per reporting month, BEFORE rendering the dashboard:
#   Rscript R/prep_data.R
#
# Reads raw DHIS2 + RADET exports from data/raw/ (never committed to GitHub),
# validates them against the lookups in data/lookups/, and writes ONE small,
# aggregated, de-identified file: data/dashboard_data.rds.
# index.qmd only reads that .rds, so it renders in seconds and the public
# repo never holds patient-level data.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(lubridate)
  library(janitor)
})

# ---- 1. Parameters: the only block to edit each month ----------------------
params <- list(
  fy                 = 2026,               # PEPFAR FY (Oct of FY-1 to Sep of FY)
  report_month       = ymd("2026-09-01"),  # last month of data to include
  program_sites_only = TRUE,               # FALSE = include every DHIS2 org unit
  files = list(
    monthly = "data/raw/prog_data_fy26_monthly.csv",  # prog_data_month_yr template
    targets = "data/raw/fy26_targets.csv",            # fyxx_targets template
    radet   = "data/raw/radet_sep26.csv"              # radet_mmmyy template
  ),
  out_rds        = "data/dashboard_data.rds",
  out_radet_csv  = "data/private/radet_clean.csv"     # patient-level, gitignored
)

fy_start   <- ymd(paste0(params$fy - 1, "-10-01"))
cutoff_dt  <- ceiling_date(params$report_month, "month") - days(1)
months_in  <- interval(fy_start, params$report_month) %/% months(1) + 1
pct_year   <- months_in / 12
issues     <- character()            # collected data-quality warnings
flag <- function(...) {
  msg <- paste0(...)
  issues <<- c(issues, msg)
  message("! ", msg)
}

# ---- 2. Lookups -------------------------------------------------------------
# Names are matched ignoring case and extra spaces ("Had a Viral Load" = "Had a viral load")
name_key <- function(x) str_to_lower(str_squish(x))

ind_map  <- read_csv("data/lookups/indicator_map.csv", show_col_types = FALSE) |>
  mutate(dhis2_name = str_squish(dhis2_name), name_key = name_key(dhis2_name))
site_map <- read_csv("data/lookups/site_map.csv", show_col_types = FALSE) |>
  mutate(across(c(org_unit, radet_name), str_squish))

# ---- 3. DHIS2 monthly program data ------------------------------------------
monthly <- read_csv(
  params$files$monthly,
  col_types = cols(.default = col_character(), Value = col_double())
) |>
  clean_names() |>
  transmute(
    dhis2_name = str_squish(data),
    name_key   = name_key(data),
    org_unit   = str_squish(organisation_unit),
    period     = my(period),           # "Oct-25" -> 2025-10-01
    value
  )

if (anyNA(monthly$period)) flag("Some Period values could not be parsed as Mon-YY.")

unmapped <- unique(monthly$dhis2_name[!monthly$name_key %in% ind_map$name_key])
if (length(unmapped)) {
  flag("Indicators missing from indicator_map.csv (ignored): ",
       paste(unmapped, collapse = " | "))
}
if (params$program_sites_only) {
  missing_sites <- setdiff(site_map$org_unit, monthly$org_unit)
  if (length(missing_sites)) flag("Program sites with no DHIS2 rows: ",
                                  paste(missing_sites, collapse = " | "))
}

prog <- monthly |>
  inner_join(select(ind_map, -dhis2_name), by = "name_key", relationship = "many-to-one") |>
  left_join(select(site_map, org_unit, site, state),
            by = "org_unit", relationship = "many-to-one")

if (params$program_sites_only) prog <- filter(prog, !is.na(site))
prog <- mutate(prog, site = coalesce(site, org_unit))

prog_fy <- filter(prog, period >= fy_start, period <= params$report_month)

# Result per data element per site:
#   cumulative indicators = sum of months so far this FY
#   snapshot indicators (agg == "latest") = value in the reporting month only
site_results <- prog_fy |>
  summarise(
    result = if (first(agg) == "latest") {
      sum(value[period == params$report_month], na.rm = TRUE)
    } else {
      sum(value, na.rm = TRUE)
    },
    .by = c(mer, dhis2_name, site, state)
  ) |>
  summarise(result = sum(result), .by = c(mer, site, state))   # e.g. OVC_SERV = CP + PR

# ---- 4. Targets -------------------------------------------------------------
targets <- read_csv(params$files$targets,
                    col_types = cols(Target = "c", Value = "d")) |>
  clean_names() |>
  transmute(
    target_key = str_squish(str_remove(target, "\\s*Target\\b.*$")),
    # TB_PREV has three target rows: the two denominator disaggregates
    # (Already on ART + New on ART) add up to the TB_PREV_D target;
    # the "completed therapy" row is the numerator.
    target_key = case_when(
      str_detect(target, regex("^TB_PREV Target:\\s*(Already|New) on ART", ignore_case = TRUE)) ~ "TB_PREV_D",
      str_detect(target, regex("^TB_PREV Target:.*completed", ignore_case = TRUE)) ~ "TB_PREV_N",
      TRUE ~ target_key
    ),
    target = value
  ) |>
  summarise(n = n(),
            target = if (all(is.na(target))) NA_real_ else sum(target, na.rm = TRUE),
            .by = target_key)

dupes <- targets |> filter(n > 1, target_key != "TB_PREV_D")
if (nrow(dupes)) flag("Duplicate target rows were added together: ",
                      paste(dupes$target_key, collapse = ", "))
targets <- select(targets, -n)

mer_info <- ind_map |> distinct(mer, target_key, agg, group, order) |>
  distinct(mer, .keep_all = TRUE)

performance <- site_results |>
  summarise(result = sum(result), .by = mer) |>
  left_join(mer_info, by = "mer") |>
  left_join(targets, by = "target_key") |>
  mutate(
    achievement = result / target,
    # Annual report: achievement is compared with the full-year target.
    # (For a mid-year report, divide achievement by `expected` to judge pace:
    #  cumulative indicators should be at months elapsed / 12 of target.)
    expected    = if_else(agg == "latest", 1, pct_year),
    status      = case_when(
      is.na(target) | target == 0 ~ "No target",
      achievement >= 0.90 ~ "90% and above",
      achievement >= 0.75 ~ "75-89%",
      achievement >= 0.50 ~ "50-74%",
      TRUE                ~ "Below 50%"
    )
  ) |>
  arrange(order)

# Monthly trend across the fiscal year (program level)
trend <- prog_fy |>
  summarise(value = sum(value, na.rm = TRUE), .by = c(mer, period)) |>
  arrange(mer, period)

# ---- 5. RADET (patient level -> aggregates only) ----------------------------
parse_dt <- function(x) {
  as_date(parse_date_time(x, orders = c("dby", "dbY", "Ymd", "dmY", "dmy"),
                          quiet = TRUE))
}
num <- function(x) suppressWarnings(parse_number(x))

radet <- read_csv(params$files$radet,
                  col_types = cols(.default = col_character()),
                  na = c("", "NA", "N/A", "null", "NULL")) |>
  clean_names() |>
  rename(sbp = systolic_blood_pressure_mm_hg,
         dbp = diastolic_blood_pressure_mm_hg) |>
  mutate(
    across(c(art_start_date, last_dispensed_date, date_of_sample_collection,
             date_of_birth), parse_dt),
    across(c(age, arv_refill_days, sbp, dbp), num),
    facility_key = str_squish(facility_name),

    # VL result: text such as "TND" / "Not detected" counts as a result of 0
    vl_num = if_else(
      str_detect(str_to_lower(coalesce(current_vl_result, "")),
                 "tnd|not detect|undetect"),
      0, num(current_vl_result)
    ),

    patient_status = case_match(
      patient_status,
      "Dead"              ~ "Known death",
      "Lost to follow up" ~ "Loss to follow up",
      "Self Transfer"     ~ "Self-transfer",
      "Treatment stop"    ~ "Stopped treatment",
      .default = patient_status
    ),

    # Active = drugs (+28-day grace) last to or beyond the cut-off date
    active = coalesce(
      last_dispensed_date <= cutoff_dt &
        last_dispensed_date + arv_refill_days + 28 >= cutoff_dt &
        !patient_status %in% c("Known death", "Transferred out", "Stopped treatment"),
      FALSE),
    vl_eligible = coalesce(active & as.numeric(cutoff_dt - art_start_date) > 180, FALSE),
    vl_current  = coalesce(vl_eligible & !is.na(vl_num) &
                             as.numeric(cutoff_dt - date_of_sample_collection) <= 365,
                           FALSE),
    supp_1000 = vl_current & coalesce(vl_num < 1000, FALSE),
    supp_50   = vl_current & coalesce(vl_num < 50, FALSE),
    vl_band = case_when(
      !vl_current   ~ NA_character_,
      vl_num < 50   ~ "<50",
      vl_num < 200  ~ "50-199",
      vl_num < 1000 ~ "200-999",
      TRUE          ~ "1000+"
    ),

    # Age band from numeric age: immune to Excel turning "5-9" into "09-May"
    age_band = cut(age, right = FALSE,
                   breaks = c(-Inf, 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, Inf),
                   labels = c("<1", "1-4", "5-9", "10-14", "15-19", "20-24",
                              "25-29", "30-34", "35-39", "40-44", "45-49", "50+")),
    sex = str_to_title(sex),

    regimen_line = case_when(
      str_detect(regimen_line, regex("1st|first", ignore_case = TRUE))  ~ "1st Line",
      str_detect(regimen_line, regex("2nd|second", ignore_case = TRUE)) ~ "2nd Line",
      str_detect(regimen_line, regex("3rd|third", ignore_case = TRUE))  ~ "3rd Line",
      TRUE ~ regimen_line
    ),
    mmd = case_when(
      is.na(arv_refill_days)  ~ NA_character_,
      arv_refill_days < 90    ~ "<3 months",
      arv_refill_days < 180   ~ "3-5 months",
      TRUE                    ~ "6+ months"
    ),
    service_status = coalesce(service_status, "Civilian"),

    # Implausible readings -> missing, then classify from most severe down
    sbp = if_else(between(sbp, 60, 300), sbp, NA_real_),
    dbp = if_else(between(dbp, 30, 200), dbp, NA_real_),
    bp_status = case_when(
      is.na(sbp) | is.na(dbp)      ~ NA_character_,
      sbp >= 180 | dbp >= 110      ~ "Grade 3 hypertension",
      sbp >= 160 | dbp >= 100      ~ "Grade 2 hypertension",
      sbp >= 140 & dbp < 90        ~ "Isolated systolic hypertension",
      sbp >= 140 | dbp >= 90       ~ "Grade 1 hypertension",
      sbp >= 130 | dbp >= 85       ~ "High normal",
      sbp >= 120 | dbp >= 80       ~ "Normal",
      TRUE                         ~ "Optimal"
    )
  ) |>
  left_join(site_map |> filter(!is.na(radet_name)) |>
              select(radet_name, site, state_site = state),
            by = c("facility_key" = "radet_name"), relationship = "many-to-one")

unmatched_fac <- radet |> filter(is.na(site)) |> distinct(facility_name) |> pull()
if (length(unmatched_fac)) flag("RADET facilities not in site_map.csv radet_name: ",
                                paste(unmatched_fac, collapse = " | "))
radet <- mutate(radet, site = coalesce(site, facility_key))

if (any(is.na(radet$last_dispensed_date) & !is.na(radet$arv_refill_days))) {
  flag(sum(is.na(radet$last_dispensed_date)), " RADET rows have an unreadable last dispensed date.")
}

cascade_cols <- function(df) {
  summarise(df,
    tx_curr     = sum(active),
    vl_eligible = sum(vl_eligible),
    vl_current  = sum(vl_current),
    supp_1000   = sum(supp_1000),
    supp_50     = sum(supp_50),
    .groups = "drop"
  ) |>
    mutate(vl_coverage    = vl_current / vl_eligible,
           vl_suppression = supp_1000 / vl_current,
           undetectable   = supp_50 / vl_current)
}

cascade       <- cascade_cols(radet)
cascade_site  <- radet |> group_by(site) |> cascade_cols() |> filter(tx_curr > 0)
vl_age_sex    <- radet |> filter(sex %in% c("Male", "Female"), !is.na(age_band)) |>
  group_by(sex, age_band) |> cascade_cols()
vl_bands_site <- radet |> filter(vl_current) |> count(site, vl_band)
mmd           <- radet |> filter(active, !is.na(mmd)) |> count(mmd)
bp            <- radet |> filter(active) |>
  summarise(n_active = n(), n_bp = sum(!is.na(bp_status))) |>
  bind_cols(radet |> filter(active, !is.na(bp_status)) |> count(bp_status) |>
              pivot_wider(names_from = bp_status, values_from = n))

# Patient-level cleaned file for local analysis only (gitignored folder).
# NAs stay NA: replacing them with 0 would invent BP and VL readings.
dir.create(dirname(params$out_radet_csv), showWarnings = FALSE, recursive = TRUE)
write_csv(radet, params$out_radet_csv, na = "")

# ---- 6. Save the aggregated bundle ------------------------------------------
saveRDS(list(
  meta = list(fy = params$fy, report_month = params$report_month,
              cutoff = cutoff_dt, pct_year = pct_year,
              generated = Sys.time(), issues = issues),
  performance   = performance,
  site_results  = site_results,
  trend         = trend,
  cascade       = cascade,
  cascade_site  = cascade_site,
  vl_age_sex    = vl_age_sex,
  vl_bands_site = vl_bands_site,
  mmd           = mmd,
  bp            = bp
), params$out_rds)

message("Saved ", params$out_rds, " | ", length(issues), " data-quality flag(s).")
