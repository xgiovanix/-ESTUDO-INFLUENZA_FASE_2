
#===========================================================================
# TABLE — VENTILATORY SUPPORT BY REGION AND EPIDEMIOLOGICAL PERIOD
#
# Purpose:
#   Generate a publication-ready (NEJM-style) table describing
#   ventilatory support according to epidemiological period and
#   Brazilian macro-region.
#
# Required input dataset:
#   influd_2013_2025.csv
#
# Dataset prerequisites:
#   Final analytic cohort including:
#     - adults (≥18 years)
#     - influenza-confirmed SARI cases
#     - no viral co-detection
#     - known hospital outcomes
#     - ICU admissions
#     - study periods 2013–2019 and 2023–2025
#
# Required variables:
#   ano
#   regiao
#   suporte_vent
#
# GIOVANI B COSTA 26/05/2026
#===========================================================================
#=========================================================
# PERIOD
#=========================================================

library(dplyr)
library(tidyr)
library(flextable)
library(officer)

srag_vm_supl <- read.csv("/Users/xgx/Desktop/DOUTORADO/_INFLUENZA/INFLUENZA/data_/data_final/influd_2013_2025.csv", sep = ",") 
names(srag_vm_supl)


tab <- srag_vm_supl %>%
  
  mutate(
    
    periodo = case_when(
      ano <= 2019 ~ "Pre-pandemic",
      ano >= 2023 ~ "Post-pandemic"
    ),
    
    suporte_vent = case_when(
      suporte_vent == "No" ~ "None",
      suporte_vent == "Yes, non-invasive" ~
        "Non-invasive ventilation",
      suporte_vent == "Yes, invasive" ~
        "Invasive MV"
    ),
    
    regiao = case_when(
      regiao == "Sul" ~ "South",
      regiao == "Sudeste" ~ "Southeast",
      regiao == "Centro-Oeste" ~ "Central-West",
      regiao == "Nordeste" ~ "Northeast",
      regiao == "Norte" ~ "North"
    )
  ) %>%
  
  filter(
    !is.na(suporte_vent),
    !is.na(periodo)
  )

#=========================================================
# BRAZIL + REGIONS
#=========================================================

brasil <- tab %>%
  mutate(regiao = "Brazil")

tab2 <- bind_rows(
  brasil,
  tab
)

#=========================================================
# DENOMINATOR
# region x period
#=========================================================

denom <- tab2 %>%
  
  group_by(
    regiao,
    periodo
  ) %>%
  
  summarise(
    total = n(),
    .groups = "drop"
  )

#=========================================================
# N + %
#=========================================================

tab_sum <- tab2 %>%
  
  group_by(
    regiao,
    periodo,
    suporte_vent
  ) %>%
  
  summarise(
    n = n(),
    .groups = "drop"
  ) %>%
  
  left_join(
    denom,
    by = c(
      "regiao",
      "periodo"
    )
  ) %>%
  
  mutate(
    
    pct = round(
      100*n/total,
      1
    ),
    
    cell = paste0(
      format(
        n,
        big.mark = ","
      ),
      " (",
      sprintf(
        "%.1f",
        pct
      ),
      "%)"
    ),
    
    col = paste(
      suporte_vent,
      periodo,
      sep = "_"
    )
  )

##=========================================================
# REAL WIDE
#=========================================================

tab_wide <- tab_sum %>%
  
  select(
    regiao,
    suporte_vent,
    periodo,
    cell
  ) %>%
  
  pivot_wider(
    names_from = c(
      suporte_vent,
      periodo
    ),
    values_from = cell
  ) %>%
  
  mutate(
    across(
      -regiao,
      ~replace_na(
        .x,
        "0 (0.0%)"
      )
    )
  ) %>%
  
  select(
    regiao,
    
    `None_Pre-pandemic`,
    `None_Post-pandemic`,
    
    `Non-invasive ventilation_Pre-pandemic`,
    `Non-invasive ventilation_Post-pandemic`,
    
    `Invasive MV_Pre-pandemic`,
    `Invasive MV_Post-pandemic`
  ) %>%
  
  arrange(
    factor(
      regiao,
      levels = c(
        "Brazil",
        "South",
        "Southeast",
        "Central-West",
        "Northeast",
        "North"
      )
    )
  )
#=========================================================
# TITLE
#=========================================================

titulo_tab <- paste(
  "Table X. Ventilatory support according to",
  "epidemiological period and region of hospitalization."
)

#=========================================================
# FOOTNOTE
#=========================================================

legenda_tab <- paste(
  "Data are presented as n (%).",
  "Percentages were calculated within each region",
  "and epidemiological period, excluding missing",
  "data for ventilatory support.",
  "IMV denotes invasive mechanical ventilation."
)

#=========================================================
# NEJM STYLE
#=========================================================

ft <- ft %>%
  
  font(
    fontname = "Times New Roman",
    part = "all"
  ) %>%
  
  fontsize(
    size = 12,
    part = "all"
  ) %>%
  
  bold(
    part = "header"
  ) %>%
  
  align(
    j = 1,
    align = "left",
    part = "all"
  ) %>%
  
  align(
    j = 2:7,
    align = "center",
    part = "all"
  ) %>%
  
  border_remove() %>%
  
  hline_top(
    border = fp_border(
      width = 1.2
    ),
    part = "all"
  ) %>%
  
  hline(
    i = 3,
    border = fp_border(
      width = 1
    ),
    part = "header"
  ) %>%
  
  hline_bottom(
    border = fp_border(
      width = 1.2
    ),
    part = "body"
  ) %>%
  
  autofit()

#=========================================================
# EXPORT WORD
#=========================================================

doc <- read_docx()

doc <- body_add_par(
  doc,
  titulo_tab,
  style = "Normal"
)

doc <- body_add_flextable(
  doc,
  ft
)

doc <- body_add_par(
  doc,
  legenda_tab,
  style = "Normal"
)

print(
  doc,
  target =
    "Table_ventilatory_support.docx"
)

ft
