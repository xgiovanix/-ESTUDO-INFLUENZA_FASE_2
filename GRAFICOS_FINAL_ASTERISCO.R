#========================================================
  install.packages("flextable")
library(tidyverse)
library(lubridate)
library(janitor)
library(stringr)
library(data.table)
library(arrow)
library(gtsummary)
library(ggplot2)
library(mice)
library(skimr)
library(dplyr)
library(summarytools)
library(lme4)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(corrplot)
library(GGally)
library(janitor)
#========================================================
#PREPARACAO FILTROS E CRIACAO VARIAVEIS PARA AANALISE E GRAFICOS
#========================================================
influ <- influd_export
names(influ)

unique(influ$estado)
#============================CODIGO IBGE=====================
# dicionário UF -> código IBGE
uf_to_ibge <- c(
  "RO"="11", "AC"="12", "AM"="13", "RR"="14", "PA"="15", "AP"="16", "TO"="17",
  "MA"="21", "PI"="22", "CE"="23", "RN"="24", "PB"="25", "PE"="26", "AL"="27", "SE"="28", "BA"="29",
  "MG"="31", "ES"="32", "RJ"="33", "SP"="35",
  "PR"="41", "SC"="42", "RS"="43",
  "MS"="50", "MT"="51", "GO"="52", "DF"="53"
)

influ <- influ %>%
  mutate(
    estado = case_when(
      estado %in% names(uf_to_ibge) ~ uf_to_ibge[estado],  # converte sigla → código
      grepl("^\\d{2}$", estado) ~ estado,                 # mantém se já for código
      TRUE ~ NA_character_                                 # qualquer outro vira NA
    )
  )


unique(influ$estado)#=================FILTRANDO INFLUENZA, CODETECCAO, ANO=====================


influ_filter <- influ %>% 
  filter(influenza == 1, codeteccao == 0, ano %in% c(2013:2019, 2023:2025)) 

unique(influ_filter$evolucao)
unique(influ_filter$influenza)
unique(influ_filter$codeteccao)
unique(influ_filter$ano)
#=========================criando comorbidade categorica

influ_filter <- influ_filter %>% 
  mutate(
    comorbidade_cat = case_when(
      is.na(qtd_comorbidade) ~ "Missings",
      qtd_comorbidade == 0 ~ "0",
      qtd_comorbidade == 1 ~ "1",
      qtd_comorbidade > 1 ~ ">1"
    )
  )

tabyl(influ_filter$periodo)
#=================CRIANDO PERIODO PRE E POS PANDEMIA========== 
influ_filter <- influ_filter %>% 
  mutate(
    periodo = case_when(
      ano <= 2019 ~ "pre",
      ano >= 2023 & ano <= 2025 ~ "pos",
      TRUE ~ NA_character_
    )
    
  )



#=================CRIANDO CATEGORIA IDADE========== 

influ_filter <- influ_filter %>%
  mutate(
    faixa_etaria = cut(
      idade,
      breaks = c(18, 41, 61, 81, Inf),
      labels = c("18-40", "41-60", "61-80", ">80"),
      right = FALSE
    )
  )


#=========================================================
# BIBLIOTECAS
#=========================================================

library(dplyr)
library(ggplot2)
library(binom)
library(scales)
library(patchwork)

#=========================================================
# FUNCAO PARA PREPARAR BASE DE CFR + IC95%
#=========================================================

preparar_cfr <- function(data, var, periodo_ref){
  
  var_quo <- rlang::enquo(var)
  
  data %>% 
    filter(
      periodo == periodo_ref,
      !is.na(!!var_quo),
      !is.na(evolucao)
    ) %>% 
    group_by(!!var_quo, faixa_etaria) %>% 
    summarise(
      n = n(),
      morte = sum(evolucao == "Death"),
      CFR = ifelse(n > 0, morte/n, NA_real_),
      .groups = "drop"
    ) %>% 
    mutate(
      ic = binom.confint(morte, n, method = "wilson"),
      CFR_low = ic$lower,
      CFR_high = ic$upper
    ) %>% 
    select(-ic)
}

#=========================================================
# FUNCAO PARA GERAR GRAFICOS
#=========================================================


#=========================================================
# TESTES ESTATISTICOS
# PRE VS POST
# DENTRO DE CADA SUPORTE E FAIXA ETARIA
#=========================================================

teste_vm <- influ_filter %>%
  
  filter(
    
    faixa_etaria %in% c(
      "18-40",
      "41-60",
      "61-80",
      ">80"
    ),
    
    suporte_vent %in% c(
      "No",
      "Yes, non-invasive",
      "Yes, invasive"
    ),
    
    periodo %in% c("pre","pos"),
    
    !is.na(evolucao)
    
  ) %>%
  
  group_by(
    faixa_etaria,
    suporte_vent
  ) %>%
  
  summarise(
    
    p = {
      
      tabela <- table(
        periodo,
        evolucao
      )
      
      if(
        nrow(tabela) == 2 &
        ncol(tabela) >= 2
      ){
        
        suppressWarnings(
          chisq.test(
            tabela,
            correct = FALSE
          )$p.value
        )
        
      } else {
        
        NA_real_
        
      }
      
    },
    
    .groups = "drop"
    
  ) %>%
  
  mutate(
    
    signif = ifelse(
      p < 0.05,
      "*",
      ""
    )
    
  )

#=========================================================
# POSICAO DOS ASTERISCOS
#=========================================================

asteriscos <- vm_periodo %>%
  
  group_by(
    faixa_etaria,
    suporte_vent
  ) %>%
  
  summarise(
    
    y = max(CFR_high) + 0.04,
    
    .groups = "drop"
    
  ) %>%
  
  left_join(
    teste_vm,
    by = c(
      "faixa_etaria",
      "suporte_vent"
    )
  ) %>%
  
  mutate(
    
    faixa_num = case_when(
      faixa_etaria == "18-40" ~ 1,
      faixa_etaria == "41-60" ~ 2,
      faixa_etaria == "61-80" ~ 3,
      faixa_etaria == ">80" ~ 4
    ),
    
    deslocamento = case_when(
      
      suporte_vent == "No" ~ -0.28,
      
      suporte_vent ==
        "Yes, non-invasive" ~ 0,
      
      suporte_vent ==
        "Yes, invasive" ~ 0.28
      
    ),
    
    x = faixa_num + deslocamento
    
  )

#=========================================================
# GRAFICO FINAL
#=========================================================

grafico_vm_periodo <- ggplot(
  
  vm_periodo,
  
  aes(
    x = faixa_etaria,
    y = CFR,
    fill = grupo
  )
  
) +
  
  geom_col(
    
    position = position_dodge(width = 0.9),
    
    color = "black",
    
    width = 0.8
    
  ) +
  
  geom_errorbar(
    
    aes(
      ymin = CFR_low,
      ymax = CFR_high
    ),
    
    position = position_dodge(width = 0.9),
    
    width = 0.2,
    
    linewidth = 0.5
    
  ) +
  
  #=======================================================
# ASTERISCOS
#=======================================================

geom_text(
  
  data = asteriscos,
  
  aes(
    x = x,
    y = y,
    label = signif
  ),
  
  inherit.aes = FALSE,
  
  size = 7,
  
  family = "Times New Roman",
  
  fontface = "bold"
  
) +
  
  scale_fill_manual(
    
    values = c(
      "No - Pre" = "#aecde1",
      "No - Post" = "#2171b5",
      "Yes, non-invasive - Pre" = "#c7e9c0",
      "Yes, non-invasive - Post" = "#238b45",
      "Yes, invasive - Pre" = "#f5b689",
      "Yes, invasive - Post" = "#eb690c"
    ),
    
    labels = c(
      "No ventilatory support, pre-pandemic",
      "No ventilatory support, post-pandemic",
      "Non-invasive ventilation, pre-pandemic",
      "Non-invasive ventilation, post-pandemic",
      "Invasive mechanical ventilation, pre-pandemic",
      "Invasive mechanical ventilation, post-pandemic"
    ),
    
    name = NULL
    
  ) +
  
  scale_y_continuous(
    
    limits = c(0,1),
    
    labels = scales::percent_format(
      accuracy = 1
    )
    
  ) +
  
  labs(
    
    title =
      "In-Hospital Case Fatality Rates according to ventilatory support and age group",
    
    subtitle =
      "Adult ICU patients hospitalized with influenza-associated SARI in Brazil",
    
    x =
      "Age group (years)",
    
    y =
      "Case fatality rate (%)",
    
    caption =
      paste0(
        "Bars represent crude in-hospital case fatality rates and error bars represent 95% confidence intervals estimated using the Wilson method.\n",
        "Pre-pandemic period includes 2013–2019 and post-pandemic period includes 2023–2025.\n",
        "* p<0.05 for comparisons between pre- and post-pandemic periods within each ventilatory support category and age stratum using Pearson's chi-square test."
      )
    
  ) +
  
  theme_classic() +
  
  theme(
    
    panel.background =
      element_rect(
        fill = "white",
        color = NA
      ),
    
    plot.background =
      element_rect(
        fill = "white",
        color = NA
      ),
    
    panel.grid.major.y =
      element_line(
        color = "grey85",
        linewidth = 0.8
      ),
    
    panel.grid.minor.y =
      element_blank(),
    
    axis.line =
      element_line(
        color = "black",
        linewidth = 0.8
      ),
    
    axis.ticks =
      element_line(
        color = "black"
      ),
    
    axis.text =
      element_text(
        size = 16,
        family = "Times New Roman"
      ),
    
    axis.title.y =
      element_text(
        size = 16,
        family = "Times New Roman"
      ),
    
    axis.title.x =
      element_text(
        size = 16,
        family = "Times New Roman"
      ),
    
    plot.title =
      element_text(
        size = 16,
        face = "bold",
        family = "Times New Roman"
      ),
    
    plot.subtitle =
      element_text(
        size = 14,
        family = "Times New Roman"
      ),
    
    plot.caption =
      element_text(
        size = 14,
        family = "Times New Roman",
        hjust = 0
      ),
    
    legend.position = "top",
    
    legend.justification = "left",
    
    legend.text =
      element_text(
        size = 12,
        family = "Times New Roman"
      ),
    
    legend.title =
      element_blank()
    
  )

grafico_vm_periodo

#=========================================================
# EXPORTAR
#=========================================================

ggsave(
  "grafico_vm_periodo_final.png",
  grafico_vm_periodo,
  width = 14,
  height = 8,
  dpi = 600,
  bg = "white"
)

#=======SEXO
#=========================================================
# SEX
# CFR + IC95% + CHI-SQUARE + ASTERISCOS
#=========================================================

sexo_periodo <- influ_filter %>% 
  filter(
    sexo %in% c("Female", "Male"),
    !is.na(sexo),
    !is.na(evolucao),
    !is.na(periodo)
  ) %>% 
  group_by(periodo, sexo, faixa_etaria) %>% 
  summarise(
    n = n(),
    morte = sum(evolucao == "Death"),
    CFR = ifelse(n > 0, morte / n, NA_real_),
    .groups = "drop"
  ) %>% 
  mutate(
    ic = binom.confint(
      morte,
      n,
      method = "wilson"
    ),
    CFR_low = ic$lower,
    CFR_high = ic$upper
  ) %>% 
  select(-ic) %>% 
  mutate(
    sexo = factor(
      sexo,
      levels = c(
        "Female",
        "Male"
      )
    ),
    
    periodo = factor(
      periodo,
      levels = c("pre", "pos"),
      labels = c("Pre", "Post")
    ),
    
    grupo = interaction(
      sexo,
      periodo,
      sep = " - "
    ),
    
    grupo = factor(
      grupo,
      levels = c(
        "Female - Pre",
        "Female - Post",
        "Male - Pre",
        "Male - Post"
      )
    )
  )

#=========================================================
# TESTES ESTATISTICOS
# PRE VS POST DENTRO DE SEXO E FAIXA ETARIA
#=========================================================

teste_sexo <- influ_filter %>%
  
  filter(
    sexo %in% c(
      "Female",
      "Male"
    ),
    
    periodo %in% c(
      "pre",
      "pos"
    ),
    
    !is.na(evolucao)
  ) %>%
  
  group_by(
    faixa_etaria,
    sexo
  ) %>%
  
  summarise(
    
    p = {
      
      tabela <- table(
        periodo,
        evolucao
      )
      
      if(
        nrow(tabela) == 2 &
        ncol(tabela) >= 2
      ){
        
        suppressWarnings(
          chisq.test(
            tabela,
            correct = FALSE
          )$p.value
        )
        
      } else {
        
        NA_real_
      }
    },
    
    .groups = "drop"
  ) %>%
  
  mutate(
    signif = ifelse(
      p < 0.05,
      "*",
      ""
    )
  )

#=========================================================
# POSICAO DOS ASTERISCOS
#=========================================================

asteriscos_sexo <- sexo_periodo %>%
  
  group_by(
    faixa_etaria,
    sexo
  ) %>%
  
  summarise(
    y = max(CFR_high) + 0.04,
    .groups = "drop"
  ) %>%
  
  left_join(
    teste_sexo,
    by = c(
      "faixa_etaria",
      "sexo"
    )
  ) %>%
  
  mutate(
    
    faixa_num = case_when(
      faixa_etaria == "18-40" ~ 1,
      faixa_etaria == "41-60" ~ 2,
      faixa_etaria == "61-80" ~ 3,
      faixa_etaria == ">80" ~ 4
    ),
    
    deslocamento = case_when(
      sexo == "Female" ~ -0.18,
      sexo == "Male" ~ 0.18
    ),
    
    x = faixa_num + deslocamento
  )

#=========================================================
# GRAFICO
#=========================================================

grafico_sexo_periodo <- ggplot(
  sexo_periodo,
  aes(
    x = faixa_etaria,
    y = CFR,
    fill = grupo
  )
) +
  
  geom_col(
    position = position_dodge(width = 0.9),
    color = "black",
    width = 0.8
  ) +
  
  geom_errorbar(
    aes(
      ymin = CFR_low,
      ymax = CFR_high
    ),
    position = position_dodge(width = 0.9),
    width = 0.2,
    linewidth = 0.5
  ) +
  
  geom_text(
    data = asteriscos_sexo,
    aes(
      x = x,
      y = y,
      label = signif
    ),
    inherit.aes = FALSE,
    size = 7,
    family = "Times New Roman",
    fontface = "bold"
  ) +
  
  scale_fill_manual(
    values = c(
      "Female - Pre" = "#f4c6d7",
      "Female - Post" = "#c94f7c",
      "Male - Pre" = "#b8d4f0",
      "Male - Post" = "#2f6db3"
    ),
    labels = c(
      "Female, pre-pandemic",
      "Female, post-pandemic",
      "Male, pre-pandemic",
      "Male, post-pandemic"
    ),
    name = NULL
  ) +
  
  scale_y_continuous(
    limits = c(0,1),
    labels = scales::percent_format(
      accuracy = 1
    )
  ) +
  
  labs(
    title = "In-Hospital Case Fatality Rates according to sex and age group",
    subtitle = "Adult ICU patients hospitalized with influenza-associated SARI in Brazil",
    x = "Age group (years)",
    y = "Case fatality rate (%)",
    caption =
      paste0(
        "Bars represent crude in-hospital case fatality rates and error bars represent 95% confidence intervals estimated using the Wilson method.\n",
        "Pre-pandemic period includes 2013–2019 and post-pandemic period includes 2023–2025.\n",
        "* p<0.05 for comparisons between pre- and post-pandemic periods within each sex and age stratum using Pearson's chi-square test."
      )
  ) +
  
  theme_classic() +
  
  theme(
    panel.background = element_rect(fill="white",color=NA),
    plot.background = element_rect(fill="white",color=NA),
    
    panel.grid.major.y =
      element_line(
        color="grey85",
        linewidth=0.8
      ),
    
    panel.grid.minor.y =
      element_blank(),
    
    axis.line =
      element_line(
        color="black",
        linewidth=0.8
      ),
    
    axis.ticks =
      element_line(
        color="black"
      ),
    
    axis.text =
      element_text(
        size=12,
        family="Times New Roman"
      ),
    
    axis.title.y =
      element_text(
        size=12,
        family="Times New Roman"
      ),
    
    axis.title.x =
      element_text(
        size=12,
        family="Times New Roman"
      ),
    
    plot.title =
      element_text(
        size=14,
        face="bold",
        family="Times New Roman"
      ),
    
    plot.subtitle =
      element_text(
        size=11,
        family="Times New Roman"
      ),
    
    plot.caption =
      element_text(
        size=10,
        family="Times New Roman",
        hjust=0
      ),
    
    legend.position="top",
    legend.justification="left",
    legend.text=
      element_text(
        size=11,
        family="Times New Roman"
      ),
    legend.title=element_blank()
  )

grafico_sexo_periodo

ggsave(
  filename = "Figure_CFR_Sex_AgeGroup.png",
  plot = grafico_sexo_periodo,
  width = 13,
  height = 7,
  units = "in",
  dpi = 600,
  bg = "white"
)

#==============
#=============
#=============

#=========================================================
# COMORBIDITY BURDEN
# CFR + IC95% + CHI-SQUARE + ASTERISCOS
#=========================================================

comorb_periodo <- influ_filter %>% 
  filter(
    comorbidade_cat %in% c("0","1",">1"),
    !is.na(comorbidade_cat),
    !is.na(evolucao),
    !is.na(periodo)
  ) %>% 
  group_by(
    periodo,
    comorbidade_cat,
    faixa_etaria
  ) %>% 
  summarise(
    n = n(),
    morte = sum(evolucao == "Death"),
    CFR = ifelse(n > 0, morte / n, NA_real_),
    .groups = "drop"
  ) %>% 
  mutate(
    ic = binom.confint(
      morte,
      n,
      method = "wilson"
    ),
    CFR_low = ic$lower,
    CFR_high = ic$upper
  ) %>% 
  select(-ic) %>% 
  mutate(
    
    comorbidade_cat = factor(
      comorbidade_cat,
      levels = c(
        "0",
        "1",
        ">1"
      )
    ),
    
    periodo = factor(
      periodo,
      levels = c("pre","pos"),
      labels = c("Pre","Post")
    ),
    
    grupo = interaction(
      comorbidade_cat,
      periodo,
      sep = " - "
    ),
    
    grupo = factor(
      grupo,
      levels = c(
        "0 - Pre",
        "0 - Post",
        "1 - Pre",
        "1 - Post",
        ">1 - Pre",
        ">1 - Post"
      )
    )
  )

#=========================================================
# TESTES ESTATISTICOS
# PRE VS POST DENTRO DE COMORBIDADE E FAIXA ETARIA
#=========================================================

teste_comorb <- influ_filter %>%
  
  filter(
    comorbidade_cat %in% c(
      "0",
      "1",
      ">1"
    ),
    
    periodo %in% c(
      "pre",
      "pos"
    ),
    
    !is.na(evolucao)
  ) %>%
  
  group_by(
    faixa_etaria,
    comorbidade_cat
  ) %>%
  
  summarise(
    
    p = {
      
      tabela <- table(
        periodo,
        evolucao
      )
      
      if(
        nrow(tabela)==2 &
        ncol(tabela)>=2
      ){
        
        suppressWarnings(
          chisq.test(
            tabela,
            correct=FALSE
          )$p.value
        )
        
      } else {
        NA_real_
      }
    },
    
    .groups="drop"
  ) %>%
  
  mutate(
    signif=ifelse(
      p<0.05,
      "*",
      ""
    )
  )

#=========================================================
# POSICAO DOS ASTERISCOS
#=========================================================

asteriscos_comorb <- comorb_periodo %>%
  
  group_by(
    faixa_etaria,
    comorbidade_cat
  ) %>%
  
  summarise(
    y=max(CFR_high)+0.04,
    .groups="drop"
  ) %>%
  
  left_join(
    teste_comorb,
    by=c(
      "faixa_etaria",
      "comorbidade_cat"
    )
  ) %>%
  
  mutate(
    
    faixa_num = case_when(
      faixa_etaria=="18-40" ~ 1,
      faixa_etaria=="41-60" ~ 2,
      faixa_etaria=="61-80" ~ 3,
      faixa_etaria==">80" ~ 4
    ),
    
    deslocamento = case_when(
      comorbidade_cat=="0" ~ -0.28,
      comorbidade_cat=="1" ~ 0,
      comorbidade_cat==">1" ~ 0.28
    ),
    
    x = faixa_num + deslocamento
  )

#=========================================================
# GRAFICO
#=========================================================

grafico_comorb_periodo <- ggplot(
  comorb_periodo,
  aes(
    x = faixa_etaria,
    y = CFR,
    fill = grupo
  )
) +
  
  geom_col(
    position = position_dodge(width = 0.9),
    color = "black",
    width = 0.8
  ) +
  
  geom_errorbar(
    aes(
      ymin = CFR_low,
      ymax = CFR_high
    ),
    position = position_dodge(width = 0.9),
    width = 0.2,
    linewidth = 0.5
  ) +
  
  geom_text(
    data = asteriscos_comorb,
    aes(
      x = x,
      y = y,
      label = signif
    ),
    inherit.aes = FALSE,
    size = 7,
    family = "Times New Roman",
    fontface = "bold"
  ) +
  
  scale_fill_manual(
    values = c(
      "0 - Pre" = "#d9d9d9",
      "0 - Post" = "#7f7f7f",
      "1 - Pre" = "#c7e9c0",
      "1 - Post" = "#41ab5d",
      ">1 - Pre" = "#74c476",
      ">1 - Post" = "#005a32"
    ),
    
    labels = c(
      "No comorbidity, pre-pandemic",
      "No comorbidity, post-pandemic",
      "One comorbidity, pre-pandemic",
      "One comorbidity, post-pandemic",
      ">1 comorbidity, pre-pandemic",
      ">1 comorbidity, post-pandemic"
    ),
    
    name = NULL
  ) +
  
  scale_y_continuous(
    limits = c(0,1),
    labels = scales::percent_format(
      accuracy = 1
    )
  ) +
  
  labs(
    
    title =
      "In-Hospital Case Fatality Rates according to comorbidity burden and age group",
    
    subtitle =
      "Adult ICU patients hospitalized with influenza-associated SARI in Brazil",
    
    x =
      "Age group (years)",
    
    y =
      "Case fatality rate (%)",
    
    caption =
      paste0(
        "Bars represent crude in-hospital case fatality rates and error bars represent 95% confidence intervals estimated using the Wilson method.\n",
        "Pre-pandemic period includes 2013–2019 and post-pandemic period includes 2023–2025.\n",
        "* p<0.05 for comparisons between pre- and post-pandemic periods within each comorbidity category and age stratum using Pearson's chi-square test."
      )
  ) +
  
  theme_classic() +
  
  theme(
    panel.background=element_rect(fill="white",color=NA),
    plot.background=element_rect(fill="white",color=NA),
    
    panel.grid.major.y=
      element_line(
        color="grey85",
        linewidth=0.8
      ),
    
    panel.grid.minor.y=
      element_blank(),
    
    axis.line=
      element_line(
        color="black",
        linewidth=0.8
      ),
    
    axis.text=
      element_text(
        size=12,
        family="Times New Roman"
      ),
    
    axis.title=
      element_text(
        size=12,
        family="Times New Roman"
      ),
    
    plot.title=
      element_text(
        size=14,
        face="bold",
        family="Times New Roman"
      ),
    
    plot.subtitle=
      element_text(
        size=11,
        family="Times New Roman"
      ),
    
    plot.caption=
      element_text(
        size=10,
        family="Times New Roman",
        hjust=0
      ),
    
    legend.position="top",
    legend.justification="left",
    legend.text=
      element_text(
        size=14,
        family="Times New Roman"
      )
  )

grafico_comorb_periodo

ggsave(
  "Figure_CFR_Comorbidity_AgeGroup.png",
  grafico_comorb_periodo,
  width=13,
  height=7,
  dpi=600,
  bg="white"
)

#==============
#=============
#=============
#=========================================================
# CFR VACCINATION
#=========================================================

vacina_periodo <- influ_filter %>% 
  
  mutate(
    vacina_cat = case_when(
      vacina == 1 ~ "Vaccinated",
      vacina == 2 ~ "Unvaccinated",
      TRUE ~ "Unknown"
    )
  ) %>% 
  
  filter(
    !is.na(evolucao),
    !is.na(periodo)
  ) %>% 
  
  group_by(
    periodo,
    vacina_cat,
    faixa_etaria
  ) %>% 
  
  summarise(
    n = n(),
    morte = sum(evolucao == "Death"),
    CFR = ifelse(
      n > 0,
      morte / n,
      NA_real_
    ),
    .groups = "drop"
  ) %>% 
  
  mutate(
    ic = binom.confint(
      morte,
      n,
      method = "wilson"
    ),
    
    CFR_low = ic$lower,
    CFR_high = ic$upper
  ) %>% 
  
  select(-ic) %>% 
  
  mutate(
    
    vacina_cat = factor(
      vacina_cat,
      levels = c(
        "Unknown",
        "Unvaccinated",
        "Vaccinated"
      )
    ),
    
    periodo = factor(
      periodo,
      levels = c("pre", "pos"),
      labels = c("Pre", "Post")
    ),
    
    grupo = interaction(
      vacina_cat,
      periodo,
      sep = " - "
    ),
    
    grupo = factor(
      grupo,
      levels = c(
        "Unknown - Pre",
        "Unknown - Post",
        "Unvaccinated - Pre",
        "Unvaccinated - Post",
        "Vaccinated - Pre",
        "Vaccinated - Post"
      )
    )
  )

#=========================================================
# TESTES ESTATISTICOS
# PRE VS POST DENTRO DE VACINACAO E FAIXA ETARIA
#=========================================================

teste_vacina <- influ_filter %>%
  
  mutate(
    vacina_cat = case_when(
      vacina == 1 ~ "Vaccinated",
      vacina == 2 ~ "Unvaccinated",
      TRUE ~ "Unknown"
    )
  ) %>%
  
  filter(
    !is.na(evolucao),
    !is.na(periodo)
  ) %>%
  
  group_by(
    faixa_etaria,
    vacina_cat
  ) %>%
  
  summarise(
    
    p = {
      
      tabela <- table(
        periodo,
        evolucao
      )
      
      if(
        nrow(tabela)==2 &
        ncol(tabela)>=2
      ){
        
        suppressWarnings(
          chisq.test(
            tabela,
            correct = FALSE
          )$p.value
        )
        
      } else {
        NA_real_
      }
    },
    
    .groups="drop"
  ) %>%
  
  mutate(
    signif = ifelse(
      p < 0.05,
      "*",
      ""
    )
  )

#=========================================================
# POSICAO DOS ASTERISCOS
#=========================================================

asteriscos_vacina <- vacina_periodo %>%
  
  group_by(
    faixa_etaria,
    vacina_cat
  ) %>%
  
  summarise(
    y = max(CFR_high) + 0.04,
    .groups = "drop"
  ) %>%
  
  left_join(
    teste_vacina,
    by = c(
      "faixa_etaria",
      "vacina_cat"
    )
  ) %>%
  
  mutate(
    
    faixa_num = case_when(
      faixa_etaria == "18-40" ~ 1,
      faixa_etaria == "41-60" ~ 2,
      faixa_etaria == "61-80" ~ 3,
      faixa_etaria == ">80" ~ 4
    ),
    
    deslocamento = case_when(
      vacina_cat == "Unknown" ~ -0.28,
      vacina_cat == "Unvaccinated" ~ 0,
      vacina_cat == "Vaccinated" ~ 0.28
    ),
    
    x = faixa_num + deslocamento
  )

#=========================================================
# GRAFICO
#=========================================================

grafico_vacina_periodo <- ggplot(
  vacina_periodo,
  aes(
    x = faixa_etaria,
    y = CFR,
    fill = grupo
  )
) +
  
  geom_col(
    position = position_dodge(width = 0.9),
    color = "black",
    width = 0.8
  ) +
  
  geom_errorbar(
    aes(
      ymin = CFR_low,
      ymax = CFR_high
    ),
    position = position_dodge(width = 0.9),
    width = 0.2,
    linewidth = 0.5
  ) +
  
  #=========================================
# ASTERISCOS
#=========================================

geom_text(
  data = asteriscos_vacina,
  aes(
    x = x,
    y = y,
    label = signif
  ),
  inherit.aes = FALSE,
  size = 7,
  family = "Times New Roman",
  fontface = "bold"
) +
  
  scale_fill_manual(
    values = c(
      "Unknown - Pre" = "#d9d9d9",
      "Unknown - Post" = "#7f7f7f",
      "Unvaccinated - Pre" = "#fdd0a2",
      "Unvaccinated - Post" = "#e6550d",
      "Vaccinated - Pre" = "#6baed6",
      "Vaccinated - Post" = "#08519c"
    ),
    
    labels = c(
      "Unknown vaccination status, pre-pandemic",
      "Unknown vaccination status, post-pandemic",
      "Unvaccinated, pre-pandemic",
      "Unvaccinated, post-pandemic",
      "Vaccinated, pre-pandemic",
      "Vaccinated, post-pandemic"
    ),
    
    name = NULL
  ) +
  
  scale_y_continuous(
    limits = c(0,1),
    labels = scales::percent_format(
      accuracy = 1
    )
  ) +
  
  labs(
    title =
      "In-Hospital Case Fatality Rates according to influenza vaccination status and age group",
    
    subtitle =
      "Adult ICU patients hospitalized with influenza-associated SARI in Brazil",
    
    x =
      "Age group (years)",
    
    y =
      "Case fatality rate (%)",
    
    caption =
      paste0(
        "Bars represent crude in-hospital case fatality rates and error bars represent 95% confidence intervals estimated using the Wilson method.\n",
        "Pre-pandemic period includes 2013–2019 and post-pandemic period includes 2023–2025.\n",
        "* p<0.05 for comparisons between pre- and post-pandemic periods within each vaccination category and age stratum using Pearson's chi-square test."
      )
  ) +
  
  theme_classic() +
  
  theme(
    panel.background = element_rect(fill="white", color=NA),
    plot.background = element_rect(fill="white", color=NA),
    
    panel.grid.major.y =
      element_line(
        color="grey85",
        linewidth=0.8
      ),
    
    panel.grid.minor.y =
      element_blank(),
    
    axis.line =
      element_line(
        color="black",
        linewidth=0.8
      ),
    
    axis.ticks =
      element_line(
        color="black"
      ),
    
    axis.text =
      element_text(
        size=12,
        family="Times New Roman"
      ),
    
    axis.title.y =
      element_text(
        size=12,
        family="Times New Roman"
      ),
    
    axis.title.x =
      element_text(
        size=12,
        family="Times New Roman"
      ),
    
    plot.title =
      element_text(
        size=14,
        face="bold",
        family="Times New Roman"
      ),
    
    plot.subtitle =
      element_text(
        size=11,
        family="Times New Roman"
      ),
    
    plot.caption =
      element_text(
        size=12,
        family="Times New Roman",
        hjust=0
      ),
    
    legend.position="top",
    legend.justification="left",
    legend.text=
      element_text(
        size=12,
        family="Times New Roman"
      ),
    
    legend.title=
      element_blank()
  )

grafico_vacina_periodo

ggsave(
  filename = "Figure_CFR_Vaccination_AgeGroup.png",
  plot = grafico_vacina_periodo,
  width = 13,
  height = 7,
  units = "in",
  dpi = 600,
  bg = "white"
)
#====








#=========================================================
# CFR REGION
#=========================================================

regiao_periodo <- influ_filter %>% 
  filter(
    !is.na(regiao),
    !is.na(evolucao),
    !is.na(periodo)
  ) %>% 
  mutate(
    regiao = case_when(
      regiao == "Norte" ~ "North",
      regiao == "Nordeste" ~ "Northeast",
      regiao == "Centro-Oeste" ~ "Central-West",
      regiao == "Sudeste" ~ "Southeast",
      regiao == "Sul" ~ "South",
      TRUE ~ NA_character_
    )
  ) %>% 
  group_by(periodo, regiao, faixa_etaria) %>% 
  summarise(
    n = n(),
    morte = sum(evolucao == "Death"),
    CFR = ifelse(n > 0, morte / n, NA_real_),
    .groups = "drop"
  )

ic_regiao <- binom.confint(
  x = regiao_periodo$morte,
  n = regiao_periodo$n,
  method = "wilson"
)

regiao_periodo <- regiao_periodo %>% 
  mutate(
    CFR_low = ic_regiao$lower,
    CFR_high = ic_regiao$upper,
    
    regiao = factor(
      regiao,
      levels = c(
        "North",
        "Northeast",
        "Central-West",
        "Southeast",
        "South"
      )
    ),
    
    periodo = factor(
      periodo,
      levels = c("pre","pos"),
      labels = c("Pre","Post")
    ),
    
    grupo = interaction(
      regiao,
      periodo,
      sep = " - "
    ),
    
    grupo = factor(
      grupo,
      levels = c(
        "North - Pre","North - Post",
        "Northeast - Pre","Northeast - Post",
        "Central-West - Pre","Central-West - Post",
        "Southeast - Pre","Southeast - Post",
        "South - Pre","South - Post"
      )
    )
  )

#=========================================================
# TESTE CHI-SQUARE
#=========================================================

teste_regiao <- influ_filter %>%
  
  filter(
    !is.na(regiao),
    !is.na(evolucao),
    !is.na(periodo)
  ) %>%
  
  mutate(
    regiao = case_when(
      regiao=="Norte" ~ "North",
      regiao=="Nordeste" ~ "Northeast",
      regiao=="Centro-Oeste" ~ "Central-West",
      regiao=="Sudeste" ~ "Southeast",
      regiao=="Sul" ~ "South",
      TRUE ~ NA_character_
    )
  ) %>%
  
  group_by(
    faixa_etaria,
    regiao
  ) %>%
  
  summarise(
    
    p = {
      
      tabela <- table(
        periodo,
        evolucao
      )
      
      if(
        nrow(tabela)==2 &
        ncol(tabela)>=2
      ){
        
        suppressWarnings(
          chisq.test(
            tabela,
            correct = FALSE
          )$p.value
        )
        
      } else {
        NA_real_
      }
    },
    
    .groups="drop"
  ) %>%
  
  mutate(
    signif = ifelse(
      p < 0.05,
      "*",
      ""
    )
  )

#=========================================================
# ASTERISCOS
#=========================================================

asteriscos_regiao <- regiao_periodo %>%
  
  group_by(
    faixa_etaria,
    regiao
  ) %>%
  
  summarise(
    y = max(CFR_high)+0.04,
    .groups="drop"
  ) %>%
  
  left_join(
    teste_regiao,
    by=c(
      "faixa_etaria",
      "regiao"
    )
  ) %>%
  
  mutate(
    
    faixa_num = case_when(
      faixa_etaria=="18-40" ~ 1,
      faixa_etaria=="41-60" ~ 2,
      faixa_etaria=="61-80" ~ 3,
      faixa_etaria==">80" ~ 4
    ),
    
    deslocamento = case_when(
      regiao=="North" ~ -0.38,
      regiao=="Northeast" ~ -0.19,
      regiao=="Central-West" ~ 0,
      regiao=="Southeast" ~ 0.19,
      regiao=="South" ~ 0.38
    ),
    
    x = faixa_num + deslocamento
  )

#=========================================================
# GRAFICO
#=========================================================

grafico_regiao_periodo <- ggplot(
  regiao_periodo,
  aes(
    x = faixa_etaria,
    y = CFR,
    fill = grupo
  )
) +
  
  geom_col(
    position = position_dodge(width = 0.95),
    color = "black",
    width = 0.8
  ) +
  
  geom_errorbar(
    aes(
      ymin = CFR_low,
      ymax = CFR_high
    ),
    position = position_dodge(width = 0.95),
    width = 0.2,
    linewidth = 0.5
  ) +
  
  geom_text(
    data = asteriscos_regiao,
    aes(
      x = x,
      y = y,
      label = signif
    ),
    inherit.aes = FALSE,
    size = 6.5,
    family = "Times New Roman",
    fontface = "bold"
  ) +
  
  scale_fill_manual(
    values = c(
      "North - Pre"="#c7e9c0",
      "North - Post"="#238b45",
      "Northeast - Pre"="#dadaeb",
      "Northeast - Post"="#756bb1",
      "Central-West - Pre"="#fdd0a2",
      "Central-West - Post"="#e6550d",
      "Southeast - Pre"="#fcbfd2",
      "Southeast - Post"="#dd1c77",
      "South - Pre"="#bdd7e7",
      "South - Post"="#3182bd"
    ),
    
    labels = c(
      "North, pre-pandemic",
      "North, post-pandemic",
      "Northeast, pre-pandemic",
      "Northeast, post-pandemic",
      "Central-West, pre-pandemic",
      "Central-West, post-pandemic",
      "Southeast, pre-pandemic",
      "Southeast, post-pandemic",
      "South, pre-pandemic",
      "South, post-pandemic"
    ),
    
    name=NULL
  ) +
  
  scale_y_continuous(
    limits=c(0,1),
    labels=scales::percent_format(
      accuracy=1
    )
  ) +
  
  labs(
    title =
      "In-Hospital Case Fatality Rates according to region and age group",
    
    subtitle =
      "Adult ICU patients hospitalized with influenza-associated SARI in Brazil",
    
    x =
      "Age group (years)",
    
    y =
      "Case fatality rate (%)",
    
    caption =
      paste0(
        "Bars represent crude in-hospital case fatality rates and error bars represent 95% confidence intervals estimated using the Wilson method.\n",
        "Pre-pandemic period includes 2013–2019 and post-pandemic period includes 2023–2025.\n",
        "* p<0.05 for comparisons between pre- and post-pandemic periods within each region and age stratum using Pearson's chi-square test."
      )
  ) +
  
  theme_classic() +
  
  theme(
    panel.background=element_rect(fill="white",color=NA),
    plot.background=element_rect(fill="white",color=NA),
    panel.grid.major.y=element_line(color="grey85",linewidth=0.8),
    panel.grid.minor.y=element_blank(),
    axis.line=element_line(color="black",linewidth=0.8),
    axis.ticks=element_line(color="black"),
    axis.text=element_text(size=12,family="Times New Roman"),
    axis.title.y=element_text(size=12,family="Times New Roman"),
    axis.title.x=element_text(size=12,family="Times New Roman"),
    plot.title=element_text(size=14,face="bold",family="Times New Roman"),
    plot.subtitle=element_text(size=11,family="Times New Roman"),
    plot.caption=element_text(size=12,family="Times New Roman",hjust=0),
    legend.position="top",
    legend.justification="left",
    legend.text=element_text(size=11,family="Times New Roman"),
    legend.title=element_blank()
  )

grafico_regiao_periodo

ggsave(
  filename="Figure_CFR_Region_AgeGroup.png",
  plot=grafico_regiao_periodo,
  width=14,
  height=8,
  units="in",
  dpi=600,
  bg="white"
)










#=========================================================
# ESCOLARIDADE
#=========================================================

escolaridade_periodo <- influ_filter %>% 
  
  filter(
    !is.na(escolaridade),
    escolaridade != "",
    !is.na(evolucao),
    !is.na(periodo)
  ) %>% 
  
  mutate(
    
    escolaridade = str_trim(escolaridade),
    
    escolaridade = case_when(
      escolaridade == "Illiterate" ~ "Illiterate",
      escolaridade == "Up to high school" ~ "Up to high school",
      escolaridade == "High school" ~ "High school",
      escolaridade == "College/University" ~ "College/University",
      TRUE ~ "Unknown"
    )
  ) %>% 
  
  group_by(
    periodo,
    escolaridade,
    faixa_etaria
  ) %>% 
  
  summarise(
    n = n(),
    morte = sum(evolucao == "Death"),
    CFR = ifelse(
      n > 0,
      morte / n,
      NA_real_
    ),
    .groups = "drop"
  )

ic_escolaridade <- binom.confint(
  x = escolaridade_periodo$morte,
  n = escolaridade_periodo$n,
  method = "wilson"
)

escolaridade_periodo <- escolaridade_periodo %>% 
  
  mutate(
    
    CFR_low = ic_escolaridade$lower,
    CFR_high = ic_escolaridade$upper,
    
    periodo = factor(
      periodo,
      levels = c("pre","pos"),
      labels = c("Pre","Post")
    ),
    
    escolaridade = factor(
      escolaridade,
      levels = c(
        "Illiterate",
        "Up to high school",
        "High school",
        "College/University"
      )
    ),
    
    grupo = interaction(
      escolaridade,
      periodo,
      sep = " - "
    ),
    
    grupo = factor(
      grupo,
      levels = c(
        "Illiterate - Pre",
        "Illiterate - Post",
        "Up to high school - Pre",
        "Up to high school - Post",
        "High school - Pre",
        "High school - Post",
        "College/University - Pre",
        "College/University - Post"
      )
    )
  )

#=========================================================
# TESTE CHI-SQUARE
#=========================================================

teste_escolaridade <- influ_filter %>%
  
  filter(
    !is.na(escolaridade),
    escolaridade != "",
    !is.na(evolucao),
    !is.na(periodo)
  ) %>%
  
  mutate(
    
    escolaridade = str_trim(escolaridade),
    
    escolaridade = case_when(
      escolaridade == "Illiterate" ~ "Illiterate",
      escolaridade == "Up to high school" ~ "Up to high school",
      escolaridade == "High school" ~ "High school",
      escolaridade == "College/University" ~ "College/University",
      TRUE ~ "Unknown"
    )
  ) %>%
  
  group_by(
    faixa_etaria,
    escolaridade
  ) %>%
  
  summarise(
    
    p = {
      
      tabela <- table(
        periodo,
        evolucao
      )
      
      if(
        nrow(tabela)==2 &
        ncol(tabela)>=2
      ){
        
        suppressWarnings(
          chisq.test(
            tabela,
            correct = FALSE
          )$p.value
        )
        
      } else {
        NA_real_
      }
    },
    
    .groups="drop"
  ) %>%
  
  mutate(
    signif = ifelse(
      p < 0.05,
      "*",
      ""
    )
  )

#=========================================================
# ASTERISCOS
#=========================================================

asteriscos_escolaridade <- escolaridade_periodo %>%
  
  group_by(
    faixa_etaria,
    escolaridade
  ) %>%
  
  summarise(
    y = max(CFR_high)+0.04,
    .groups="drop"
  ) %>%
  
  left_join(
    teste_escolaridade,
    by=c(
      "faixa_etaria",
      "escolaridade"
    )
  ) %>%
  
  mutate(
    
    faixa_num = case_when(
      faixa_etaria=="18-40" ~ 1,
      faixa_etaria=="41-60" ~ 2,
      faixa_etaria=="61-80" ~ 3,
      faixa_etaria==">80" ~ 4
    ),
    
    deslocamento = case_when(
      escolaridade=="Illiterate" ~ -0.30,
      escolaridade=="Up to high school" ~ -0.10,
      escolaridade=="High school" ~ 0.10,
      escolaridade=="College/University" ~ 0.30
    ),
    
    x = faixa_num + deslocamento
  )

#=========================================================
# GRAFICO
#=========================================================

grafico_escolaridade_periodo <- ggplot(
  escolaridade_periodo,
  aes(
    x = faixa_etaria,
    y = CFR,
    fill = grupo
  )
) +
  
  geom_col(
    position = position_dodge(width = 0.95),
    color = "black",
    width = 0.8
  ) +
  
  geom_errorbar(
    aes(
      ymin = CFR_low,
      ymax = CFR_high
    ),
    position = position_dodge(width = 0.95),
    width = 0.2,
    linewidth = 0.5
  ) +
  
  geom_text(
    data = asteriscos_escolaridade,
    aes(
      x = x,
      y = y,
      label = signif
    ),
    inherit.aes = FALSE,
    size = 6.5,
    family = "Times New Roman",
    fontface = "bold"
  ) +
  
  scale_fill_manual(
    values = c(
      "Illiterate - Pre"="#d9d9d9",
      "Illiterate - Post"="#7f7f7f",
      "Up to high school - Pre"="#bdd7e7",
      "Up to high school - Post"="#3182bd",
      "High school - Pre"="#fcbfd2",
      "High school - Post"="#dd1c77",
      "College/University - Pre"="#c7e9c0",
      "College/University - Post"="#238b45"
    ),
    
    labels = c(
      "Illiterate, pre-pandemic",
      "Illiterate, post-pandemic",
      "Up to high school, pre-pandemic",
      "Up to high school, post-pandemic",
      "High school, pre-pandemic",
      "High school, post-pandemic",
      "College/University, pre-pandemic",
      "College/University, post-pandemic"
    ),
    
    name=NULL
  ) +
  
  scale_y_continuous(
    limits=c(0,1),
    labels=scales::percent_format(
      accuracy=1
    )
  ) +
  
  labs(
    title =
      "In-Hospital Case Fatality Rates according to educational level and age group",
    
    subtitle =
      "Adult ICU patients hospitalized with influenza-associated SARI in Brazil",
    
    x =
      "Age group (years)",
    
    y =
      "Case fatality rate (%)",
    
    caption =
      paste0(
        "Bars represent crude in-hospital case fatality rates and error bars represent 95% confidence intervals estimated using the Wilson method.\n",
        "Pre-pandemic period includes 2013–2019 and post-pandemic period includes 2023–2025.\n",
        "* p<0.05 for comparisons between pre- and post-pandemic periods within each educational level and age stratum using Pearson's chi-square test."
      )
  ) +
  
  theme_classic() +
  
  theme(
    panel.background=element_rect(fill="white",color=NA),
    plot.background=element_rect(fill="white",color=NA),
    panel.grid.major.y=element_line(color="grey85",linewidth=0.8),
    panel.grid.minor.y=element_blank(),
    axis.line=element_line(color="black",linewidth=0.8),
    axis.ticks=element_line(color="black"),
    axis.text=element_text(size=12,family="Times New Roman"),
    axis.title.y=element_text(size=12,family="Times New Roman"),
    axis.title.x=element_text(size=12,family="Times New Roman"),
    plot.title=element_text(size=14,face="bold",family="Times New Roman"),
    plot.subtitle=element_text(size=11,family="Times New Roman"),
    plot.caption=element_text(size=12,family="Times New Roman",hjust=0),
    legend.position="top",
    legend.justification="left",
    legend.text=element_text(size=10,family="Times New Roman"),
    legend.title=element_blank()
  )

grafico_escolaridade_periodo

ggsave(
  filename="Figure_CFR_EducationalLevel_AgeGroup.png",
  plot=grafico_escolaridade_periodo,
  width=15,
  height=8,
  units="in",
  dpi=600,
  bg="white"
)











#=========================================================
# CFR RACE
#=========================================================

raca_periodo <- influ_filter %>% 
  
  filter(
    !is.na(raca),
    raca != "",
    !is.na(evolucao),
    !is.na(periodo)
  ) %>% 
  
  mutate(
    
    raca = str_trim(raca),
    
    raca = case_when(
      raca == "White" ~ "White",
      raca == "Black/Brown" ~ "Black/Brown",
      raca == "Asian" ~ "Asian",
      raca == "Indigenous" ~ "Indigenous",
      TRUE ~ "Unknown"
    )
  ) %>% 
  
  group_by(
    periodo,
    raca,
    faixa_etaria
  ) %>% 
  
  summarise(
    n = n(),
    morte = sum(evolucao == "Death"),
    CFR = ifelse(
      n > 0,
      morte / n,
      NA_real_
    ),
    .groups = "drop"
  )

ic_raca <- binom.confint(
  x = raca_periodo$morte,
  n = raca_periodo$n,
  method = "wilson"
)

raca_periodo <- raca_periodo %>% 
  
  mutate(
    
    CFR_low = ic_raca$lower,
    CFR_high = ic_raca$upper,
    
    periodo = factor(
      periodo,
      levels = c("pre","pos"),
      labels = c("Pre","Post")
    ),
    
    raca = factor(
      raca,
      levels = c(
        "White",
        "Black/Brown",
        "Asian",
        "Indigenous"
      )
    ),
    
    grupo = interaction(
      raca,
      periodo,
      sep = " - "
    ),
    
    grupo = factor(
      grupo,
      levels = c(
        "White - Pre",
        "White - Post",
        "Black/Brown - Pre",
        "Black/Brown - Post",
        "Asian - Pre",
        "Asian - Post",
        "Indigenous - Pre",
        "Indigenous - Post"
      )
    )
  )

#=========================================================
# TESTE CHI-SQUARE
#=========================================================

teste_raca <- influ_filter %>%
  
  filter(
    !is.na(raca),
    raca != "",
    !is.na(evolucao),
    !is.na(periodo)
  ) %>%
  
  mutate(
    
    raca = str_trim(raca),
    
    raca = case_when(
      raca == "White" ~ "White",
      raca == "Black/Brown" ~ "Black/Brown",
      raca == "Asian" ~ "Asian",
      raca == "Indigenous" ~ "Indigenous",
      TRUE ~ "Unknown"
    )
  ) %>%
  
  group_by(
    faixa_etaria,
    raca
  ) %>%
  
  summarise(
    
    p = {
      
      tabela <- table(
        periodo,
        evolucao
      )
      
      if(
        nrow(tabela)==2 &
        ncol(tabela)>=2
      ){
        
        suppressWarnings(
          chisq.test(
            tabela,
            correct = FALSE
          )$p.value
        )
        
      } else {
        NA_real_
      }
    },
    
    .groups="drop"
  ) %>%
  
  mutate(
    signif = ifelse(
      p < 0.05,
      "*",
      ""
    )
  )

#=========================================================
# ASTERISCOS
#=========================================================

asteriscos_raca <- raca_periodo %>%
  
  group_by(
    faixa_etaria,
    raca
  ) %>%
  
  summarise(
    y = max(CFR_high)+0.04,
    .groups="drop"
  ) %>%
  
  left_join(
    teste_raca,
    by=c(
      "faixa_etaria",
      "raca"
    )
  ) %>%
  
  mutate(
    
    faixa_num = case_when(
      faixa_etaria=="18-40" ~ 1,
      faixa_etaria=="41-60" ~ 2,
      faixa_etaria=="61-80" ~ 3,
      faixa_etaria==">80" ~ 4
    ),
    
    deslocamento = case_when(
      raca=="White" ~ -0.30,
      raca=="Black/Brown" ~ -0.10,
      raca=="Asian" ~ 0.10,
      raca=="Indigenous" ~ 0.30
    ),
    
    x = faixa_num + deslocamento
  )

#=========================================================
# GRAFICO
#=========================================================

grafico_raca_periodo <- ggplot(
  raca_periodo,
  aes(
    x = faixa_etaria,
    y = CFR,
    fill = grupo
  )
) +
  
  geom_col(
    position = position_dodge(width = 0.95),
    color = "black",
    width = 0.8
  ) +
  
  geom_errorbar(
    aes(
      ymin = CFR_low,
      ymax = CFR_high
    ),
    position = position_dodge(width = 0.95),
    width = 0.2,
    linewidth = 0.5
  ) +
  
  geom_text(
    data = asteriscos_raca,
    aes(
      x = x,
      y = y,
      label = signif
    ),
    inherit.aes = FALSE,
    size = 6.5,
    family = "Times New Roman",
    fontface = "bold"
  ) +
  
  scale_fill_manual(
    values = c(
      "White - Pre"="#bdd7e7",
      "White - Post"="#3182bd",
      "Black/Brown - Pre"="#fcbfd2",
      "Black/Brown - Post"="#dd1c77",
      "Asian - Pre"="#c7e9c0",
      "Asian - Post"="#238b45",
      "Indigenous - Pre"="#fdd0a2",
      "Indigenous - Post"="#e6550d"
    ),
    
    labels = c(
      "White, pre-pandemic",
      "White, post-pandemic",
      "Black/Brown, pre-pandemic",
      "Black/Brown, post-pandemic",
      "Asian, pre-pandemic",
      "Asian, post-pandemic",
      "Indigenous, pre-pandemic",
      "Indigenous, post-pandemic"
    ),
    
    name = NULL
  ) +
  
  scale_y_continuous(
    limits = c(0,1),
    labels = scales::percent_format(
      accuracy = 1
    )
  ) +
  
  labs(
    title =
      "In-Hospital Case Fatality Rates according to race and age group",
    
    subtitle =
      "Adult ICU patients hospitalized with influenza-associated SARI in Brazil",
    
    x =
      "Age group (years)",
    
    y =
      "Case fatality rate (%)",
    
    caption =
      paste0(
        "Bars represent crude in-hospital case fatality rates and error bars represent 95% confidence intervals estimated using the Wilson method.\n",
        "Pre-pandemic period includes 2013–2019 and post-pandemic period includes 2023–2025.\n",
        "* p<0.05 for comparisons between pre- and post-pandemic periods within each race category and age stratum using Pearson's chi-square test."
      )
  ) +
  
  theme_classic() +
  
  theme(
    panel.background = element_rect(fill="white",color=NA),
    plot.background = element_rect(fill="white",color=NA),
    
    panel.grid.major.y =
      element_line(
        color="grey85",
        linewidth=0.8
      ),
    
    panel.grid.minor.y =
      element_blank(),
    
    axis.line =
      element_line(
        color="black",
        linewidth=0.8
      ),
    
    axis.ticks =
      element_line(
        color="black"
      ),
    
    axis.text =
      element_text(
        size=12,
        family="Times New Roman"
      ),
    
    axis.title.y =
      element_text(
        size=12,
        family="Times New Roman"
      ),
    
    axis.title.x =
      element_text(
        size=12,
        family="Times New Roman"
      ),
    
    plot.title =
      element_text(
        size=14,
        face="bold",
        family="Times New Roman"
      ),
    
    plot.subtitle =
      element_text(
        size=11,
        family="Times New Roman"
      ),
    
    plot.caption =
      element_text(
        size=12,
        family="Times New Roman",
        hjust=0
      ),
    
    legend.position="top",
    legend.justification="left",
    legend.text=
      element_text(
        size=10,
        family="Times New Roman"
      ),
    
    legend.title=
      element_blank()
  )

grafico_raca_periodo

ggsave(
  filename = "Figure_CFR_Race_AgeGroup.png",
  plot = grafico_raca_periodo,
  width = 15,
  height = 8,
  units = "in",
  dpi = 600,
  bg = "white"
)


#
#
#
#
#
#
#

#=========================================================
# CFR COMORBIDITY × VENTILATORY SUPPORT
#=========================================================

vm_comorb_periodo <- influ_filter %>% 
  
  filter(
    !is.na(comorbidade_cat),
    !is.na(suporte_vent),
    !is.na(evolucao),
    !is.na(periodo)
  ) %>% 
  
  filter(
    suporte_vent %in% c(
      "No",
      "Yes, non-invasive",
      "Yes, invasive"
    )
  ) %>% 
  
  group_by(
    periodo,
    comorbidade_cat,
    suporte_vent
  ) %>% 
  
  summarise(
    n = n(),
    morte = sum(evolucao=="Death"),
    CFR = ifelse(
      n>0,
      morte/n,
      NA_real_
    ),
    .groups="drop"
  )

ic_vm_comorb <- binom.confint(
  x=vm_comorb_periodo$morte,
  n=vm_comorb_periodo$n,
  method="wilson"
)

vm_comorb_periodo <- vm_comorb_periodo %>% 
  
  mutate(
    
    CFR_low = ic_vm_comorb$lower,
    CFR_high = ic_vm_comorb$upper,
    
    comorbidade_cat = factor(
      comorbidade_cat,
      levels=c(
        "0",
        "1",
        ">1"
      )
    ),
    
    periodo = factor(
      periodo,
      levels=c("pre","pos"),
      labels=c("Pre","Post")
    ),
    
    grupo = interaction(
      suporte_vent,
      periodo,
      sep=" - "
    ),
    
    grupo = factor(
      grupo,
      levels=c(
        "No - Pre",
        "No - Post",
        "Yes, non-invasive - Pre",
        "Yes, non-invasive - Post",
        "Yes, invasive - Pre",
        "Yes, invasive - Post"
      )
    )
  )

#=========================================================
# TESTE CHI-SQUARE
#=========================================================

teste_vm_comorb <- influ_filter %>%
  
  filter(
    !is.na(comorbidade_cat),
    !is.na(suporte_vent),
    !is.na(evolucao),
    !is.na(periodo)
  ) %>%
  
  filter(
    suporte_vent %in% c(
      "No",
      "Yes, non-invasive",
      "Yes, invasive"
    )
  ) %>%
  
  group_by(
    comorbidade_cat,
    suporte_vent
  ) %>%
  
  summarise(
    
    p={
      
      tabela <- table(
        periodo,
        evolucao
      )
      
      if(
        nrow(tabela)==2 &
        ncol(tabela)>=2
      ){
        
        suppressWarnings(
          chisq.test(
            tabela,
            correct=FALSE
          )$p.value
        )
        
      } else {
        NA_real_
      }
    },
    
    .groups="drop"
  ) %>%
  
  mutate(
    signif=ifelse(
      p<0.05,
      "*",
      ""
    )
  )

#=========================================================
# ASTERISCOS
#=========================================================

asteriscos_vm_comorb <- vm_comorb_periodo %>%
  
  group_by(
    comorbidade_cat,
    suporte_vent
  ) %>%
  
  summarise(
    y=max(CFR_high)+0.04,
    .groups="drop"
  ) %>%
  
  left_join(
    teste_vm_comorb,
    by=c(
      "comorbidade_cat",
      "suporte_vent"
    )
  ) %>%
  
  mutate(
    
    x_num = case_when(
      comorbidade_cat=="0" ~ 1,
      comorbidade_cat=="1" ~ 2,
      comorbidade_cat==">1" ~ 3
    ),
    
    deslocamento = case_when(
      suporte_vent=="No" ~ -0.28,
      suporte_vent=="Yes, non-invasive" ~ 0,
      suporte_vent=="Yes, invasive" ~ 0.28
    ),
    
    x = x_num + deslocamento
  )

#=========================================================
# GRAFICO
#=========================================================

grafico_vm_comorb <- ggplot(
  vm_comorb_periodo,
  aes(
    x=comorbidade_cat,
    y=CFR,
    fill=grupo
  )
) +
  
  geom_col(
    position=position_dodge(width=0.9),
    color="black",
    width=0.8
  ) +
  
  geom_errorbar(
    aes(
      ymin=CFR_low,
      ymax=CFR_high
    ),
    position=position_dodge(width=0.9),
    width=0.2,
    linewidth=0.5
  ) +
  
  geom_text(
    data=asteriscos_vm_comorb,
    aes(
      x=x,
      y=y,
      label=signif
    ),
    inherit.aes=FALSE,
    size=6.5,
    family="Times New Roman",
    fontface="bold"
  ) +
  
  scale_fill_manual(
    values=c(
      "No - Pre"="#aecde1",
      "No - Post"="#2171b5",
      "Yes, non-invasive - Pre"="#c7e9c0",
      "Yes, non-invasive - Post"="#238b45",
      "Yes, invasive - Pre"="#f5b689",
      "Yes, invasive - Post"="#eb690c"
    ),
   
    labels=c(
      "No ventilatory support, pre-pandemic",
      "No ventilatory support, post-pandemic",
      "Non-invasive ventilation, pre-pandemic",
      "Non-invasive ventilation, post-pandemic",
      "Invasive mechanical ventilation, pre-pandemic",
      "Invasive mechanical ventilation, post-pandemic"
    ),
    
    name=NULL
  ) +
  
  scale_y_continuous(
    limits=c(0,1),
    labels=scales::percent_format(
      accuracy=1
    )
  ) +
  
  labs(
    title=
      "In-Hospital Case Fatality Rates according to comorbidity burden and ventilatory support",
    
    subtitle=
      "Adult ICU patients hospitalized with influenza-associated SARI in Brazil",
    
    x=
      "Number of comorbidities",
    
    y=
      "Case fatality rate (%)",
    
    caption=
      paste0(
        "Bars represent crude in-hospital case fatality rates and error bars represent 95% confidence intervals estimated using the Wilson method.\n",
        "Pre-pandemic period includes 2013–2019 and post-pandemic period includes 2023–2025.\n",
        "* p<0.05 for comparisons between pre- and post-pandemic periods within each ventilatory support category and comorbidity stratum using Pearson's chi-square test."
      )
  ) +
  
  theme_classic() +
  
  theme(
    panel.background=element_rect(fill="white",color=NA),
    plot.background=element_rect(fill="white",color=NA),
    panel.grid.major.y=element_line(color="grey85",linewidth=0.8),
    panel.grid.minor.y=element_blank(),
    axis.line=element_line(color="black",linewidth=0.8),
    axis.ticks=element_line(color="black"),
    axis.text=element_text(size=12,family="Times New Roman"),
    axis.title=element_text(size=12,family="Times New Roman"),
    plot.title=element_text(size=14,face="bold",family="Times New Roman"),
    plot.subtitle=element_text(size=11,family="Times New Roman"),
    plot.caption=element_text(size=12,family="Times New Roman",hjust=0),
    legend.position="top",
    legend.justification="left",
    legend.text=element_text(size=11,family="Times New Roman"),
    legend.title=element_blank()
  )

grafico_vm_comorb

ggsave(
  filename="Figure_CFR_Comorbidity_VM.png",
  plot=grafico_vm_comorb,
  width=14,
  height=8,
  units="in",
  dpi=600,
  bg="white"
)







#=========================================================
# CFR REGION × VENTILATORY SUPPORT
#=========================================================

vm_regiao_periodo <- influ_filter %>% 
  
  filter(
    !is.na(regiao),
    !is.na(suporte_vent),
    !is.na(evolucao),
    !is.na(periodo)
  ) %>% 
  
  mutate(
    regiao = case_when(
      regiao=="Norte" ~ "North",
      regiao=="Nordeste" ~ "Northeast",
      regiao=="Centro-Oeste" ~ "Central-West",
      regiao=="Sudeste" ~ "Southeast",
      regiao=="Sul" ~ "South",
      TRUE ~ NA_character_
    )
  ) %>% 
  
  filter(
    suporte_vent %in% c(
      "No",
      "Yes, non-invasive",
      "Yes, invasive"
    )
  ) %>% 
  
  group_by(
    periodo,
    regiao,
    suporte_vent
  ) %>% 
  
  summarise(
    n = n(),
    morte = sum(evolucao=="Death"),
    CFR = ifelse(
      n>0,
      morte/n,
      NA_real_
    ),
    .groups="drop"
  )

ic_vm_regiao <- binom.confint(
  x=vm_regiao_periodo$morte,
  n=vm_regiao_periodo$n,
  method="wilson"
)

vm_regiao_periodo <- vm_regiao_periodo %>% 
  
  mutate(
    
    CFR_low = ic_vm_regiao$lower,
    CFR_high = ic_vm_regiao$upper,
    
    regiao = factor(
      regiao,
      levels=c(
        "North",
        "Northeast",
        "Central-West",
        "Southeast",
        "South"
      )
    ),
    
    periodo = factor(
      periodo,
      levels=c("pre","pos"),
      labels=c("Pre","Post")
    ),
    
    grupo = interaction(
      suporte_vent,
      periodo,
      sep=" - "
    ),
    
    grupo = factor(
      grupo,
      levels=c(
        "No - Pre",
        "No - Post",
        "Yes, non-invasive - Pre",
        "Yes, non-invasive - Post",
        "Yes, invasive - Pre",
        "Yes, invasive - Post"
      )
    )
  )

#=========================================================
# TESTE CHI-SQUARE
#=========================================================

teste_vm_regiao <- influ_filter %>%
  
  filter(
    !is.na(regiao),
    !is.na(suporte_vent),
    !is.na(evolucao),
    !is.na(periodo)
  ) %>%
  
  mutate(
    regiao = case_when(
      regiao=="Norte" ~ "North",
      regiao=="Nordeste" ~ "Northeast",
      regiao=="Centro-Oeste" ~ "Central-West",
      regiao=="Sudeste" ~ "Southeast",
      regiao=="Sul" ~ "South",
      TRUE ~ NA_character_
    )
  ) %>%
  
  filter(
    suporte_vent %in% c(
      "No",
      "Yes, non-invasive",
      "Yes, invasive"
    )
  ) %>%
  
  group_by(
    regiao,
    suporte_vent
  ) %>%
  
  summarise(
    
    p={
      
      tabela <- table(
        periodo,
        evolucao
      )
      
      if(
        nrow(tabela)==2 &
        ncol(tabela)>=2
      ){
        
        suppressWarnings(
          chisq.test(
            tabela,
            correct=FALSE
          )$p.value
        )
        
      } else {
        NA_real_
      }
    },
    
    .groups="drop"
  ) %>%
  
  mutate(
    signif=ifelse(
      p<0.05,
      "*",
      ""
    )
  )

#=========================================================
# ASTERISCOS
#=========================================================

asteriscos_vm_regiao <- vm_regiao_periodo %>%
  
  group_by(
    regiao,
    suporte_vent
  ) %>%
  
  summarise(
    y=max(CFR_high)+0.04,
    .groups="drop"
  ) %>%
  
  left_join(
    teste_vm_regiao,
    by=c(
      "regiao",
      "suporte_vent"
    )
  ) %>%
  
  mutate(
    
    x_num = case_when(
      regiao=="North" ~ 1,
      regiao=="Northeast" ~ 2,
      regiao=="Central-West" ~ 3,
      regiao=="Southeast" ~ 4,
      regiao=="South" ~ 5
    ),
    
    deslocamento = case_when(
      suporte_vent=="No" ~ -0.28,
      suporte_vent=="Yes, non-invasive" ~ 0,
      suporte_vent=="Yes, invasive" ~ 0.28
    ),
    
    x = x_num + deslocamento
  )

#=========================================================
# GRAFICO
#=========================================================

grafico_vm_regiao <- ggplot(
  vm_regiao_periodo,
  aes(
    x=regiao,
    y=CFR,
    fill=grupo
  )
) +
  
  geom_col(
    position=position_dodge(width=0.9),
    color="black",
    width=0.8
  ) +
  
  geom_errorbar(
    aes(
      ymin=CFR_low,
      ymax=CFR_high
    ),
    position=position_dodge(width=0.9),
    width=0.2,
    linewidth=0.5
  ) +
  
  geom_text(
    data=asteriscos_vm_regiao,
    aes(
      x=x,
      y=y,
      label=signif
    ),
    inherit.aes=FALSE,
    size=6.5,
    family="Times New Roman",
    fontface="bold"
  ) +
  
  scale_fill_manual(
    values=c(
      "No - Pre"="#aecde1",
      "No - Post"="#2171b5",
      "Yes, non-invasive - Pre"="#c7e9c0",
      "Yes, non-invasive - Post"="#238b45",
      "Yes, invasive - Pre"="#f5b689",
      "Yes, invasive - Post"="#eb690c"
    ),
    
    labels=c(
      "No ventilatory support, pre-pandemic",
      "No ventilatory support, post-pandemic",
      "Non-invasive ventilation, pre-pandemic",
      "Non-invasive ventilation, post-pandemic",
      "Invasive mechanical ventilation, pre-pandemic",
      "Invasive mechanical ventilation, post-pandemic"
    ),
    
    name=NULL
  ) +
  
  scale_y_continuous(
    limits=c(0,1),
    labels=scales::percent_format(
      accuracy=1
    )
  ) +
  
  labs(
    title=
      "In-Hospital Case Fatality Rates according to region and ventilatory support",
    
    subtitle=
      "Adult ICU patients hospitalized with influenza-associated SARI in Brazil",
    
    x=
      "Region",
    
    y=
      "Case fatality rate (%)",
    
    caption=
      paste0(
        "Bars represent crude in-hospital case fatality rates and error bars represent 95% confidence intervals estimated using the Wilson method.\n",
        "Pre-pandemic period includes 2013–2019 and post-pandemic period includes 2023–2025.\n",
        "* p<0.05 for comparisons between pre- and post-pandemic periods within each ventilatory support category and region using Pearson's chi-square test."
      )
  ) +
  
  theme_classic() +
  
  theme(
    panel.background=element_rect(fill="white",color=NA),
    plot.background=element_rect(fill="white",color=NA),
    panel.grid.major.y=element_line(color="grey85",linewidth=0.8),
    panel.grid.minor.y=element_blank(),
    axis.line=element_line(color="black",linewidth=0.8),
    axis.ticks=element_line(color="black"),
    axis.text=element_text(size=12,family="Times New Roman"),
    axis.title=element_text(size=12,family="Times New Roman"),
    plot.title=element_text(size=14,face="bold",family="Times New Roman"),
    plot.subtitle=element_text(size=11,family="Times New Roman"),
    plot.caption=element_text(size=12,family="Times New Roman",hjust=0),
    legend.position="top",
    legend.justification="left",
    legend.text=element_text(size=11,family="Times New Roman"),
    legend.title=element_blank()
  )

grafico_vm_regiao

ggsave(
  filename="Figure_CFR_Region_VM.png",
  plot=grafico_vm_regiao,
  width=15,
  height=8,
  units="in",
  dpi=600,
  bg="white"
)
