#Load libraries
library(dplyr)
library(targets)
library(ggplot2)
library(tidyr)
library(readxl)
library(stringr)
library(here)
library(readr)

#Load data
drugs_metabolites <- read_excel(
  here("Input", "Kreuztabelle Chemie Mittelwerte Arzneistoffe und Metaboliten 2022-2025.xlsx"))
EC10_Fish_drugs_metabolites_GRO <- read_csv(
  here("Input", "TRIDENT_Fish_EC10_GRO_Arzneistoffe und Metaboliten.csv"))
Substance_names <- read_delim(
  here("Input", "ELWAS Arzneistoffe und Metaboliten Stoffe.csv"),
  delim = ";",
  col_names = FALSE)

#Clean data
View(drugs_metabolites)
new_col_names_drugs <-  as.character(unlist(drugs_metabolites[3,]))
names(drugs_metabolites) <-  new_col_names_drugs
drugs_metabolites <-  drugs_metabolites [-c(1:3),]
drugs_metabolites <-  drugs_metabolites [,-c(9)]
colnames(drugs_metabolites)
drugs_metabolites <-  pivot_longer(drugs_metabolites, 
                                   cols =  "Clofibrinsäure [µg/l]" : "Amantadin [µg/l]" , 
                                   names_to = 'Stoffname', 
                                   values_to = 'Stoffkonzentration')
drugs_metabolites$Stoffkonzentration <- as.numeric(gsub(",",".",drugs_metabolites$Stoffkonzentration))
drugs_metabolites$Stoffkonzentration <- drugs_metabolites$Stoffkonzentration * 0.001 ##μg/L in mg/L

#Descriptive Data
number_locations <-length(unique(drugs_metabolites$Messstellenname)) 
number_locations # 817
number_substances <- unique(drugs_metabolites$Stoffname)
number_substances # 126


#Two tables, one Limit of determination  = 0, one Limit of determination = #### muss noch gefunden werden
drugs_metabolites_B0 <- drugs_metabolites
drugs_metabolites_B0[is.na(drugs_metabolites_B0)] <- 0

drugs_metabolites_B0_001 <- drugs_metabolites
drugs_metabolites_B0_001[is.na(drugs_metabolites_B0_001)] <- ### muss noch gefunden werden

#Tables for every year
drugs_metabolites2022_B0 <-  drugs_metabolites_B0 |> filter(Kalenderjahr == 2022)
drugs_metabolites2023_B0 <-  drugs_metabolites_B0 |> filter(Kalenderjahr == 2023)
drugs_metabolites2024_B0 <-  drugs_metabolites_B0 |> filter(Kalenderjahr == 2024)
drugs_metabolites2025_B0 <-  drugs_metabolites_B0 |> filter(Kalenderjahr == 2025)

#Descriptive data Tables per year
table_22 <- table(drugs_metabolites2022_B0$Messstellennummer)
length(unique(names(table_22))) # 345 Messstellen

table_23 <- table(drugs_metabolites2023_B0$Messstellennummer)
length(unique(names(table_23))) # 272 Messstellen

table_24<- table(drugs_metabolites2024_B0$Messstellennummer)
length(unique(names(table_24))) # 332 Messstellen

table_25 <- table(drugs_metabolites2025_B0$Messstellennummer)
length(unique(names(table_25))) # 99 Messstellen

normalize_names <- function(x) {
  trimws(tolower(x))  # entfernt Leerzeichen, setzt auf Kleinschreibung
}

names_22 <- normalize_names(names(table_22))
names_23 <- normalize_names(names(table_23))
names_24 <- normalize_names(names(table_24))
names_25 <- normalize_names(names(table_25))

setdiff(names(table_22), names(table_23)) ## 251 Messstellen, die nicht in Tabelle 23 vorkommen --> 94 gleiche Stellen
setdiff(names(table_23), names(table_22)) ## 178 Messstellen, die nicht in Tabelle 22 vorkommen --> 94 gleiche Stellen 

setdiff(names(table_24), names(table_22)) ## 248 Messstellen, die nicht in Tabelle 22 vorkommen --> 84 gleiche Stellen
setdiff(names(table_22), names(table_24)) ## 261 Messstellen, die nicht in Tabelle 24 vorkommen --> 84 gleiche Stellen

setdiff(names(table_22), names(table_25)) ## 318 Messstellen, die nicht in Tabelle 25 vorkommen --> 27 gleiche Stellen
setdiff(names(table_25), names(table_22)) ## 72 Messstellen, die nicht in Tabelle 22 vorkommen --> 27 gleiche Stellen 

setdiff(names(table_23), names(table_24)) ## 212 Messstellen, die nicht in Tabelle 24 vorkommen --> 60 gleiche Stellen
setdiff(names(table_24), names(table_23)) ## 272 Messstellen, die nicht in Tabelle 23 vorkommen --> 60 gleiche Stellen

setdiff(names(table_23), names(table_25)) ## 255 Messstellen, die nicht in Tabelle 25 vorkommen --> 17 gleiche Stellen
setdiff(names(table_25), names(table_23)) ## 82 Messstellen, die nicht in Tabelle 24 vorkommen --> 17 gleiche Stellen 

setdiff(names(table_24), names(table_25)) ## 303 Messstellen, die nicht in Tabelle 25 vorkommen --> 29 gleiche Stellen
setdiff(names(table_25), names(table_24)) ## 70 Messstellen, die nicht in Tabelle 24 vorkommen --> 29 gleiche Stellen 


#Barplot
drugs_metabolites_top10 <- drugs_metabolites |>
  group_by(Stoffname) |>
  slice_max(Stoffkonzentration, n = 1, with_ties = FALSE) |>
  ungroup() |>
  slice_max(Stoffkonzentration, n = 10)

drugs_metabolites_top10 <- mutate(drugs_metabolites_top10, 
                                  Messstellenname = paste(Messstellenname, Kalenderjahr, sep = " "))

plot_drugs_metabolites_top10 <- ggplot(drugs_metabolites_top10, 
                                       aes(x =  reorder(Stoffname, -Stoffkonzentration), 
                                           y = Stoffkonzentration, 
                                           fill = Messstellenname))+
  geom_bar(stat="identity")+
  labs(title = "Top 10 Substances",
       x = "Substancename",
       y = "Substance concentration [mg/l]",
       fill = "Measuring point")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_drugs_metabolites_top10


#Calculate Risk
#Fish
View(EC10_Fish_drugs_metabolites_GRO)
EC10_Fish_drugs_metabolites_GRO <- EC10_Fish_drugs_metabolites_GRO[ , !(names(EC10_Fish_drugs_metabolites_GRO) %in% c(
  "most.similar.chemical",
  "SMILES_Canonical_RDKit",
  "endpoint",
  "effect",
  "Chemical.Alert",
  "SMILES.Alert"
))]
Substance_names <- Substance_names[,!(names(Substance_names)%in% c(
  "V1",
  "V3"
))]
Substance_names
EC10_Fish_drugs_metabolites_GRO <-mutate(EC10_Fish_drugs_metabolites_GRO, Stoffname = Substance_names)

# 1. Clean Stoffnamen in beiden Tabellen
drugs_metabolites_clean <- drugs_metabolites_B0 %>%
  mutate(Stoffname_clean = str_remove(Stoffname, "\\s*\\[.*\\]") %>% str_trim())

EC10_clean <- EC10_Fish_drugs_metabolites_GRO %>%
  mutate(Stoffname_clean = str_trim(Stoffname)) 

# 2. Join über bereinigten Namen
EC10_Fish_Risk_drugs_metabolites <- drugs_metabolites_clean %>%
  left_join(EC10_clean, by = "Stoffname_clean") %>%
  mutate(Risk = Stoffkonzentration / predictions..mg.L.)
View(EC10)

EC10_Fish_Risk_drugs_metabolites <- EC10_Fish_Risk_drugs_metabolites[,(names(EC10_Fish_Risk_drugs_metabolites)%in% c(
  "Risk",
  "Stoffname_clean",
  "Stoffkonzentration",
  "predictions..mg.L.",
  "Messstellenname",
  "Kalenderjahr"
))]

View(EC10_Fish_Risk_drugs_metabolites)

EC10_Fish_Risk_drugs_metabolites <- EC10_Fish_Risk_drugs_metabolites |>
  group_by(Messstellenname,Kalenderjahr) |>
  mutate(Risk_Mesuringpoint = sum(Risk, na.rm = TRUE)) |>
  ungroup()

EC10_Fish_Risk_drugs_metabolites_topRisk <- EC10_Fish_Risk_drugs_metabolites |>
  group_by(Messstellenname) |>
  filter(Risk_Mesuringpoint>= 1)

EC10_Fish_Risk_drugs_metabolites_topRisk <-  mutate(EC10_Fish_Risk_drugs_metabolites_topRisk, 
                                                    Messstellenname = paste(Messstellenname, Kalenderjahr, sep = " "))

Fish_Risk_drugs_metabolites_plot_topRisk <- ggplot(EC10_Fish_Risk_drugs_metabolites_topRisk,
                                                   aes(x = Messstellenname, y = Risk, fill = Stoffname_clean)) +
  geom_bar(stat = "identity", position = "stack") + 
  labs(title = "Risk per Mesuring point",
       x = "Mesring Point",
       y = "Risk",
       fill = "Substance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "none")


Fish_Risk_drugs_metabolites_plot_topRisk

View(EC10_Fish_Risk_drugs_metabolites_topRisk)
EC10_Fish_Risk_drugs_metabolites_topRisk <- EC10_Fish_Risk_drugs_metabolites_topRisk |>
  arrange(desc(Risk_Mesuringpoint))

top10_Fish_Risk_drugs_metabolites <- EC10_Fish_Risk_drugs_metabolites |>
  arrange(desc(Risk)) |>
  slice_head(n = 10)

Fish_Risk_drugs_metabolites_plot_topRisk_10 <- ggplot(top10_Fish_Risk_drugs_metabolites,
                                                      aes(x = reorder(Stoffname_clean, -Risk),
                                                          y = Risk,
                                                          fill = Messstellenname)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Top 10 Risks",
       x = "Substance",
       y = "Risk",
       fill = "Messuring Point") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

Fish_Risk_drugs_metabolites_plot_topRisk_10
