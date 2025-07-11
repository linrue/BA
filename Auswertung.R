#Load libraries
library(dplyr)
library(targets)
library(ggplot2)
library(tidyr)
library(readxl)
library(stringr)
library(here)
library(readr)

#####################################Load data########################################################################################################################
#Durgs and Metabolites
drugs_metabolites <- read_excel(
  here("Input", "Einzelwerte Chemie_Arzneimittel_Metaboliten_2022_2025.xlsx"))
EC10_Fish_drugs_metabolites_GRO <- read_csv(
  here("Input", "TRIDENT_prediction_Arzneistoffe_fish_EC10_GRO_720h_PubChem.csv"))
EC10_aquatic_invertebrates_drugs_metabolites_POP <- read_csv(
  here("Input", "TRIDENT_prediction_Arzneistoffe_aquatic_invertebrates_EC10_POP_720h_PubChem.csv"))
EC10_algae_drugs_metabolites_POP <- read_csv(
  here("Input", "TRIDENT_prediction_Arzneistoffe_algae_EC10_POP_720h_PubChem.csv"))
Substance_names_drugs <- read_delim(
  here("Input", "ELWAS Arzneistoffe und Metaboliten Stoffe.csv"),
  delim = ";",
  col_names = FALSE)

#PSM and Metabolites
PSM_metabolites_2022 <- read_excel(
  here("Input", "Einzelwerte Chemie_PSM_Metaboliten_2022.xlsx"))
PSM_metabolites_2023 <- read_excel(
  here("Input", "Einzelwerte Chemie_PSM_Metaboliten_2023.xlsx"))
PSM_metabolites_2024 <- read_excel(
  here("Input", "Einzelwerte Chemie_PSM_Metaboliten_2024.xlsx"))
PSM_metabolites_2025 <- read_excel(
  here("Input", "Einzelwerte Chemie_PSM_Metaboliten_2025.xlsx"))
EC10_Fish_PSM_metabolites_GRO <- read_csv(
  here("Input", "TRIDENT_prediction_PSM_Fish_EC10_GRO_SMILES_PubChem_720h.csv"))
EC10_aquatic_invertebrates_PSM_metabolites_POP <- read_csv(
  here("Input", "TRIDENT_prediction_PSM_aquatic invertebrates_EC10_POP_720h_SMILES_PubChem.csv"))
EC10_algae_PSM_metabolites_POP <- read_csv(
  here("Input", "TRIDENT_prediction_PSM_algea_EC10_POP_720h_PubChem.csv"))
Substance_names_PSM <- read_delim(
  here("Input", "ELWAS PSM und Metaboliten Stoffe.csv"),
  delim = ";",
  col_names = FALSE)

####################################################Clean data##################################################################################################################
#Drugs and Metabolites
View(drugs_metabolites)
new_col_names_drugs <-  as.character(unlist(drugs_metabolites[3,]))
names(drugs_metabolites) <-  new_col_names_drugs
drugs_metabolites <-  drugs_metabolites [-c(1:3),]
drugs_metabolites_B0 <- drugs_metabolites |> mutate(Messergebnis = ifelse(grepl("<",Messergebnis),0,Messergebnis))  #Werte an der Bestimmungsgrenze werden 0 gesetzt
View(drugs_metabolites_B0)
drugs_metabolites_BB <- drugs_metabolites |> mutate(Messergebnis = ifelse(grepl("<",Messergebnis),Bestimmungsgrenze,Messergebnis))  #Werte an der Bestimmungsgrenze werden auf die Bestimmungsgrenze gesetzt
View(drugs_metabolites_BB)
drugs_metabolites_B0$Messergebnis <- as.numeric(gsub(",",".",drugs_metabolites_B0$Messergebnis))
drugs_metabolites_B0$Bestimmungsgrenze <- as.numeric(gsub(",",".",drugs_metabolites_B0$Bestimmungsgrenze))
drugs_metabolites_BB$Messergebnis <- as.numeric(gsub(",",".",drugs_metabolites_BB$Messergebnis))
drugs_metabolites_BB$Bestimmungsgrenze <- as.numeric(gsub(",",".",drugs_metabolites_BB$Bestimmungsgrenze))
unique(drugs_metabolites_B0$Einheit)
drugs_metabolites_B0 <- drugs_metabolites_B0 |>   mutate(
  Messergebnis = ifelse(Einheit == "ng/l", Messergebnis / 1000, Messergebnis),
  Einheit = ifelse(Einheit == "ng/l", "µg/l", Einheit)
) #Einheitliche Einheiten
drugs_metabolites_BB <- drugs_metabolites_BB |>   mutate(
  Messergebnis = ifelse(Einheit == "ng/l", Messergebnis / 1000, Messergebnis),
  Einheit = ifelse(Einheit == "ng/l", "µg/l", Einheit)
) #Einheitliche Einheiten


#PSM and Metabolites
#2022
View(PSM_metabolites_2022)
new_col_names_PSM_2022 <-  as.character(unlist(PSM_metabolites_2022[3,]))
names(PSM_metabolites_2022) <-  new_col_names_PSM_2022
PSM_metabolites_2022 <-  PSM_metabolites_2022 [-c(1:3),]
PSM_metabolites_2022_B0 <- PSM_metabolites_2022 |> mutate(Messergebnis = ifelse(grepl("<",Messergebnis),0,Messergebnis))  #Werte an der Bestimmungsgrenze werden 0 gesetzt
View(PSM_metabolites_2022_B0)
PSM_metabolites_2022_BB <- PSM_metabolites_2022 |> mutate(Messergebnis = ifelse(grepl("<",Messergebnis),Bestimmungsgrenze,Messergebnis))  #Werte an der Bestimmungsgrenze werden auf die Bestimmungsgrenze gesetzt
View(PSM_metabolites_2022_BB)
PSM_metabolites_2022_B0$Messergebnis <- as.numeric(gsub(",",".",PSM_metabolites_2022_B0$Messergebnis))
PSM_metabolites_2022_B0$Bestimmungsgrenze <- as.numeric(gsub(",",".",PSM_metabolites_2022_B0$Bestimmungsgrenze))
PSM_metabolites_2022_BB$Messergebnis <- as.numeric(gsub(",",".",PSM_metabolites_2022_BB$Messergebnis))
PSM_metabolites_2022_BB$Bestimmungsgrenze <- as.numeric(gsub(",",".",PSM_metabolites_2022_BB$Bestimmungsgrenze))
unique(PSM_metabolites_2022_B0$Einheit)
PSM_metabolites_2022_B0 <- PSM_metabolites_2022_B0 |>   mutate(
  Messergebnis = ifelse(Einheit == "ng/l", Messergebnis / 1000, Messergebnis),
  Einheit = ifelse(Einheit == "ng/l", "µg/l", Einheit)
)
PSM_metabolites_2022_BB <- PSM_metabolites_2022_BB |>   mutate(
  Messergebnis = ifelse(Einheit == "ng/l", Messergebnis / 1000, Messergebnis),
  Einheit = ifelse(Einheit == "ng/l", "µg/l", Einheit)
)
#2023
View(PSM_metabolites_2023)
new_col_names_PSM_2023 <-  as.character(unlist(PSM_metabolites_2023[3,]))
names(PSM_metabolites_2023) <-  new_col_names_PSM_2023
PSM_metabolites_2023 <-  PSM_metabolites_2023 [-c(1:3),]
PSM_metabolites_2023_B0 <- PSM_metabolites_2023 |> mutate(Messergebnis = ifelse(grepl("<",Messergebnis),0,Messergebnis))  #Werte an der Bestimmungsgrenze werden 0 gesetzt
View(PSM_metabolites_2023_B0)
PSM_metabolites_2023_BB <- PSM_metabolites_2023 |> mutate(Messergebnis = ifelse(grepl("<",Messergebnis),Bestimmungsgrenze,Messergebnis))  #Werte an der Bestimmungsgrenze werden auf die Bestimmungsgrenze gesetzt
View(PSM_metabolites_2023_BB)
PSM_metabolites_2023_B0$Messergebnis <- as.numeric(gsub(",",".",PSM_metabolites_2023_B0$Messergebnis))
PSM_metabolites_2023_B0$Bestimmungsgrenze <- as.numeric(gsub(",",".",PSM_metabolites_2023_B0$Bestimmungsgrenze))
PSM_metabolites_2023_BB$Messergebnis <- as.numeric(gsub(",",".",PSM_metabolites_2023_BB$Messergebnis))
PSM_metabolites_2023_BB$Bestimmungsgrenze <- as.numeric(gsub(",",".",PSM_metabolites_2023_BB$Bestimmungsgrenze))
unique(PSM_metabolites_2023_B0$Einheit)
PSM_metabolites_2023_B0 <- PSM_metabolites_2023_B0 |>   mutate(
  Messergebnis = ifelse(Einheit == "ng/l", Messergebnis / 1000, Messergebnis),
  Einheit = ifelse(Einheit == "ng/l", "µg/l", Einheit)
)
PSM_metabolites_2023_BB <- PSM_metabolites_2023_BB |>   mutate(
  Messergebnis = ifelse(Einheit == "ng/l", Messergebnis / 1000, Messergebnis),
  Einheit = ifelse(Einheit == "ng/l", "µg/l", Einheit)
)
#2024
View(PSM_metabolites_2024)
new_col_names_PSM_2024 <-  as.character(unlist(PSM_metabolites_2024[3,]))
names(PSM_metabolites_2024) <-  new_col_names_PSM_2024
PSM_metabolites_2024 <-  PSM_metabolites_2024 [-c(1:3),]
PSM_metabolites_2024_B0 <- PSM_metabolites_2024 |> mutate(Messergebnis = ifelse(grepl("<",Messergebnis),0,Messergebnis))  #Werte an der Bestimmungsgrenze werden 0 gesetzt
View(PSM_metabolites_2024_B0)
PSM_metabolites_2024_BB <- PSM_metabolites_2024 |> mutate(Messergebnis = ifelse(grepl("<",Messergebnis),Bestimmungsgrenze,Messergebnis))  #Werte an der Bestimmungsgrenze werden auf die Bestimmungsgrenze gesetzt
View(PSM_metabolites_2024_BB)
PSM_metabolites_2024_B0$Messergebnis <- as.numeric(gsub(",",".",PSM_metabolites_2024_B0$Messergebnis))
PSM_metabolites_2024_B0$Bestimmungsgrenze <- as.numeric(gsub(",",".",PSM_metabolites_2024_B0$Bestimmungsgrenze))
PSM_metabolites_2024_BB$Messergebnis <- as.numeric(gsub(",",".",PSM_metabolites_2024_BB$Messergebnis))
PSM_metabolites_2024_BB$Bestimmungsgrenze <- as.numeric(gsub(",",".",PSM_metabolites_2024_BB$Bestimmungsgrenze))
unique(PSM_metabolites_2024_B0$Einheit)
PSM_metabolites_2024_B0 <- PSM_metabolites_2024_B0 |>   mutate(
  Messergebnis = ifelse(Einheit == "ng/l", Messergebnis / 1000, Messergebnis),
  Einheit = ifelse(Einheit == "ng/l", "µg/l", Einheit)
)
PSM_metabolites_2024_BB <- PSM_metabolites_2024_BB |>   mutate(
  Messergebnis = ifelse(Einheit == "ng/l", Messergebnis / 1000, Messergebnis),
  Einheit = ifelse(Einheit == "ng/l", "µg/l", Einheit)
)
#2025
View(PSM_metabolites_2025_BB)
new_col_names_PSM_2025 <-  as.character(unlist(PSM_metabolites_2025[3,]))
names(PSM_metabolites_2025) <-  new_col_names_PSM_2025
PSM_metabolites_2025 <-  PSM_metabolites_2025 [-c(1:3),]
PSM_metabolites_2025_B0 <- PSM_metabolites_2025 |> mutate(Messergebnis = ifelse(grepl("<",Messergebnis),0,Messergebnis))  #Werte an der Bestimmungsgrenze werden 0 gesetzt
View(PSM_metabolites_2025_B0)
PSM_metabolites_2025_BB <- PSM_metabolites_2025 |> mutate(Messergebnis = ifelse(grepl("<",Messergebnis),Bestimmungsgrenze,Messergebnis))  #Werte an der Bestimmungsgrenze werden auf die Bestimmungsgrenze gesetzt
View(PSM_metabolites_2025_BB)
PSM_metabolites_2025_B0$Messergebnis <- as.numeric(gsub(",",".",PSM_metabolites_2025_B0$Messergebnis))
PSM_metabolites_2025_B0$Bestimmungsgrenze <- as.numeric(gsub(",",".",PSM_metabolites_2025_B0$Bestimmungsgrenze))
PSM_metabolites_2025_BB$Messergebnis <- as.numeric(gsub(",",".",PSM_metabolites_2025_BB$Messergebnis))
PSM_metabolites_2025_BB$Bestimmungsgrenze <- as.numeric(gsub(",",".",PSM_metabolites_2025_BB$Bestimmungsgrenze))
unique(PSM_metabolites_2025_B0$Einheit)
PSM_metabolites_2025_B0 <- PSM_metabolites_2025_B0 |>   mutate(
  Messergebnis = ifelse(Einheit == "ng/l", Messergebnis / 1000, Messergebnis),
  Einheit = ifelse(Einheit == "ng/l", "µg/l", Einheit)
)
PSM_metabolites_2025_BB <- PSM_metabolites_2025_BB |>   mutate(
  Messergebnis = ifelse(Einheit == "ng/l", Messergebnis / 1000, Messergebnis),
  Einheit = ifelse(Einheit == "ng/l", "µg/l", Einheit))
#2022-2025
PSM_metabolites_B0 <- rbind(PSM_metabolites_2022_B0,PSM_metabolites_2023_B0,PSM_metabolites_2024_B0,PSM_metabolites_2025_B0)
View(PSM_metabolites_B0)
PSM_metabolites_BB <- rbind(PSM_metabolites_2022_BB,PSM_metabolites_2023_BB,PSM_metabolites_2024_BB,PSM_metabolites_2025_BB)
View(PSM_metabolites_BB)

#####################################Tables for every year##############################################################################################################
#Drugs and Metabolites
drugs_metabolites_2022_B0 <-  drugs_metabolites_B0 |> filter(grepl("2022$", Datum))
drugs_metabolites_2023_B0 <-  drugs_metabolites_B0 |> filter(grepl("2023$", Datum))
drugs_metabolites_2024_B0 <-  drugs_metabolites_B0 |> filter(grepl("2024$", Datum))
drugs_metabolites_2025_B0 <-  drugs_metabolites_B0 |> filter(grepl("2025$", Datum))

drugs_metabolites_2022_BB <-  drugs_metabolites_BB |> filter(grepl("2022$", Datum))
drugs_metabolites_2023_BB <-  drugs_metabolites_BB |> filter(grepl("2023$", Datum))
drugs_metabolites_2024_BB <-  drugs_metabolites_BB |> filter(grepl("2024$", Datum))
drugs_metabolites_2025_BB <-  drugs_metabolites_BB |> filter(grepl("2025$", Datum))
#Bei PSM schon durch Download vorhanden

########################################Descriptive Data########################################################################################################################################
###Drugs and Metabolites###
number_locations_drugs_metabolites <-length(unique(drugs_metabolites$Messstellenname)) 
number_locations_drugs_metabolites # 652 locations
number_substances_drugs_metabolites <- length(unique(drugs_metabolites$Stoffname))
number_substances_drugs_metabolites # 85 substances

#Descriptive data Tables per year
table_22_drugs <- table(drugs_metabolites2022_B0$Messstellennummer)
length(unique(names(table_22_drugs))) # 262 locations 2022
table_23_drugs <- table(drugs_metabolites2023_B0$Messstellennummer)
length(unique(names(table_23_drugs))) # 206 locations 2023
table_24_drugs<- table(drugs_metabolites2024_B0$Messstellennummer)
length(unique(names(table_24_drugs))) # 293 locations 2024
table_25_drugs <- table(drugs_metabolites2025_B0$Messstellennummer)
length(unique(names(table_25_drugs))) # 105 locations 2025

normalize_names <- function(x) {
  trimws(tolower(x))  # entfernt Leerzeichen, setzt auf Kleinschreibung
}
names_22_drugs <- normalize_names(names(table_22_drugs))
names_23_drugs <- normalize_names(names(table_23_drugs))
names_24_drugs <- normalize_names(names(table_24_drugs))
names_25_drugs <- normalize_names(names(table_25_drugs))

setdiff(names(table_22_drugs), names(table_23_drugs)) ## 207 Messstellen, die nicht in Tabelle 23 vorkommen --> 55 gleiche Messstellen
setdiff(names(table_23_drugs), names(table_22_drugs)) ## 151 Messstellen, die nicht in Tabelle 22 vorkommen --> 55 gleiche Messstellen 

setdiff(names(table_24_drugs), names(table_22_drugs)) ## 237 Messstellen, die nicht in Tabelle 22 vorkommen --> 56 gleiche Messstellen
setdiff(names(table_22_drugs), names(table_24_drugs)) ## 206 Messstellen, die nicht in Tabelle 24 vorkommen --> 56 gleiche Messstellen

setdiff(names(table_22_drugs), names(table_25_drugs)) ## 237 Messstellen, die nicht in Tabelle 25 vorkommen --> 25 gleiche Messstellen
setdiff(names(table_25_drugs), names(table_22_drugs)) ## 80 Messstellen, die nicht in Tabelle 22 vorkommen --> 25 gleiche Messstellen 

setdiff(names(table_23_drugs), names(table_24_drugs)) ## 174 Messstellen, die nicht in Tabelle 24 vorkommen --> 32 gleiche Messstellen
setdiff(names(table_24_drugs), names(table_23_drugs)) ## 261 Messstellen, die nicht in Tabelle 23 vorkommen --> 32 gleiche Messstellen

setdiff(names(table_23_drugs), names(table_25_drugs)) ## 192 Messstellen, die nicht in Tabelle 25 vorkommen --> 14 gleiche Messstellen
setdiff(names(table_25_drugs), names(table_23_drugs)) ## 91 Messstellen, die nicht in Tabelle 23 vorkommen --> 14 gleiche Messstellen 

setdiff(names(table_24_drugs), names(table_25_drugs)) ## 267 Messstellen, die nicht in Tabelle 25 vorkommen --> 26 gleiche Messstellen
setdiff(names(table_25_drugs), names(table_24_drugs)) ## 79 Messstellen, die nicht in Tabelle 24 vorkommen --> 26 gleiche Messstellen 

number_substances_drugs_metabolites_2022 <- length(unique(drugs_metabolites_2022_B0$Stoffname))
number_substances_drugs_metabolites_2022 #79 substances 2025
number_substances_drugs_metabolites_2023 <- length(unique(drugs_metabolites_2023_B0$Stoffname))
number_substances_drugs_metabolites_2023 #82 substances 2025
number_substances_drugs_metabolites_2024 <- length(unique(drugs_metabolites_2024_B0$Stoffname))
number_substances_drugs_metabolites_2024 #79 substances 2025
number_substances_drugs_metabolites_2025 <- length(unique(drugs_metabolites_2025_B0$Stoffname))
number_substances_drugs_metabolites_2025 #79 substaces 2025

###PSM and Metabolites###
number_locations_PSM_metabolites <-length(unique(PSM_metabolites_B0$Messstellenname)) 
number_locations_PSM_metabolites # 736 locations
number_substances_PSM_metabolites <- length(unique(PSM_metabolites_B0$Stoffname))
number_substances_PSM_metabolites # 217 substances

#Descriptive Data per Year
number_locations_PSM_metabolites_2022 <-length(unique(PSM_metabolites_2022_B0$Messstellenname)) 
number_locations_PSM_metabolites_2022 # 279 locations 2022
number_substances_PSM_metabolites_2022 <- length(unique(PSM_metabolites_2022_B0$Stoffname))
number_substances_PSM_metabolites_2022 # 171 substances 2022

number_locations_PSM_metabolites_2023 <-length(unique(PSM_metabolites_2023_B0$Messstellenname)) 
number_locations_PSM_metabolites_2023 # 267 locations 2023
number_substances_PSM_metabolites_2023 <- length(unique(PSM_metabolites_2023_B0$Stoffname))
number_substances_PSM_metabolites_2023 # 210 substances 2023

number_locations_PSM_metabolites_2024 <-length(unique(PSM_metabolites_2024_B0$Messstellenname)) 
number_locations_PSM_metabolites_2024 # 291 locations 2024
number_substances_PSM_metabolites_2024 <- length(unique(PSM_metabolites_2024_B0$Stoffname))
number_substances_PSM_metabolites_2024 # 205 substances 2024

number_locations_PSM_metabolites_2025 <-length(unique(PSM_metabolites_2025_B0$Messstellenname)) 
number_locations_PSM_metabolites_2025 # 168 locations 2025
number_substances_PSM_metabolites_2025 <- length(unique(PSM_metabolites_2025_B0$Stoffname))
number_substances_PSM_metabolites_2025 # 205 substance 2025

###########################Barplot#################################################################################################################################
#Drugs and Metabolites
drugs_metabolites_top10 <- drugs_metabolites_B0 |>
  group_by(Stoffname) |>
  slice_max(Messergebnis, n = 1, with_ties = FALSE) |>
  ungroup() |>
  slice_max(Messergebnis, n = 10)

drugs_metabolites_top10 <- mutate(drugs_metabolites_top10, 
                                  Messstellenname = paste(Messstellenname, Datum, sep = " "))

plot_drugs_metabolites_top10 <- ggplot(drugs_metabolites_top10, 
                                       aes(x =  reorder(Stoffname, -Messergebnis), 
                                           y = Messergebnis, 
                                           fill = Messstellenname))+
  geom_bar(stat="identity")+
  labs(title = "Top 10 Substances",
       x = "Substancename",
       y = "Substance concentration [µg/l]",
       fill = "Measuring point")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_drugs_metabolites_top10
#Barplot per year
#2022
drugs_metabolites_top10_2022 <- drugs_metabolites_2022_B0 |>
  group_by(Stoffname) |>
  slice_max(Messergebnis, n = 1, with_ties = FALSE) |>
  ungroup() |>
  slice_max(Messergebnis, n = 10)

drugs_metabolites_top10_2022 <- mutate(drugs_metabolites_top10_2022, 
                                  Messstellenname = paste(Messstellenname, Datum, sep = " "))

plot_drugs_metabolites_top10_2022 <- ggplot(drugs_metabolites_top10_2022, 
                                       aes(x =  reorder(Stoffname, -Messergebnis), 
                                           y = Messergebnis, 
                                           fill = Messstellenname))+
  geom_bar(stat="identity")+
  labs(title = "Top 10 Substances 2022",
       x = "Substancename",
       y = "Substance concentration [µg/l]",
       fill = "Measuring point")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot_drugs_metabolites_top10_2022

#2023
drugs_metabolites_top10_2023 <- drugs_metabolites_2023_B0 |>
  group_by(Stoffname) |>
  slice_max(Messergebnis, n = 1, with_ties = FALSE) |>
  ungroup() |>
  slice_max(Messergebnis, n = 10)

drugs_metabolites_top10_2023 <- mutate(drugs_metabolites_top10_2023, 
                                       Messstellenname = paste(Messstellenname, Datum, sep = " "))

plot_drugs_metabolites_top10_2023 <- ggplot(drugs_metabolites_top10_2023, 
                                            aes(x =  reorder(Stoffname, -Messergebnis), 
                                                y = Messergebnis, 
                                                fill = Messstellenname))+
  geom_bar(stat="identity")+
  labs(title = "Top 10 Substances 2023",
       x = "Substancename",
       y = "Substance concentration [µg/l]",
       fill = "Measuring point")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot_drugs_metabolites_top10_2023

#2024
drugs_metabolites_top10_2024 <- drugs_metabolites_2024_B0 |>
  group_by(Stoffname) |>
  slice_max(Messergebnis, n = 1, with_ties = FALSE) |>
  ungroup() |>
  slice_max(Messergebnis, n = 10)

drugs_metabolites_top10_2024 <- mutate(drugs_metabolites_top10_2024, 
                                       Messstellenname = paste(Messstellenname, Datum, sep = " "))

plot_drugs_metabolites_top10_2024 <- ggplot(drugs_metabolites_top10_2024, 
                                            aes(x =  reorder(Stoffname, -Messergebnis), 
                                                y = Messergebnis, 
                                                fill = Messstellenname))+
  geom_bar(stat="identity")+
  labs(title = "Top 10 Substances 2024",
       x = "Substancename",
       y = "Substance concentration [µg/l]",
       fill = "Measuring point")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot_drugs_metabolites_top10_2024

#2025
drugs_metabolites_top10_2025 <- drugs_metabolites_2025_B0 |>
  group_by(Stoffname) |>
  slice_max(Messergebnis, n = 1, with_ties = FALSE) |>
  ungroup() |>
  slice_max(Messergebnis, n = 10)

drugs_metabolites_top10_2025 <- mutate(drugs_metabolites_top10_2025, 
                                       Messstellenname = paste(Messstellenname, Datum, sep = " "))

plot_drugs_metabolites_top10_2025 <- ggplot(drugs_metabolites_top10_2025, 
                                            aes(x =  reorder(Stoffname, -Messergebnis), 
                                                y = Messergebnis, 
                                                fill = Messstellenname))+
  geom_bar(stat="identity")+
  labs(title = "Top 10 Substances 2025",
       x = "Substancename",
       y = "Substance concentration [µg/l]",
       fill = "Measuring point")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot_drugs_metabolites_top10_2025


###PSM and Metabolites###
PSM_metabolites_top10 <- PSM_metabolites_B0 |>
  group_by(Stoffname) |>
  slice_max(Messergebnis, n = 1, with_ties = FALSE) |>
  ungroup() |>
  slice_max(Messergebnis, n = 10)

PSM_metabolites_top10 <- mutate(PSM_metabolites_top10, 
                                  Messstellenname = paste(Messstellenname, Datum, sep = " "))

plot_PSM_metabolites_top10 <- ggplot(PSM_metabolites_top10, 
                                       aes(x =  reorder(Stoffname, -Messergebnis), 
                                           y = Messergebnis, 
                                           fill = Messstellenname))+
  geom_bar(stat="identity")+
  labs(title = "Top 10 Substances",
       x = "Substancename",
       y = "Substance concentration [µg/l]",
       fill = "Measuring point")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_PSM_metabolites_top10

#Barplot per year
#2022
PSM_metabolites_top10_2022 <- PSM_metabolites_2022_B0 |>
  group_by(Stoffname) |>
  slice_max(Messergebnis, n = 1, with_ties = FALSE) |>
  ungroup() |>
  slice_max(Messergebnis, n = 10)

PSM_metabolites_top10_2022 <- mutate(PSM_metabolites_top10_2022, 
                                       Messstellenname = paste(Messstellenname, Datum, sep = " "))

plot_PSM_metabolites_top10_2022 <- ggplot(PSM_metabolites_top10_2022, 
                                            aes(x =  reorder(Stoffname, -Messergebnis), 
                                                y = Messergebnis, 
                                                fill = Messstellenname))+
  geom_bar(stat="identity")+
  labs(title = "Top 10 Substances 2022",
       x = "Substancename",
       y = "Substance concentration [µg/l]",
       fill = "Measuring point")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot_PSM_metabolites_top10_2022

#2023
PSM_metabolites_top10_2023 <- PSM_metabolites_2023_B0 |>
  group_by(Stoffname) |>
  slice_max(Messergebnis, n = 1, with_ties = FALSE) |>
  ungroup() |>
  slice_max(Messergebnis, n = 10)

PSM_metabolites_top10_2023 <- mutate(PSM_metabolites_top10_2023, 
                                     Messstellenname = paste(Messstellenname, Datum, sep = " "))

plot_PSM_metabolites_top10_2023 <- ggplot(PSM_metabolites_top10_2023, 
                                          aes(x =  reorder(Stoffname, -Messergebnis), 
                                              y = Messergebnis, 
                                              fill = Messstellenname))+
  geom_bar(stat="identity")+
  labs(title = "Top 10 Substances 2023",
       x = "Substancename",
       y = "Substance concentration [µg/l]",
       fill = "Measuring point")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot_PSM_metabolites_top10_2023

#2024
PSM_metabolites_top10_2024 <- PSM_metabolites_2024_B0 |>
  group_by(Stoffname) |>
  slice_max(Messergebnis, n = 1, with_ties = FALSE) |>
  ungroup() |>
  slice_max(Messergebnis, n = 10)

PSM_metabolites_top10_2024 <- mutate(PSM_metabolites_top10_2024, 
                                     Messstellenname = paste(Messstellenname, Datum, sep = " "))

plot_PSM_metabolites_top10_2024 <- ggplot(PSM_metabolites_top10_2024, 
                                          aes(x =  reorder(Stoffname, -Messergebnis), 
                                              y = Messergebnis, 
                                              fill = Messstellenname))+
  geom_bar(stat="identity")+
  labs(title = "Top 10 Substances 2024",
       x = "Substancename",
       y = "Substance concentration [µg/l]",
       fill = "Measuring point")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot_PSM_metabolites_top10_2024

#2025
PSM_metabolites_top10_2025 <- PSM_metabolites_2025_B0 |>
  group_by(Stoffname) |>
  slice_max(Messergebnis, n = 1, with_ties = FALSE) |>
  ungroup() |>
  slice_max(Messergebnis, n = 10)

PSM_metabolites_top10_2025 <- mutate(PSM_metabolites_top10_2025, 
                                     Messstellenname = paste(Messstellenname, Datum, sep = " "))

plot_PSM_metabolites_top10_2025 <- ggplot(PSM_metabolites_top10_2025, 
                                          aes(x =  reorder(Stoffname, -Messergebnis), 
                                              y = Messergebnis, 
                                              fill = Messstellenname))+
  geom_bar(stat="identity")+
  labs(title = "Top 10 Substances 2025",
       x = "Substancename",
       y = "Substance concentration [µg/l]",
       fill = "Measuring point")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot_PSM_metabolites_top10_2025

##################################Calculate Risk#######################################################################################################################
###Drugs and Metabolites###
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
