#### Intro to R with the CAL FIRE MAJOR FIRES dataset (2013 - 2019) ####
# Nicholas Baetge
# 10/8/2020

#### Load Packages ####

# install.packages("tidyverse")
library(tidyverse)
?tidyverse
# install.packages("dplyr")
library(dplyr)
?dplyr
# install.packages("readxl")
library(readxl)
# install.packages("praise")
library(praise)

# to clear your environments you can type in the following:
# rm(list = ls())

#### Load Dataset ####

excel_sheets("Input_Data/week1/2013_2019_CALFIRE_Redbook.xlsx")

metadata <- read_excel("Input_Data/week1/2013_2019_CALFIRE_Redbook.xlsx", sheet = "Metadata") 
view(metadata)
data <- read_excel("Input_Data/week1/2013_2019_CALFIRE_Redbook.xlsx", sheet = "Data")
view(data)

#### Initial Data Exploration ####

names(data)
?names
dim(data)
class(data)
head(data)
tail(data)
str(data)
glimpse(data)
typeof(data$Total_Acres_Burned) # single columns can be referred to using '$'

max(data$Total_Acres_Burned)
max(data$Structures_Destroyed)
?max
max(data$Structures_Destroyed, na.rm = T)

summary(data)

#### Basic data wrangling (dplyr functions) ####

df1 <- select(data, County_Unit:Controlled_Date, Total_Acres_Burned, Cause:Structures_Damaged)
unique(df1$County_Unit)

df2 <- filter(df1, County_Unit %in% c("SANTA BARBARA", "VENTURA", "LOS ANGELES", "SAN DIEGO", "ORANGE", "VENTURA/SANTA BARBARA") & Total_Acres_Burned >= 500)

# | = "or" 
# == = "equals"/"matches" , %in% c() 
# & = "and"

df3 <- arrange(df2, desc(Total_Acres_Burned))

df4 <- mutate_at(df3, vars("Structures_Destroyed", "Structures_Damaged"), replace_na, 0)

df5 <- mutate(df4, struc_impact = Structures_Damaged + Structures_Destroyed)

#mess with time
library(lubridate)

df6 <- mutate(df5, interv = interval(Start_Date, Controlled_Date), 
              dur = as.duration(interv),
              days = as.numeric(dur, "days"))


#### Introduction to Piping ####

socal.fires <- data %>% # mac : cmd + shift + m, pc, ctrl + shift + m
  select(County_Unit:Controlled_Date, Total_Acres_Burned, Cause:Structures_Damaged) %>% 
  filter(County_Unit %in% c("SANTA BARBARA", "VENTURA", "LOS ANGELES", "SAN DIEGO", "ORANGE", "VENTURA/SANTA BARBARA") & Total_Acres_Burned >= 500) %>% 
  arrange(desc(Total_Acres_Burned)) %>% 
  mutate_at(vars("Structures_Destroyed", "Structures_Damaged"), replace_na, 0) %>% 
  mutate(struc_impact = Structures_Damaged + Structures_Destroyed,
         interv = interval(Start_Date, Controlled_Date), 
         dur = as.duration(interv),
         days = as.numeric(dur, "days")) 

##### Our first graphs in ggplot #####

# We're going to make a graph of acres burned in the South Coast from 2013 - 2018, with the color dependent on which county we're showing. 

# Three things you must tell R to make a graph in ggplot: 
# (1) That you're using ggplot
# (2) What data you're using (including what should be x and what should be y)
# (3) What type of graph you want to create
# Everything after that is extra to make it beautiful

ggplot(socal.fires, aes(x = Start_Date, y = Total_Acres_Burned)) +
  geom_point(aes(color = County_Unit)) +
  ggtitle("CA South Coast Major Fires \n2014 - 2018") + 
  labs(x = "", y = "Total Acres Burned", color = "County") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  facet_grid(rows = "County_Unit", scales = "free")

plot.data <- socal.fires %>% 
  rename(county = County_Unit,
         acres = Total_Acres_Burned,
         start = Start_Date,
         end = Controlled_Date) %>% 
  mutate(year = year(start),
         county = ifelse(county == "VENTURA/SANTA BARBARA", "VENTURA", county)) 

incidents <- plot.data %>% 
  group_by(county, year) %>% 
  tally()
  ungroup()
  
incidents.plot <- incidents %>% 
  ggplot(aes(x = year, y = n)) +
  # geom_point(aes(color = county)) +
  # geom_line(aes(color = county)) +
  geom_point(color = "blue") +
  geom_line(color = "blue") +
  labs(title = "CA South Coast Major Fire Incidents \n 2014 -2018", x = "", y = "Incidents", color = "County") + 
  theme_bw() +
  facet_grid(rows = "county", scales = "free") 
  # guides(color = F)


all_incidents <- plot.data %>% 
  group_by(year) %>% 
  tally() %>% 
  ungroup()

all_incidents.plot <- all_incidents %>% 
  ggplot(aes(x = year, y = n)) +
  geom_point(color = "blue") +
  geom_line(color = "blue") +
  labs(title = "CA South Coast Major Fire Incidents \n 2014 -2018", x = "", y = "Incidents") + 
  theme_bw() 

##### Save Data and Plots ####

saveRDS(socal.fires, file = "Output_Data/week1/socal_fires_data.rds")
write_csv(socal.fires, "Output_Data/week1/socal_fires_data.csv")

ggsave(filename = "Fire_Incidents", all_incidents.plot, device = "jpeg", "Output_Data/week1/" )











