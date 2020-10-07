##### TITLE #####
# EEMB 144L Intro to R with CAL FIRE Major Fires 2013 - 2019 Dataset
# Nicholas Baetge
# 10/2/20202

##### Projects #####

### Create a new project with File > New Project. 

### When you open the project, R will recognize that as your working directory. Notice that the pathway shows up in the top bar of RStudio, and any files that you save/add are automatically added to that working directory. A folder for the entire project is also created.

# We're working in a SCRIPT. It's like a text editor (as opposed to the always-active console) and only runs lines of code when you ask it to. It also allows you to keep a clear record of the work you have done. 

### Load libraries ###

library(tidyverse)
library(readxl)

##### Loading data from .csv or xlsx files #####

### Loading files into RStudio is EASY if you're working in a project. All you have to do is drag and drop the file into your project folder.

# Drag and drop into the project folder, and note that it appears in the 'Files' tab in RStudio automatically. That's because R knows that it's been added to our working directory. 

# Once it's in the Project, then it's easy to load using read_csv() or read_excel().

# unlike .csv files, .xlsx files can have multiple sheets
excel_sheets("Input_Data/week1/2013_2019_CalFire_Redbook.xlsx") #let's see what the excel sheets are called

calfire.data <- read_excel("~/GITHUB/144l_students/Input_Data/week1/2013_2019_CALFIRE_Redbook.xlsx", sheet = "Data") # we can store each sheet separately as a dataframe! this one is the data
view(calfire.data)

# shorcut for "<-" is "alt -" on a pc and "alt/option -" on a mac 

calfire.metadata <- read_excel("Input_Data/week1/2013_2019_CalFire_Redbook.xlsx", sheet = "Metadata")  # this one is the metadata

##### Initial data exploration ####

names(calfire.data) #shows variable (column) names
dim(calfire.data) #dimensions of dataset (rows, columns)
class(calfire.data) #data class
head(calfire.data) #shows first 6 lines of dataset
tail(calfire.data) #shows last 6 lines of dataset

# Want to know how a function works?

?names # single ? brings up R doc for that function
??names # double ?? brings up every function that might contain "names"

# Single columns can be referred to using a '$'
county <- calfire.data$County_Unit

max_acres <- max(calfire.data$Total_Acres_Burned, na.rm = T)
max(calfire.data$Structures_Destroyed)
max(calfire.data$Structures_Destroyed, na.rm = T)

##### Basic data wrangling (dplyr functions) ####

df1 <- select(calfire.data, County_Unit:Controlled_Date, Total_Acres_Burned:Civil_Fatalities )
view(df1)

df2 <- filter(df1, County_Unit %in% c("SANTA BARBARA", "VENTURA", "LOS ANGELES", "ORANGE", "SAN DIEGO") & Total_Acres_Burned >= 500 | Fire_Name == "THOMAS")

df3 <- arrange(df2, desc(Start_Date), Total_Acres_Burned)

df4 <- mutate_at(df3, vars(Structures_Destroyed:Civil_Fatalities), replace_na, 0)

df5 <- mutate(df4, Fatalities = Fire_Fatalities + Civil_Fatalities)

#mess with time
install.packages("lubridate")
library(lubridate)

df6 <- mutate(df5,
              interv = interval(Start_Date, Controlled_Date),
              dur = as.duration(interv),
              days = as.numeric(dur, "days"))

## We used 15 lines to do all of that! Now we have 5 dataframes! This seems a little inefficient. There is a better way - it's called "piping"

### Intro to piping ####

# We want to restrict our data to the SoCal coast, exclude fires that burned less than 500 acres, add column that sums the number of fatalities, change NAs to 0s, arrange data

# The MAGICAL pipe operator: %>% (command + shift + m on a mac, control + shift + m on a windows)

# Think of this as a code way of saying "and then..."

socal_fires <- calfire.data %>% 
  filter(County_Unit %in% c("SANTA BARBARA", "VENTURA", "LOS ANGELES", "ORANGE", "SAN DIEGO") & Total_Acres_Burned >= 500 | Fire_Name == "THOMAS") %>% 
  mutate_at(vars(Structures_Destroyed:Civil_Fatalities), replace_na, 0) %>% 
  mutate(Fatalities = Fire_Fatalities + Civil_Fatalities,
         interv = interval(Start_Date, Controlled_Date),
         dur = as.duration(interv),
         days = as.numeric(dur, "days"),
         County_Unit = ifelse(County_Unit == "VENTURA/SANTA BARBARA", "VENTURA", County_Unit)
         ) %>% 
  arrange(desc(Start_Date), Total_Acres_Burned)

view(socal_fires)

#### Our first graphs in ggplot ####

# Three things you must tell R to make a graph in ggplot
# (1) That you're using ggplot
# (2) What data you're using (including what should be x and what should be y)
# (3) Type of graph that you want to create

socal.plot <- socal_fires %>%
  rename(start = Start_Date,
         acres = Total_Acres_Burned) %>% 
  ggplot(aes(x = start, y = acres)) +
  geom_point()

socal.plot <- socal_fires %>%
  rename(start = Start_Date,
         acres = Total_Acres_Burned,
         county = County_Unit) %>% 
  ggplot(aes(x = start, y = acres)) +
  geom_point(aes(color = county)) +
  ggtitle("California SoCal Major Fires 2013 - 2018") +
  xlab("Date") +
  ylab("Acres Burned") +
  theme(panel.grid.major = element_blank())   

socal.plot + facet_grid(~county)

  
  
