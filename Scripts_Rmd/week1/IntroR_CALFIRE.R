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

calfire.data <- read_excel("Input_Data/week1/2013_2019_CalFire_Redbook.xlsx", sheet = "Data") # we can store each sheet separately as a dataframe! this one is the data

# shorcut for "<-" is "alt -" on a pc and "alt/option -" on a mac 

calfire.metadata <- read_excel("Input_Data/week1/2013_2019_CalFire_Redbook.xlsx", sheet = "Metadata")  # this one is the metadata





