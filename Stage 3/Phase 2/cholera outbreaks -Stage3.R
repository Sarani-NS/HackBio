library(ggplot2)
library(plotly)
library(dplyr)
library(tidyr)
library(maps)
library(htmlwidgets)

#Data loading & checking
Cases_of_cholera <- read.csv('CHOLERA_0000000001.csv')
glimpse(Cases_of_cholera)
Cases_unique_values <- lapply(Cases_of_cholera, unique)

Deaths_from_cholera <- read.csv('CHOLERA_0000000002.csv')
glimpse(Deaths_from_cholera)
Deaths_unique_values <- lapply(Deaths_from_cholera, unique)

Fatality_rate <- read.csv('CHOLERA_0000000003.csv')
glimpse(Fatality_rate)
Fatality_rate_unique_values <- lapply(Fatality_rate, unique)

#Checking missing values
sum(is.na(Cases_of_cholera))
sum(is.na(Deaths_from_cholera))
sum(is.na(Fatality_rate))

#Editing data types
Deaths_from_cholera$Number.of.reported.deaths.from.cholera <- as.integer(Deaths_from_cholera$Number.of.reported.deaths.from.cholera)
glimpse(Deaths_from_cholera)

Fatality_rate$Cholera.case.fatality.rate <- as.double(Fatality_rate$Cholera.case.fatality.rate)
glimpse(Fatality_rate)

#Visualization
#Report and visualize cholera outbreaks since 1949
#Interactive bar plot for Countries with highest cases
highest_cases <- Cases_of_cholera %>%
  group_by(Countries..territories.and.areas) %>%
  summarise(total_cases = sum(Number.of.reported.cases.of.cholera, na.rm = TRUE)) %>%
  arrange(desc(total_cases))

bar_plot <- plot_ly(highest_cases, 
                    x = ~Countries..territories.and.areas, 
                    y = ~total_cases, 
                    type = 'bar',
                    marker = list(color = 'darkred')) %>%
  layout(title = "Total Cholera Cases by Country",
         xaxis = list(title = "Country"),
         yaxis = list(title = "Total Cases"))

bar_plot  
saveWidget(bar_plot, "total_cholera_cases_by_country.html")

#Top 5 Countries with the Highest Total Number of Reported Cases
highest_cases_top5 <- Cases_of_cholera %>%
  group_by(Countries..territories.and.areas) %>%
  summarise(total_cases = sum(Number.of.reported.cases.of.cholera, na.rm = TRUE)) %>%
  arrange(desc(total_cases)) %>%
  top_n(5, total_cases)  

bar_plot_top5 <- plot_ly(highest_cases_top5, 
                         x = ~Countries..territories.and.areas, 
                         y = ~total_cases, 
                         type = 'bar',
                         marker = list(color = 'darkred')) %>%
  layout(title = "Top 5 Countries with the Highest Cholera Cases",
         xaxis = list(title = "Country"),
         yaxis = list(title = "Total Cases"))

bar_plot_top5  
saveWidget(bar_plot_top5, "top_5_cholera_cases_by_country.html")

#LOWEST
#Top 10 Countries with the lowest Total Number of Reported Cases
lowest_10_cases <- Cases_of_cholera %>%
  group_by(Countries..territories.and.areas) %>%
  summarise(total_cases = sum(Number.of.reported.cases.of.cholera, na.rm = TRUE)) %>%
  arrange(total_cases) %>%
  slice_head(n = 10)  

bar_plot_top3 <- plot_ly(lowest_10_cases, 
                         x = ~Countries..territories.and.areas, 
                         y = ~total_cases, 
                         type = 'bar',
                         marker = list(color = 'darkred')) %>%
  layout(title = "Top 10 Countries with the Lowest Cholera Cases",
         xaxis = list(title = "Country"),
         yaxis = list(title = "Total Cases"))

bar_plot_top3  
saveWidget(bar_plot_top3, "Top_10_Countries_with_Lowest_Cholera_Cases.html")
# Find the country and year with the highest number of cholera Cases
highest_Cases_record <- Cases_of_cholera %>%
  filter(Number.of.reported.cases.of.cholera == max(Number.of.reported.cases.of.cholera, na.rm = TRUE))
# Print the result (country and year)
print(highest_Cases_record %>% select(Countries..territories.and.areas, Year, Number.of.reported.cases.of.cholera))

# List of all countries from the map data
all_countries_Cases <- data.frame(Countries = unique(map_data("world")$region))

highest_cases_top5_full <- all_countries_Cases %>%
  left_join(highest_cases_top5, by = c("Countries" = "Countries..territories.and.areas")) %>%
  replace_na(list(total_cases = 0)) 

# Interactive World Map 
map_plot_1 <- plot_ly(highest_cases_top5_full, 
                    type = 'choropleth', 
                    locations = ~Countries, 
                    locationmode = 'country names', 
                    z = ~total_cases, 
                    colorscale = list(c(0, 'rgb(240,240,240)'),  
                                      c(1, 'rgb(255,0,0)')),     
                    marker = list(line = list(color = 'rgb(0,0,0)', width = 0.5))) %>%  
  layout(title = "Top 5 Countries with the Highest Cholera Cases",
         geo = list(showframe = FALSE,
                    showcoastlines = TRUE,  
                    projection = list(type = 'natural earth')))

map_plot_1 
saveWidget(map_plot_1, "Map_top_5_cholera_cases_by_country.html")
# 2. First Country to Have the Highest Cases in the Oldest Year
first_highest_case <- Cases_of_cholera %>%
  arrange(Year) %>%
  filter(Year == min(Year)) %>%
  arrange(desc(Number.of.reported.cases.of.cholera)) %>%
  top_n(1, Number.of.reported.cases.of.cholera)

print(first_highest_case)
#line plot for total cases over years
line_plot_1 <- Cases_of_cholera %>%
  group_by(Year) %>%
  summarise(total_cases = sum(Number.of.reported.cases.of.cholera, na.rm = TRUE)) %>%
  plot_ly(x = ~Year, 
          y = ~total_cases, 
          type = 'scatter', 
          mode = 'lines+markers', 
          line = list(color = 'darkred')) %>%
  layout(title = "Cholera Cases Over Time",
         xaxis = list(title = "Year"),
         yaxis = list(title = "Total Cases"))

line_plot_1 
saveWidget(line_plot_1, "Line_plot_cholera_cases_by_country.html")
#Deaths_from_cholera visualization 

# Visualization
# Report and visualize cholera deaths since 1949

# Interactive bar plot for Countries with highest deaths
highest_deaths <- Deaths_from_cholera %>%
  group_by(Countries..territories.and.areas) %>%
  summarise(total_deaths = sum(Number.of.reported.deaths.from.cholera, na.rm = TRUE)) %>%
  arrange(desc(total_deaths))

bar_plot_d1 <- plot_ly(highest_deaths, 
                    x = ~Countries..territories.and.areas, 
                    y = ~total_deaths, 
                    type = 'bar',
                    marker = list(color = 'darkgray')) %>%
  layout(title = "Total Cholera Deaths by Country",
         xaxis = list(title = "Country"),
         yaxis = list(title = "Total Deaths"))

bar_plot_d1  
saveWidget(bar_plot_d1, "total_cholera_deaths_by_country.html")
#Top 5 Countries with the Highest Total Number of Reported Deaths
highest_deaths_top5 <- Deaths_from_cholera %>%
  group_by(Countries..territories.and.areas) %>%
  summarise(total_deaths = sum(Number.of.reported.deaths.from.cholera, na.rm = TRUE)) %>%
  arrange(desc(total_deaths)) %>%
  top_n(5, total_deaths)  

bar_plot_top5_d <- plot_ly(highest_deaths_top5, 
                         x = ~Countries..territories.and.areas, 
                         y = ~total_deaths, 
                         type = 'bar',
                         marker = list(color = 'darkgray')) %>%
  layout(title = "Top 5 Countries with the Highest Cholera Deaths",
         xaxis = list(title = "Country"),
         yaxis = list(title = "Total Deaths"))

bar_plot_top5_d  
saveWidget(bar_plot_top5_d, "top_5_countries_cholera_deaths.html")

# Find the country and year with the highest number of cholera deaths
highest_death_record <- Deaths_from_cholera %>%
  filter(Number.of.reported.deaths.from.cholera == max(Number.of.reported.deaths.from.cholera, na.rm = TRUE))
# Print the result (country and year)
print(highest_death_record %>% select(Countries..territories.and.areas, Year, Number.of.reported.deaths.from.cholera))

# List of all countries from the map data
all_countries_deaths <- data.frame(Countries = unique(map_data("world")$region))

highest_deaths_top5_full <- all_countries_deaths %>%
  left_join(highest_deaths_top5, by = c("Countries" = "Countries..territories.and.areas")) %>%
  replace_na(list(total_deaths = 0)) 

# Interactive World Map 
map_plot_2 <- plot_ly(highest_deaths_top5_full, 
                    type = 'choropleth', 
                    locations = ~Countries, 
                    locationmode = 'country names', 
                    z = ~total_deaths, 
                    colorscale = list(c(0, 'rgb(220,220,220)'),  # Lighter gray for lower values
                                      c(0.25, 'rgb(180,180,180)'),  
                                      c(0.5, 'rgb(140,140,140)'),  
                                      c(0.75, 'rgb(100,100,100)'),  
                                      c(1, 'rgb(50,50,50)')),     # Darkest gray for highest values
                    marker = list(line = list(color = 'rgb(0,0,0)', width = 0.5))) %>%  
  layout(title = "Top 5 Countries with the Highest Cholera Deaths",
         geo = list(showframe = FALSE,
                    showcoastlines = TRUE,  
                    projection = list(type = 'natural earth')))

map_plot_2  
saveWidget(map_plot_2 , "Map_of_deaths.html")
# First Country to Have the Highest Deaths in the Oldest Year
first_highest_death <- Deaths_from_cholera %>%
  arrange(Year) %>%
  filter(Year == min(Year)) %>%
  arrange(desc(Number.of.reported.deaths.from.cholera)) %>%
  top_n(1, Number.of.reported.deaths.from.cholera)

print(first_highest_death)

# Line plot showing cholera deaths over time
line_plot_2 <- Deaths_from_cholera %>%
  group_by(Year) %>%
  summarise(total_deaths = sum(Number.of.reported.deaths.from.cholera, na.rm = TRUE)) %>%
  plot_ly(x = ~Year, 
          y = ~total_deaths, 
          type = 'scatter', 
          mode = 'lines+markers', 
          line = list(color = 'darkgray')) %>%
  layout(title = "Cholera Deaths Over Time",
         xaxis = list(title = "Year"),
         yaxis = list(title = "Total Deaths"))

line_plot_2  
saveWidget(line_plot_2, "Line_plot_showing_cholera_deaths_over_time.html")
#Fatality rate visualization
#Interactive bar plot for Countries with highest Fatality rate
highest_rates <- Fatality_rate %>%
  group_by(Countries..territories.and.areas) %>%
  summarise(total_rates = sum(Cholera.case.fatality.rate, na.rm = TRUE)) %>%
  arrange(desc(total_rates))

bar_plot_f <- plot_ly(highest_rates, 
                    x = ~Countries..territories.and.areas, 
                    y = ~total_rates, 
                    type = 'bar',
                    marker = list(color = 'black')) %>%
  layout(title = "Total Cholera fatility rates by Country",
         xaxis = list(title = "Country"),
         yaxis = list(title = "Total rates"))

bar_plot_f  
saveWidget(bar_plot_f , "Countries with highest Fatality rate.html")
# Top 5 Countries with the Highest fatality rates
highest_rates_top5 <- Fatality_rate %>%
  group_by(Countries..territories.and.areas) %>%
  summarise(total_rates = sum(Cholera.case.fatality.rate, na.rm = TRUE)) %>%
  arrange(desc(total_rates)) %>%
  top_n(5, total_rates)  

bar_plot_top5 <- plot_ly(highest_rates_top5, 
                         x = ~Countries..territories.and.areas, 
                         y = ~total_rates, 
                         type = 'bar',
                         marker = list(color = 'black')) %>%
  layout(title = "Top 5 Countries with the Highest fatality rate",
         xaxis = list(title = "Country"),
         yaxis = list(title = "rates"))

bar_plot_top5  
saveWidget(bar_plot_top5, "Top 5 Countries with the Highest fatality rates.html")
#############################
#LOWEST
#40 Countries with the lowest Fatality rate
lowest_40_rates <- Fatality_rate %>%
  group_by(Countries..territories.and.areas) %>%
  summarise(total_rates = sum(Cholera.case.fatality.rate, na.rm = TRUE)) %>%
  arrange(total_rates) %>%
  slice_head(n = 40)  

bar_plot_40rates <- plot_ly(lowest_40_rates, 
                         x = ~Countries..territories.and.areas, 
                         y = ~total_rates, 
                         type = 'bar',
                         marker = list(color = 'black')) %>%
  layout(title = "40 Countries with the Lowest Fatality rate",
         xaxis = list(title = "Country"),
         yaxis = list(title = "rates"))

bar_plot_40rates 
saveWidget(bar_plot_40rates, "40 Countries with the Lowest Fatality rate.html")
# Find the country and year with the highest fatality rate
highest_rates_record <- Fatality_rate %>%
  filter(Cholera.case.fatality.rate == max(Cholera.case.fatality.rate, na.rm = TRUE))
# Print the result (country and year)
print(highest_rates_record %>% select(Countries..territories.and.areas, Year, Cholera.case.fatality.rate))

# List of all countries from the map data
all_countries_fatality <- data.frame(Countries = unique(map_data("world")$region))

highest_rates_top5_full <- all_countries_fatality %>%
  left_join(highest_rates_top5, by = c("Countries" = "Countries..territories.and.areas")) %>%
  replace_na(list(total_rates = 0)) 

# Interactive World Map 
map_plot_3 <- plot_ly(highest_rates_top5_full, 
                    type = 'choropleth', 
                    locations = ~Countries, 
                    locationmode = 'country names', 
                    z = ~total_rates, 
                    colorscale = list(c(0, 'rgb(240,240,240)'),  
                                      c(1, 'rgb(255,0,0)')),     
                    marker = list(line = list(color = 'rgb(0,0,0)', width = 0.5))) %>%  
  layout(title = "Top 5 Countries with the Highest fatality rates",
         geo = list(showframe = FALSE,
                    showcoastlines = TRUE,  
                    projection = list(type = 'natural earth')))

map_plot_3  
saveWidget(map_plot_3, "fatility rate map.html")

#First Country to Have the Highest fatality rate in the Oldest Year
first_highest_rate <- Fatality_rate %>%
  arrange(Year) %>%
  filter(Year == min(Year)) %>%
  arrange(desc(Cholera.case.fatality.rate)) %>%
  top_n(1, Cholera.case.fatality.rate)

print(first_highest_rate)
#interactive line plot
line_plot_3 <- Fatality_rate %>%
  group_by(Year) %>%
  summarise(total_rates = sum(Cholera.case.fatality.rate, na.rm = TRUE)) %>%
  plot_ly(x = ~Year, 
          y = ~total_rates, 
          type = 'scatter', 
          mode = 'lines+markers', 
          line = list(color = 'black')) %>%
  layout(title = "Cholera fatality rate Over Time",
         xaxis = list(title = "Year"),
         yaxis = list(title = "Total Cases"))

line_plot_3 
saveWidget(line_plot_3, "fatility rate line plot.html")
##################################################################3
#SPAIN Visualization
#SPAIN total reported cases
SPAIN_Reported_cases <- Cases_of_cholera %>%
  filter(Countries..territories.and.areas == 'Spain')
print(sum(SPAIN_Reported_cases$Number.of.reported.cases.of.cholera, na.rm = TRUE))
#SPAIN total reported deaths
SPAIN_Reported_deaths <- Deaths_from_cholera %>%
  filter(Countries..territories.and.areas == 'Spain')
print(sum(SPAIN_Reported_deaths$Number.of.reported.deaths.from.cholera, na.rm = TRUE))
#SPAIN fatality rate
SPAIN_Fatality_Rate <- Fatality_rate %>% 
  filter(Countries..territories.and.areas == 'Spain')
print(sum(SPAIN_Fatality_Rate$Cholera.case.fatality.rate, na.rm = TRUE))
#SPAIN reported cases over time interactive line plot
line_plot_Spain_cases <- Cases_of_cholera %>%
  filter(Countries..territories.and.areas == 'Spain') %>%
  group_by(Year) %>%
  summarise(total_cases = sum(Number.of.reported.cases.of.cholera, na.rm = TRUE)) %>%
  plot_ly(x = ~Year, 
          y = ~total_cases, 
          type = 'scatter', 
          mode = 'lines+markers', 
          line = list(color = 'orange')) %>%
  layout(title = "Spain reported cases over time",
         xaxis = list(title = "Year"),
         yaxis = list(title = "Total Cases"))

line_plot_Spain_cases 
saveWidget(line_plot_Spain_cases, "Spain cases.html")

#SPAIN reported deaths over time
line_plot_Spain_deaths <- Deaths_from_cholera %>%
  filter(Countries..territories.and.areas == 'Spain') %>%
  group_by(Year) %>%
  summarise(total_deaths = sum(Number.of.reported.deaths.from.cholera, na.rm = TRUE)) %>%
  plot_ly(x = ~Year, 
          y = ~total_deaths, 
          type = 'scatter', 
          mode = 'lines+markers', 
          line = list(color = 'orange')) %>%
  layout(title = "Spain reported deaths over time",
         xaxis = list(title = "Year"),
         yaxis = list(title = "Total Cases"))

line_plot_Spain_deaths 
saveWidget(line_plot_Spain_deaths, "Spain deaths.html")
#SPAIN fatality rate over time
line_plot_Spain_FR <- Fatality_rate %>%
  filter(Countries..territories.and.areas == 'Spain') %>%
  group_by(Year) %>%
  summarise(F_R = sum(Cholera.case.fatality.rate, na.rm = TRUE)) %>%
  plot_ly(x = ~Year, 
          y = ~F_R, 
          type = 'scatter', 
          mode = 'lines+markers', 
          line = list(color = 'orange')) %>%
  layout(title = "Spain fatality rate over time",
         xaxis = list(title = "Year"),
         yaxis = list(title = "Fatality Rate"))

line_plot_Spain_FR 
saveWidget(line_plot_Spain_FR, "Spain fatality rate.html")
###############################################################3
#India total reported cases
India_Reported_cases <- Cases_of_cholera %>%
  filter(Countries..territories.and.areas == 'India')
print(sum(India_Reported_cases$Number.of.reported.cases.of.cholera, na.rm = TRUE))
#India total reported deaths
India_Reported_deaths <- Deaths_from_cholera %>%
  filter(Countries..territories.and.areas == 'India')
print(sum(India_Reported_deaths$Number.of.reported.deaths.from.cholera, na.rm = TRUE))
#India fatality rate
India_Fatality_Rate <- Fatality_rate %>% 
  filter(Countries..territories.and.areas == 'India')
print(sum(India_Fatality_Rate$Cholera.case.fatality.rate, na.rm = TRUE))
#India reported cases over time interactive line plot
line_plot_India_cases <- Cases_of_cholera %>%
  filter(Countries..territories.and.areas == 'India') %>%
  group_by(Year) %>%
  summarise(total_cases = sum(Number.of.reported.cases.of.cholera, na.rm = TRUE)) %>%
  plot_ly(x = ~Year, 
          y = ~total_cases, 
          type = 'scatter', 
          mode = 'lines+markers', 
          line = list(color = 'orange')) %>%
  layout(title = "India reported cases over time",
         xaxis = list(title = "Year"),
         yaxis = list(title = "Total Cases"))

line_plot_India_cases 
saveWidget(line_plot_India_cases, "India cases.html")

#India reported deaths over time
line_plot_India_deaths <- Deaths_from_cholera %>%
  filter(Countries..territories.and.areas == 'India') %>%
  group_by(Year) %>%
  summarise(total_deaths = sum(Number.of.reported.deaths.from.cholera, na.rm = TRUE)) %>%
  plot_ly(x = ~Year, 
          y = ~total_deaths, 
          type = 'scatter', 
          mode = 'lines+markers', 
          line = list(color = 'orange')) %>%
  layout(title = "India reported deaths over time",
         xaxis = list(title = "Year"),
         yaxis = list(title = "Total Cases"))

line_plot_India_deaths 
saveWidget(line_plot_India_deaths, "India deaths.html")
#India fatality rate over time
line_plot_India_FR <- Fatality_rate %>%
  filter(Countries..territories.and.areas == 'India') %>%
  group_by(Year) %>%
  summarise(F_R = sum(Cholera.case.fatality.rate, na.rm = TRUE)) %>%
  plot_ly(x = ~Year, 
          y = ~F_R, 
          type = 'scatter', 
          mode = 'lines+markers', 
          line = list(color = 'orange')) %>%
  layout(title = "India fatality rate over time",
         xaxis = list(title = "Year"),
         yaxis = list(title = "Fatality Rate"))

line_plot_India_FR 
saveWidget(line_plot_India_FR, "India fatality rate.html")
