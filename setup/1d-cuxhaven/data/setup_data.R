# SPDX-FileCopyrightText: 2025 Helmholtz-Zentrum hereon GmbH
# SPDX-License-Identifier: CC0-1.0
# SPDX-FileContributor Ovidio Garcia-Oliva <ovidio.garcia@hereon.de>

################################################################################
## definition of local functions

get.year = function(...) format(as.Date(..., format = "%d/%m/%Y"), "%Y")
as.POSIXct = function(...) base::as.POSIXct(..., format = "%Y-%m-%d %H:%M:%S")

################################################################################
## loading observed values for temperature (temp) wind velocity (wind) and wind
## direction (dirw) to create forcing file

temp = read.delim("Cuxhaven_DWD!Lufttemperatur.txt", 
                  header = FALSE, 
                  comment.char = "#"
)

wind = read.delim("Cuxhaven_DWD!Windgeschwindigkeit.txt", 
                  header = FALSE, 
                  comment.char = "#"
)

dirw = read.delim("Cuxhaven_DWD!Windrichtung.txt", 
                  header = FALSE, 
                  comment.char = "#"
)

temp$V1 = as.POSIXct(temp$V1)
wind$V1 = as.POSIXct(wind$V1)
dirw$V1 = as.POSIXct(dirw$V1)

## transforming temperatures at 9m above ground to 2 meter values
## assuming linear temperature profile
h_station = 9
h_model = 2
T_gradient = -0.0065
bias = 5
temp$temp = temp$V2 + (h_model - h_station) * T_gradient + bias + rnorm(length(temp$V2), 0, 0.1)

## transforming wind speed at 9m above ground to 10 meter values
## assuming logarithmic wind profile
h_station = 9
h_model = 10
w_exponent = 0.14
wind$V3 = wind$V2 * (h_model / h_station)^w_exponent 

## merging wind files
wind = merge(wind, dirw, by = "V1", suffixes = c(".vel", ".dir"))
  
## calculating wind components
wind$u10 = wind$V3 * cos(pi * wind$V2.dir / 180) + rnorm(length(wind$V2 ), 0, 0.01)
wind$v10 = wind$V3 * sin(pi * wind$V2.dir / 180) + rnorm(length(wind$V2 ), 0, 0.01)

## creating the meteofile
meteofile = merge(temp, wind)
meteofile = subset(meteofile, 
                   get.year(V1) >= 2004 & !is.na(temp) & !is.na(v10) & !is.na(u10),
                   select=c(V1, temp, u10, v10)
)

meteofile$V1 = meteofile$V1 + 1 # adding 1 sec to avoid a bug with gotm config
meteofile = meteofile[get.year(meteofile$V1)%in%2019:2020,]

colnames(meteofile) = paste0("#", colnames(meteofile))
write.csv(meteofile, "meteofile.csv", row.names=FALSE ,quote=FALSE)
