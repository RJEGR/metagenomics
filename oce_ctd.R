# library(devtools)
library(oce)
# https://github.com/dankelley/oce
# https://dankelley.github.io/oce/articles/ctd.html
# https://github.com/RJEGR/BM/blob/master/Procesando_datos_CTD.md

# install_github("dankelley/oce", ref="develop") # bugs
# install.packages("oce")

x <- read.oce(system.file("extdata", "ctd.cnv", package="oce"))
plot(x) # summary with TS and profiles
plotTS(x)
plotScan(x)

setwd("/Users/cigom/metagenomics/CTD/Lances_XIXIMI6")

lan <- read.ctd.sbe("lan023b.cnv", 
                    missingValue = NULL,
                    columns=list(
                      Fluor= list(name = "haardtY")))


summary(lan)

names(lan@data)
str(lan@data$Fluor)

plot(lan)

plotTS(lan, col = "gray")
plotProfile(lan, which="density+N2")
plotProfile(lan, xtype = "fluorescence+oxygen",
            ytype=="pressure")
#



library(oce)
# http://cchdo.ucsd.edu/data/7971/ar18_58JH19941029_ct1.zip
# setwd("~/Downloads/ar18_58JH19941029_ct1")
files <- system("ls *.cnv", intern=TRUE)

for (i in 1:length(files)) {
  x <- read.ctd.sbe(files[i],
                    missingValue = NULL,
                    columns=list(
                      Fluor= list(name = "haardtY")))
  if (i == 1) {
    plotTS(x, type='o')
  } else {
    points(x[["salinity"]], x[["potential temperature"]])
    lines(x[["salinity"]], x[["potential temperature"]])
  }
}

# 

n <- length(files)
ctds <- vector("list", n) # to hold the CTD objects
station <- vector("list", n)

for (i in 1:n) {
  ctds[[i]] <- read.ctd.sbe(files[i],
                            missingValue = NULL,
                            columns=list(
                              Fluor= list(name = "haardtY")))
  station[[i]] <- ctds[[i]][["station"]]
}

S <- unlist(lapply(1:n, function(i) ctds[[i]][["salinity"]]))
T <- unlist(lapply(1:n, function(i) ctds[[i]][["temperature"]]))
p <- unlist(lapply(1:n, function(i) ctds[[i]][["pressure"]]))
overall <- as.ctd(S, T, p)
# png("ar18_%02d.png")
for (i in 1:n) {
  plotTS(overall, col='gray')
  lines(ctds[[i]][["salinity"]], ctds[[i]][["potential temperature"]])
  mtext(station[i], side=3, line=0)
}



#
library(tidyverse)

ctd.tb = lan@data%>%
  as_tibble()

temp = ggplot(data = ctd.tb%>%na.omit(), 
              aes(x = temperature, y = pressure))+
  geom_path( col = "red")+
  scale_y_reverse()+
  scale_x_continuous(position = "top")+
  theme_bw()+
  theme(axis.text = element_text(size = 12, colour = 1),
        axis.title = element_text(size = 14, colour = 1))+
  labs(x = expression(~Temperature~(degree~C)), y = "Pressure[dbar]")


salinity = ggplot(data = ctd.tb%>%na.omit(), 
                  aes(x = salinity, y = pressure))+
  geom_path( col = "darkgreen")+
  scale_y_reverse()+
  scale_x_continuous(position = "top")+
  theme_bw()+
  theme(axis.text = element_text(size = 12, colour = 1),
        axis.title = element_text(size = 14, colour = 1))+
  labs(x = expression(~Practical~salinity),
       y = expression(~Pressure~(dbar)))

oxygen = ggplot(data = ctd.tb%>%na.omit(), 
                aes(x = oxygen, y = pressure))+
  geom_path( col = "blue")+
  scale_y_reverse()+
  scale_x_continuous(position = "top")+
  theme_bw()+
  theme(axis.text = element_text(size = 12, colour = 1),
        axis.title = element_text(size = 14, colour = 1))+
  labs(x = expression(~DO~(mgL^{-3})),
       y = expression(~Pressure~(dbar)))

fluor = ggplot(data = ctd.tb%>%na.omit(), 
                aes(x = fluorescence, y = pressure))+
  geom_path( col = "black")+
  scale_y_reverse()+
  scale_x_continuous(position = "top")+
  theme_bw()+
  theme(axis.text = element_text(size = 12, colour = 1),
        axis.title = element_text(size = 14, colour = 1))+
  labs(x = expression(~Fluorescence~(mg/m^{3})),
       y = expression(~Pressure~(dbar)))

cowplot::plot_grid(temp, salinity, 
                   oxygen,
                   fluor, nrow = 1)
