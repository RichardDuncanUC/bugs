
  library(tidyverse)
  library(lubridate)
  library(mgcv)
  library(corrplot)
  library(RColorBrewer)
  library(scales)
  library(ozmaps)
  library(parallel)


# read in all data
  dat <- read.csv("All data final.csv")
  glimpse(dat)
  
# slope values are in radians
# convert to log(horizontal distance / vertical distance)
  dat$slope <- log(1/tan(dat$slope + 0.001))

# convert WorldClim temps to degress
  dat$av_temp <- dat$av_temp/10
  
################################################################################
# data summaries
# number of sites
  length(table(dat$id))
    
# number of sites for each habitat
  colSums(table(dat$id, dat$SampleHabitat))

# number of entries per site  
  table(table(dat$id))
  
################################################################################
# plot ept richness
  
  cols <- brewer.pal(n = 10, name = "RdBu") 
  cols <- rev(cols)
  
  state <- ozmap_states %>%
    filter(NAME %in% c("New South Wales", "Victoria", "Queensland", "Tasmania"))
  
  # ept richness
  ggplot(data = state) +
    geom_sf(fill = "white") +
    geom_point(data = dat, aes(y = Latitude, x = Longitude, colour = ept_rich), size = 0.8) +
    scale_colour_gradientn(colours = cols) +
    facet_wrap(~ SampleHabitat) +
    theme_classic() +
    theme(legend.position = "top")
  
################################################################################
# correlations among variables
  
  par(mfrow = c(1, 1))
  mm <- dat %>%
    ungroup() %>%
    dplyr::select(log_ec, log_turb, scale_temp, av_temp, slope, ept_rich) 

  mmm <- cor(as.matrix(mm))
  corrplot(mmm, is.corr = T, diag = F, type = "upper", method = "number")
  
##############################################################################
# fit GAM models to full data set

# extract variables  
  dat$var1 <- c(scale(dat$scale_temp))
  dat$var2 <- c(scale(dat$slope))
  dat$var3 <- c(scale(dat$log_turb))
  dat$var4 <- c(scale(dat$log_ec))
  dat$e <- round(dat$ept_rich)

# data for each habitat  
  rif <- filter(dat, SampleHabitat == "Riffle")
  edg <- filter(dat, SampleHabitat == "Edge")
  

################################################################################  
# Riffle model  
  mrif <- gam(e ~ s(var1) + s(var2) + s(var3) + s(var4) +
             ti(var1, var2) + ti(var1, var3) + ti(var1, var4) + ti(var2, var3) +
             ti(var2, var4) + ti(var3, var4) +
             ti(var1, var2, var3) + ti(var1, var3, var4) + ti(var2, var3, var4) +
             ti(var1, var2, var3, var4) +
             s(Longitude, Latitude),
             family = nb(), method = "ML", data = rif, select = TRUE, control = list(nthreads = 4))  
  
################################################################################  
# Edge model  
  medg <- gam(e ~ s(var1) + s(var2) + s(var3) + s(var4) +
             ti(var1, var2) + ti(var1, var3) + ti(var1, var4) + ti(var2, var3) +
             ti(var2, var4) + ti(var3, var4) +
             ti(var1, var2, var3) + ti(var1, var3, var4) + ti(var2, var3, var4) +
             ti(var1, var2, var3, var4) +
             s(Longitude, Latitude),
             family = nb(), method = "ML", data = edg, select = TRUE, control = list(nthreads = 4))
    
###############################################################################  
# save the results of model fitting  
  save.image(file = "Bugs GAM final.RData")
  
