
  library(tidyverse)
  library(lubridate)
  library(mgcv)
  library(corrplot)
  library(RColorBrewer)
  library(scales)
  library(ozmaps)
  library(cowplot)

# load previous analysis  
  load("Bugs GAM final.RData")
  
# data summaries
  glimpse(dat)
  
# number of sites
  length(table(dat$id))
# number of habitats sampled  
  colSums(table(dat$id, dat$SampleHabitat))

# one or both habitats at a site  
  table(rowSums(table(dat$id, dat$SampleHabitat)))

################################################################################  
# Figure S1
# scale_temp versus av_temp
  
  pdf("Figure S1.pdf", width = 7, height = 4)
  
  par(mfrow = c(1, 2))
  
  plot(scale_temp ~ av_temp, data = filter(dat, SampleHabitat == "Riffle"), bty = "l",
       xlab = "WorldClim average annual temperature", ylab = "Scaled temperature",
       main = "Riffle", pch = 19, col = rgb(0, 0, 0, 0.1), cex.lab = 0.8)
  mtext(bquote(bold((a))), adj = 0, xpd = NA, cex = 0.8)
  
  plot(scale_temp ~ av_temp, data = filter(dat, SampleHabitat == "Edge"), bty = "l",
       xlab = "WorldClim average annual temperature", ylab = "",
       main = "Edge", pch = 19, col = rgb(0, 0, 0, 0.1), cex.lab = 0.8)
  mtext(bquote(bold((b))), adj = 0, xpd = NA, cex = 0.8)
  
  dev.off()
  
  cor.test(dat$scale_temp[dat$SampleHabitat == "Riffle"], dat$av_temp[dat$SampleHabitat == "Riffle"])
  cor.test(dat$scale_temp[dat$SampleHabitat == "Edge"], dat$av_temp[dat$SampleHabitat == "Edge"])
  

################################################################################
# var1 = scale_temp
# var2 = slope
# var3 = log_turbidity
# var4 = log_ec
  
# summary of model fits
  
  r <- summary(mrif)
  e <- summary(medg)

  a <- as.data.frame(r$s.table)
  b <- as.data.frame(e$s.table)
  
  out.m <- data.frame(var = rownames(a),
                      r.edf = round(a[, 1], 1), 
                      r.shrink = ifelse(r$s.pv >0.05, "*", ""),
                      r.chi = round(a[, 3], 1),
                      r.p = round(a[, 4], 3),
                      
                      e.edf = round(b[, 1], 1),
                      e.shrink = ifelse(e$s.pv >0.05, "*", ""),
                      e.chi = round(b[, 3], 1),
                      e.p = round(b[, 4], 3))
  
  out.m$var <- gsub("var1", "temperature", out.m$var)
  out.m$var <- gsub("var2", "slope", out.m$var)
  out.m$var <- gsub("var3", "turbidity", out.m$var)
  out.m$var <- gsub("var4", "EC", out.m$var)
  
  out.m
  write.csv(out.m, "c:\\users\\s429217\\onedrive\\data\\bugs\\Model results.csv", row.names = FALSE)
  
##########################################################################################
# Figure 1
  table(dat$SampleHabitat)
  
  dat_text <- data.frame(
    label = c("n = 4339", "n = 2533"),
    SampleHabitat   = c("Edge", "Riffle"),
    x     = c(150, 150),
    y     = c(-10, -10))
  
  cols <- brewer.pal(n = 10, name = "RdBu") 
  cols <- rev(cols)
  
  state <- ozmap_states %>%
    filter(NAME %in% c("New South Wales", "Victoria", "Queensland", "Tasmania"))
  
  # ept richness
  p <- ggplot(data = state) +
    geom_sf(fill = "white") +
    geom_point(data = dat, aes(y = Latitude, x = Longitude, colour = ept_rich), size = 0.2) +
    scale_colour_gradientn(colours = cols) +
    labs(color = "EPT family richness") +
    facet_wrap(~ SampleHabitat) +
    geom_text(data = dat_text, mapping = aes(x = x, y = y, label = label)) +
    theme_bw() +
    theme(legend.position = "top",
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  
  p1 <- ggplot(dat, aes(x = ept_rich, after_stat(density))) +
    geom_histogram(binwidth = 1, fill = "grey80", col = "black") +
    facet_wrap(~ SampleHabitat) +
    xlab("EPT family richness") +
    ylab("Density") +
    theme_classic() +
    theme(strip.background = element_blank(),
          strip.text.x = element_blank())


  pdf("Figure 1.pdf", width = 7, height = 8)
  
  ggdraw() +
    draw_plot(p, x = 0, y = .15, width = 1, height = .85) +
    draw_plot(p1, x = 0.145, y = 0, width = 0.72, height = 0.15) +
    draw_plot_label(label = c("(a)", "(b)"), size = 15,
                  x = c(0.13, 0.13), y = c(0.97, 0.18))
  
  dev.off()
  
################################################################################
# Figure S4
  
  pdf("Figure S4.pdf", width = 7, height = 8)
  
  # scale_temp  
  ggplot(data = state) +
    geom_sf(fill = "white") +
    geom_point(data = dat, aes(y = Latitude, x = Longitude, colour = scale_temp), size = 0.5) +
    scale_colour_gradientn(colours = cols) +
    theme_classic() +
    theme(legend.position = "none")
  
  dev.off()
  
  
################################################################################
# univariate plots
# function for plotting
  plot.rel <- function(param, let1, let2, y.lab1 = "", y.lab2 = "", x.lab = "", xtick, xval, m = FALSE) {
    
    dat$var <- param
    dat$ept_rich <- round(dat$ept_rich)
    
    dat$svar <- c(scale(param))
    
    rif <- filter(dat, SampleHabitat == "Riffle")
    edg <- filter(dat, SampleHabitat == "Edge")
    
  # models for plotting fit against data
    
    m0r <- gam(ept_rich ~ s(var), data = rif, family = nb(), method = "ML")
    m0e <- gam(ept_rich ~ s(var), data = edg, family = nb(), method = "ML")
    
    # linear fit
    m1r <- gam(ept_rich ~ var, data = rif, family = nb(), method = "ML")
    m1e <- gam(ept_rich ~ var, data = edg, family = nb(), method = "ML")
    
    xr <- data.frame(var = seq(min(rif$var), max(rif$var), 0.01))
    mod <- predict(m0r, newdata = xr, type = "response", se.fit = TRUE)
    fitr <- data.frame(rich = mod$fit, se = mod$se.fit, 
                       var = xr$var)
    fitr$SampleHabitat <- "Riffle"
    fitr$lcl <- fitr$rich - 2*fitr$se
    fitr$ucl <- fitr$rich + 2*fitr$se
    
    xe <- data.frame(var = seq(min(edg$var), max(edg$var), 0.01))
    mod <- predict(m0e, newdata = xe, type = "response", se.fit = TRUE)
    fite <- data.frame(rich = mod$fit, se = mod$se.fit, 
                       var = xe$var)
    fite$SampleHabitat <- "Edge"
    fite$lcl <- fite$rich - 2*fite$se
    fite$ucl <- fite$rich + 2*fite$se
    
    fitd <- bind_rows(fitr, fite)
    
    temp <- dat %>%
      group_by(SampleHabitat) %>%
      mutate(tc = cut_interval(var, 20)) %>%
      group_by(SampleHabitat, tc) %>%
      summarise(temp = mean(var),
                sd_rich = sd(ept_rich),
                mean_rich = mean(ept_rich),
                n = n()) %>%
      mutate(se_rich = sd_rich / sqrt(n),
             lcl = mean_rich - 2*se_rich,
             ucl = mean_rich + 2*se_rich)
    
    yl <- range(dat$ept_rich)
    xl <- range(param)
    
    m1 <- ""
    m2 <- ""
    if(m == TRUE) {
      m1 <- "Riffle"
      m2 <- "Edge"
    }

  # models for comparing fit taking account of spatial autocorrelation
    
    m0r.sp <- gam(ept_rich ~ s(var) + s(Latitude, Longitude), data = rif, family = nb(), method = "ML")
    m0e.sp <- gam(ept_rich ~ s(var) + s(Latitude, Longitude), data = edg, family = nb(), method = "ML")
    
  # linear fit
    m1r.sp <- gam(ept_rich ~ var + s(Latitude, Longitude), data = rif, family = nb(), method = "ML")
    m1e.sp <- gam(ept_rich ~ var + s(Latitude, Longitude), data = edg, family = nb(), method = "ML")    
    
  # AIC comparison (positive values = non-linear is better fitting)      
    c1 <- round(AIC(m1r.sp) - AIC(m0r.sp), 1)
    c2 <- round(AIC(m1e.sp) - AIC(m0e.sp), 1)
    
  # scaled fit to get marginal effect
    m0r.sc <- gam(ept_rich ~ s(svar) + s(Latitude, Longitude), data = rif, family = nb(), method = "ML")
    m0e.sc <- gam(ept_rich ~ s(svar) + s(Latitude, Longitude), data = edg, family = nb(), method = "ML")
    
    pvar <- c(-1.6, 1.6)
    pred.dat <- data.frame(svar = pvar,
                           Latitude = rep(-30, 2),
                           Longitude = rep(150, 2))  
    
    predr <- predict.gam(m0r.sc, newdata = pred.dat, type = "response")  
    prede <- predict.gam(m0e.sc, newdata = pred.dat, type = "response")  
    pr.r <- round(1 - predr[2]/predr[1], 2)
    pr.e <- round(1 - prede[2]/prede[1], 2)
    
    
    plot(ept_rich ~ var, data = filter(dat, SampleHabitat == "Riffle"), pch = 19, 
         col = rgb(0, 0, 0, 0.02), ylim = yl, xlim = xl,
         ylab = y.lab1, xlab = x.lab, cex.lab = 1.2, bty = "l", xaxt = "n", main = m1, cex.main = 1.4)
    lines(rich ~ var, data = filter(fitd, SampleHabitat == "Riffle"), col = "blue", lwd = 2)
    lines(lcl ~ var, data = filter(fitd, SampleHabitat == "Riffle"), col = "blue", lwd = 2, lty = 3)
    lines(ucl ~ var, data = filter(fitd, SampleHabitat == "Riffle"), col = "blue", lwd = 2, lty = 3)
    points(mean_rich ~ temp, data = filter(temp, SampleHabitat == "Riffle"), pch = 19, col = "tomato", cex = 1.2)
    
    text(xl[2], 18, paste("AIC dif. = ", c1), cex = 1, xpd = NA, pos = 2)
    text(xl[2], 15.5, paste("Deviance expl. = ", round(summary(m0r.sp)$dev.expl, 2)), cex = 1, xpd = NA, pos = 2)
    text(xl[2], 13, paste("Reduction = ", pr.r), cex = 1, xpd = NA, pos = 2)
    mtext(bquote(bold((.(let1)))), adj = 0, xpd = NA, cex = 0.8)
    
    axis(1, at = xtick, labels = xval)
    
    plot(ept_rich ~ var, data = filter(dat, SampleHabitat == "Edge"), pch = 19, 
         col = rgb(0, 0, 0, 0.02), ylim = yl, xlim = xl,
         ylab = y.lab2, xlab = x.lab, cex.lab = 1.2, bty = "l", xaxt = "n", main = m2, cex.main = 1.4)
    lines(rich ~ var, data = filter(fitd, SampleHabitat == "Edge"), col = "blue", lwd = 2)
    lines(lcl ~ var, data = filter(fitd, SampleHabitat == "Edge"), col = "blue", lwd = 2, lty = 3)
    lines(ucl ~ var, data = filter(fitd, SampleHabitat == "Edge"), col = "blue", lwd = 2, lty = 3)
    points(mean_rich ~ temp, data = filter(temp, SampleHabitat == "Edge"), pch = 19, col = "tomato", cex = 1.2)
    
    text(xl[2], 18, paste("AIC dif. = ", c2), cex = 1, xpd = NA, pos = 2)
    text(xl[2], 15.5, paste("Deviance expl. = ", round(summary(m0e.sp)$dev.expl, 2)), cex = 1, xpd = NA, pos = 2)
    text(xl[2], 13, paste("Reduction = ", pr.e), cex = 1, xpd = NA, pos = 2)
    mtext(bquote(bold((.(let2)))), adj = 0, xpd = NA, cex = 0.8)
    
    axis(1, at = xtick, labels = xval)
    
    return(list(m1r.sp, m0r.sp, m1e.sp, m0e.sp))  
  }  
# end function

################################################################################
# Figure 2
  pdf("Figure 2.pdf", width = 7, height = 7)
  
  par(mfrow = c(4, 2), mar = c(4, 4, 2, 1))
  
  temp <- plot.rel(dat$scale_temp, "a", "b", y.lab1 = "EPT richness", x.lab = "Scaled temperature", 
           xtick = seq(-1, 1.5, 0.5), xval = seq(-1, 1.5, 0.5), m = TRUE)
  
  slope <- plot.rel(dat$slope, "c", "d", y.lab1 = "EPT richness", x.lab = "Slope (log(horizontal/vertical))",
           xtick = 2:7, xval = 2:7)  
  
  turb <- plot.rel(dat$log_turb, "e", "f", y.lab1 = "EPT richness", x.lab = "Turbidity (NTU)",
           xtick = log(c(0.1, 1, 10, 100, 1000)), xval = c(0.1, 1, 10, 100, 1000))
  
  ec <- plot.rel(dat$log_ec, "g", "h", y.lab1 = "EPT richness", x.lab = "Electircal conductivity (uS/cm)",
           xtick = log(c(1, 10, 100, 1000, 10000)), xval = c(1, 10, 100, 1000, 10000))
  
  dev.off()
  
###############################################################################  
# Figure S3
# observed versus predicted

  rif$pred <- predict(mrif, type = "response")
  edg$pred <- predict(medg, type = "response")

# Riffle
# plot observed versus predicted
  a <- tapply(rif$ept_rich, round(rif$pred), mean)

  pdf("Figure S3.pdf", width = 7, height = 4)
  
  par(mfrow = c(1, 2), mar = c(4, 4, 4, 4))
  plot(rif$ept_rich ~ rif$pred, pch = 19, col = rgb(0, 0, 0, 0.02), main = "Riffle",
       xlab = "Predictied richness", ylab = "Observed EPT family richness",
       ylim = c(0, 15), xlim = c(0, 15), bty = "l", cex.lab = 0.8)
    points(a ~ as.numeric(names(a)), col = "red", cex = 1.5, pch = 19)
  abline(0, 1, col = "blue", lwd = 2)
  text(0, 17, paste("Deviance explained = ", round(summary(mrif)$dev.expl, 2)), xpd = NA, pos = 4, cex = 0.8)
  mtext(bquote(bold((a))), adj = 0, xpd = NA, cex = 0.8, line = 2)
  
  cor(rif$pred, rif$ept_rich)

# Edge
# plot observed versus predicted
  a <- tapply(edg$ept_rich, round(edg$pred), mean)
  
  plot(edg$ept_rich ~ edg$pred, pch = 19, col = rgb(0, 0, 0, 0.02), main = "Edge",
       xlab = "Predictied richness", ylab = "Observed EPT family richness",
       ylim = c(0, 15), xlim = c(0, 15), bty = "l", cex.lab = 0.8)
  points(a ~ as.numeric(names(a)), col = "red", cex = 1.5, pch = 19)
  abline(0, 1, col = "blue", lwd = 2)
  text(0, 17, paste("Deviance explained = ", round(summary(medg)$dev.expl, 2)), xpd = NA, pos = 4, cex = 0.8)
  mtext(bquote(bold((b))), adj = 0, xpd = NA, cex = 0.8, line = 2)
  
  cor(edg$pred, edg$ept_rich)
  
  dev.off()
  
###############################################################################
# number of sites by variables
  
# function to draw box
  mep = 1.6
  dbox <- function(me = mep) {
    lines(x = c(-me, -me), y = c(-me, me), lwd = 2)
    lines(x = c(-me, me), y = c(me, me), lwd = 2)
    lines(x = c(me, me), y = c(-me, me), lwd = 2)
    lines(x = c(-me, me), y = c(-me, -me), lwd = 2)
  }
  
  x <- seq(-2.5, 2.5, 0.5)
  y <- seq(-2.5, 2.5, 0.5)
  
  rx <- cut(scale(rif$scale_temp), breaks = x)
  ry <- cut(scale(rif$slope), breaks = y)
  tr <- table(rx, ry)  
  tr <- ifelse(tr > 0 & tr <= 10, 1, tr)
  tr <- ifelse(tr > 10 & tr <= 50, 2, tr)
  tr <- ifelse(tr > 50 & tr <= 100, 3, tr)
  tr <- ifelse(tr > 100, 4, tr)
  
  ex <- cut(scale(edg$scale_temp), breaks = x)
  ey <- cut(scale(edg$slope), breaks = y)
  te <- table(ex, ey)  
  te <- ifelse(te > 0 & te <= 10, 1, te)
  te <- ifelse(te > 10 & te <= 50, 2, te)
  te <- ifelse(te > 50 & te <= 100, 3, te)
  te <- ifelse(te > 100, 4, te)
  
  
  cl <- hcl.colors(5, "YlOrRd", rev = TRUE)
  cl[1] <- "#FFFFFF"
  
  pdf("Figure S2.pdf", width = 7, height = 7)

  par(mfrow = c(2, 2))
  image(x = x, y = y, z = tr, ylab = "Slope", xlab = "Temperature", col = cl)
  mtext(expression(bold(Riffle)), line = 1, adj = 0)
  mtext(bquote(bold((a))), adj = -0.2, xpd = NA, cex = 0.8)
  dbox()
  dbox(2.5)

  legend(x = 0, y = 4, legend = c("0", "1-10", "11-50", "51-100", ">100"),
         fill = cl, horiz = TRUE, xpd = NA, cex = 0.8, title = "Number of sites sampled")

  image(x = x, y = y, z = te, ylab = "", xlab = "Temperature", col = cl)
  mtext(expression(bold(Edge)), line = 1, adj = 1)
  mtext(bquote(bold((b))), adj = -0.2, xpd = NA, cex = 0.8)
  dbox()
  dbox(2.5)
  
  
  rx <- cut(scale(rif$log_turb), breaks = x)
  ry <- cut(scale(rif$log_ec), breaks = y)
  tr <- table(rx, ry)  
  tr <- ifelse(tr > 0 & tr <= 10, 1, tr)
  tr <- ifelse(tr > 10 & tr <= 50, 2, tr)
  tr <- ifelse(tr > 50 & tr <= 100, 3, tr)
  tr <- ifelse(tr > 100, 4, tr)
  
  ex <- cut(scale(edg$log_turb), breaks = x)
  ey <- cut(scale(edg$log_ec), breaks = y)
  te <- table(ex, ey)  
  te <- ifelse(te > 0 & te <= 10, 1, te)
  te <- ifelse(te > 10 & te <= 50, 2, te)
  te <- ifelse(te > 50 & te <= 100, 3, te)
  te <- ifelse(te > 100, 4, te)
  
  image(x = x, y = y, z = tr, ylab = "EC", xlab = "Turbidity", col = cl)
  mtext(bquote(bold((c))), adj = -0.2, xpd = NA, cex = 0.8)
  dbox()
  dbox(2.5)

  image(x = x, y = y, z = te, ylab = "", xlab = "Turbidity", col = cl)
  mtext(bquote(bold((d))), adj = -0.2, xpd = NA, cex = 0.8)
  dbox()
  dbox(2.5)
  
  dev.off()
  
################################################################################################
# plot the interactions
# me specifies the extremes to plot
  
  # var1 = temp
  # var2 = slope
  # var3 = log_turbidity
  # var4 = log_ec
  
# 2 dimension plots
# function to plot - specify x axis variable (xv) and interaction variable (iv)

# specify max/min value for the extremes
  me <- 1.6
  
  val <- round(seq(-me, me, 0.1), 1)
  n <- length(val)
  yxl <- c(0, 15)

# start function  
  pl <- function(let, v1, v2, v3, v4, xv, iv, mod, yaxlim = yxl,
                 lab.x = TRUE, lab.y = TRUE) {

  if(xv == 1) xvar <- v1
  if(xv == 2) xvar <- v2
  if(xv == 3) xvar <- v3
  if(xv == 4) xvar <- v4
  
  if(iv == 1) ivar <- v1 
  if(iv == 2) ivar <- v2 
  if(iv == 3) ivar <- v3
  if(iv == 4) ivar <- v4
  
  if(xv == 1) xl <- "Temperature"
  if(xv == 2) xl <- "Slope"
  if(xv == 3) xl <- "Turbidity"
  if(xv == 4) xl <- "EC"
  
  if(lab.x == FALSE) xl = ""
  
  pred.dat <- data.frame(var1 = v1,
                         var2 = v2,
                         var3 = v3,
                         var4 = v4,
                         Latitude = rep(-30, n*3),
                         Longitude = rep(150, n*3))  
  
  pred <- predict.gam(mod, newdata = pred.dat, type = "response", se.fit = TRUE)
  pred.dat$ept <- pred$fit
  pred.dat$se <- pred$se.fit
  
  pred.dat <- pred.dat %>%
    mutate(lcl = ept - 2*se,
           ucl = ept + 2*se)
  
  cols <- t(col2rgb(lc <- brewer.pal(n = 3, name = "Accent")))
  cols <- cols/255

  
    pd <- filter(pred.dat, ivar == -me)
    cl <- rgb(cols[1,1], cols[1,2], cols[1,3], 0.1)

    yl <- "EPT family richness"
    if(lab.y == FALSE) yl <- ""
    
    plot(ept ~ val, data = pd, type = "n", bty = "l", ylim = yaxlim, xlab = xl, ylab = yl)
    mtext(bquote(bold((.(let)))), adj = 0, xpd = NA, cex = 0.8)

    polygon(x = c(val, rev(val)), y = c(pd$lcl, rev(pd$ucl)), 
            col = cl, border = NA)
    lines(pd$ept ~ val, col = lc[1], lwd = 3)
    
    pd <- filter(pred.dat, ivar == 0)
    cl <- rgb(cols[2,1], cols[2,2], cols[2,3], 0.1)
    polygon(x = c(val, rev(val)), y = c(pd$lcl, rev(pd$ucl)), 
            col = cl, border = NA)
    lines(pd$ept ~ val, col = lc[2], lwd = 3)
    
    pd <- filter(pred.dat, ivar == me)
    cl <- rgb(cols[3,1], cols[3,2], cols[3,3], 0.1)
    polygon(x = c(val, rev(val)), y = c(pd$lcl, rev(pd$ucl)), 
            col = cl, border = NA)
    lines(pd$ept ~ val, col = lc[3], lwd = 3)
    
  # measure of the strength of interaction = ratio of max to min difference
    pd1 <- filter(pred.dat, ivar == me)          # high (worst) condition
    pd2 <- filter(pred.dat, ivar == -me)         # low (best) condition
    max.dif <- abs(max(pd2$ept - pd1$ept))
    min.dif <- abs(min(pd2$ept - pd1$ept))
    
    dif <- abs(round(max.dif / min.dif, 1))
    
  # dif at extremes
    e1 <- abs(pred.dat$ept[ivar == me & xvar == me] - pred.dat$ept[ivar == -me & xvar == me])
    e2 <- abs(pred.dat$ept[ivar == me & xvar == -me] - pred.dat$ept[ivar == -me & xvar == -me])
    emid <- abs(pred.dat$ept[ivar == me & xvar == 0] - pred.dat$ept[ivar == -me & xvar == 0])

  # define type
    type <- "Additive"
    if(dif >= 2 & (e2 > e1)) type = "Antagonistic"
    if(dif >= 2 & (e2 < e1)) type = "Synergistic"
    if(dif >= 2 & (emid > e2 & emid > e1)) type = "Antagonistic"

    text(me, yaxlim[2], type, xpd = NA, cex = 1, pos = 2)

  # another measure to evaluate consistency
  # proportional reduction in richness from low to high stress due to cumulative effects
    exdif <- 1 - round(pred.dat$ept[ivar == me & xvar == me] / pred.dat$ept[ivar == -me & xvar == -me], 2)
    text(me, yaxlim[2] - 1, paste("C =", exdif), xpd = NA, cex = 1, pos = 2)
    
  # marginal effects
  # reduction from low to high averaged across all values of the other variable for EC and turbidity
  # EC
    
    x <- seq(-me, me, 0.01)
    n <- length(x)
    
    pred.dat <- data.frame(var1 = rep(v1[1], n),
                           var2 = rep(v2[1], n),
                           Latitude = rep(-30, n),
                           Longitude = rep(150, n),
                           var3 = x,
                           var4 = rep(-me, n))      
    
    pred1 <- predict.gam(mod, newdata = pred.dat, type = "response", se.fit = TRUE)
    
    pred.dat <- data.frame(var1 = rep(v1[1], n),
                           var2 = rep(v2[1], n),
                           Latitude = rep(-30, n),
                           Longitude = rep(150, n),
                           var3 = x,
                           var4 = rep(me, n))      
    
    pred2 <- predict.gam(mod, newdata = pred.dat, type = "response", se.fit = TRUE)

    exdif1 <- 1 - round(mean(pred2$fit)/mean(pred1$fit), 2)
    text(me, yaxlim[2] - 2, paste("E =", exdif1), xpd = NA, cex = 1, pos = 2)

    # Turbidity
    
    pred.dat <- data.frame(var1 = rep(v1[1], n),
                           var2 = rep(v2[1], n),
                           Latitude = rep(-30, n),
                           Longitude = rep(150, n),
                           var3 = rep(-me, n),
                           var4 = x)      
    
    pred3 <- predict.gam(mod, newdata = pred.dat, type = "response", se.fit = TRUE)
    
    pred.dat <- data.frame(var1 = rep(v1[1], n),
                           var2 = rep(v2[1], n),
                           Latitude = rep(-30, n),
                           Longitude = rep(150, n),
                           var3 = rep(me, n),
                           var4 = x)      
    
    pred4 <- predict.gam(mod, newdata = pred.dat, type = "response", se.fit = TRUE)
    
    exdif2 <- 1 - round(mean(pred4$fit)/mean(pred3$fit), 2)
    text(me, yaxlim[2] - 3, paste("T =", exdif2), xpd = NA, cex = 1, pos = 2)
    
    # ratio of additive mariginal effects to cumulative effect
    rat <- exdif / (exdif1 + exdif2) 
#    text(me, yaxlim[2], round(rat, 2), xpd = NA, cex = 1, pos = 2)
  }
# end function
  
###################################################################################
# RIFFLE 
  pdf("Figure 3.pdf", width = 7, height = 7)
  
  yxl <- c(2, 12)
  par(mfrow = c(3, 3), mar = c(4,4,2,0), oma = c(0, 0, 5, 5))
  
# first row = steep slope (-2)  
  
  pl("a", v1 = rep(-me, n*3), v2 = rep(-me, n*3), v3 = rep(val, 3), v4 = rep(c(-me, 0, me), each = n),
     xv = 3, iv = 4, mod = mrif, lab.x = FALSE)
  pl("b", v1 = rep(0, n*3), v2 = rep(-me, n*3), v3 = rep(val, 3), v4 = rep(c(-me, 0, me), each = n),
     xv = 3, iv = 4, mod = mrif, lab.y = FALSE, lab.x = FALSE)
  pl("c", v1 = rep(me, n*3), v2 = rep(-me, n*3), v3 = rep(val, 3), v4 = rep(c(-me, 0, me), each = n),
     xv = 3, iv = 4, mod = mrif, lab.y = FALSE, lab.x = FALSE)
  

# second row = mid slope (0)  
  pl("d", v1 = rep(-me, n*3), v2 = rep(0, n*3), v3 = rep(val, 3), v4 = rep(c(-me, 0, me), each = n),
     xv = 3, iv = 4, mod = mrif, lab.x = FALSE)
  pl("e", v1 = rep(0, n*3), v2 = rep(0, n*3), v3 = rep(val, 3), v4 = rep(c(-me, 0, me), each = n),
     xv = 3, iv = 4, mod = mrif, lab.x = FALSE, lab.y = FALSE)
  pl("f", v1 = rep(me, n*3), v2 = rep(0, n*3), v3 = rep(val, 3), v4 = rep(c(-me, 0, me), each = n),
     xv = 3, iv = 4, mod = mrif, lab.x = FALSE, lab.y = FALSE)
  
# third row = mid slope (0)  
  pl("g", v1 = rep(-me, n*3), v2 = rep(me, n*3), v3 = rep(val, 3), v4 = rep(c(-me, 0, me), each = n),
     xv = 3, iv = 4, mod = mrif)
  pl("h", v1 = rep(0, n*3), v2 = rep(me, n*3), v3 = rep(val, 3), v4 = rep(c(-me, 0, me), each = n),
     xv = 3, iv = 4, mod = mrif, lab.y = FALSE)
  pl("i", v1 = rep(me, n*3), v2 = rep(me, n*3), v3 = rep(val, 3), v4 = rep(c(-me, 0, me), each = n),
     xv = 3, iv = 4, mod = mrif, lab.y = FALSE)
  
  text(-4.5, 52, "Temperature", xpd = NA, cex = 2)
  text(-9, 50.5, "Cool", xpd = NA, cex = 1.8)
  text(-4.5, 50.5, "Intermediate", xpd = NA, cex = 1.8)
  text(0, 50.5, "Warm", xpd = NA, cex = 1.8)
  
  text(2.6, 24, "Slope", xpd = NA, srt = 270, cex = 2)
  text(2.2, 24, "Intermediate", xpd = NA, srt = 270, cex = 1.8)
  text(2.2, 40, "Steep", xpd = NA, srt = 270, cex = 1.8)
  text(2.2, 6, "Shallow", xpd = NA, srt = 270, cex = 1.8)
  
  dev.off()
  
########################################################################################  
# EDGE
  pdf("Figure 4.pdf", width = 7, height = 7)
  
  yxl <- c(2, 10)
  par(mfrow = c(3, 3), mar = c(4,4,2,0), oma = c(0, 0, 5, 5))
  
# first row = steep slope (-2)  
  pl("a", v1 = rep(-me, n*3), v2 = rep(-me, n*3), v3 = rep(val, 3), v4 = rep(c(-me, 0, me), each = n),
     xv = 3, iv = 4, mod = medg, lab.x = FALSE)
  pl("b", v1 = rep(0, n*3), v2 = rep(-me, n*3), v3 = rep(val, 3), v4 = rep(c(-me, 0, me), each = n),
     xv = 3, iv = 4, mod = medg, lab.x = FALSE, lab.y = FALSE)
  pl("c", v1 = rep(me, n*3), v2 = rep(-me, n*3), v3 = rep(val, 3), v4 = rep(c(-me, 0, me), each = n),
     xv = 3, iv = 4, mod = medg, lab.x = FALSE, lab.y = FALSE)
  
  
# second row = mid slope (0)  
  pl("d", v1 = rep(-me, n*3), v2 = rep(0, n*3), v3 = rep(val, 3), v4 = rep(c(-me, 0, me), each = n),
     xv = 3, iv = 4, mod = medg, lab.x = FALSE)
  pl("e", v1 = rep(0, n*3), v2 = rep(0, n*3), v3 = rep(val, 3), v4 = rep(c(-me, 0, me), each = n),
     xv = 3, iv = 4, mod = medg, lab.x = FALSE, lab.y = FALSE)
  pl("f", v1 = rep(me, n*3), v2 = rep(0, n*3), v3 = rep(val, 3), v4 = rep(c(-me, 0, me), each = n),
     xv = 3, iv = 4, mod = medg, lab.x = FALSE, lab.y = FALSE)
  
# third row = mid slope (0)  
  pl("g", v1 = rep(-me, n*3), v2 = rep(me, n*3), v3 = rep(val, 3), v4 = rep(c(-me, 0, me), each = n),
     xv = 3, iv = 4, mod = medg)
  pl("h", v1 = rep(0, n*3), v2 = rep(me, n*3), v3 = rep(val, 3), v4 = rep(c(-me, 0, me), each = n),
     xv = 3, iv = 4, mod = medg, lab.y = FALSE)
  pl("i", v1 = rep(me, n*3), v2 = rep(me, n*3), v3 = rep(val, 3), v4 = rep(c(-me, 0, me), each = n),
     xv = 3, iv = 4, mod = medg, lab.y = FALSE)
  
  text(-4.5, 42, "Temperature", xpd = NA, cex = 2)
  text(-9, 40.5, "Cool", xpd = NA, cex = 1.8)
  text(-4.5, 40.5, "Intermediate", xpd = NA, cex = 1.8)
  text(0, 40.5, "Warm", xpd = NA, cex = 1.8)
  
  text(2.6, 20, "Slope", xpd = NA, srt = 270, cex = 2)
  text(2.2, 20, "Intermediate", xpd = NA, srt = 270, cex = 1.8)
  text(2.2, 34, "Steep", xpd = NA, srt = 270, cex = 1.8)
  text(2.2, 6, "Shallow", xpd = NA, srt = 270, cex = 1.8)
  
  dev.off()  

############################################################################################
# marginal effects
# reduction from low to high when others are held constant

# EC
# the nine values for temp and slope  
  me <- 1.6
  
  x <- seq(-me, me, 0.01)
  n <- length(x)
  
  pred.dat <- data.frame(var1 = rep(-me, n),
                         var2 = rep(-me, n),
                         Latitude = rep(-30, n),
                         Longitude = rep(150, n),
                         var3 = x,
                         var4 = -me)  
  
  pred1 <- predict.gam(mrif, newdata = pred.dat, type = "response", se.fit = TRUE)

  pred.dat <- data.frame(var1 = rep(-me, n),
                         var2 = rep(-me, n),
                         Latitude = rep(-30, n),
                         Longitude = rep(150, n),
                         var3 = x,
                         var4 = me)  
  
  pred2 <- predict.gam(mrif, newdata = pred.dat, type = "response", se.fit = TRUE)
  
  1 - mean(pred2$fit)/mean(pred1$fit)
  
