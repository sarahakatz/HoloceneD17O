## SCRIPT S1: DETERMINISTIC LAKE WATER ISOTOPE MASS BALANCE MODEL
## "Holocene temperature and water stress in the Peruvian Andes: insights from lake carbonate clumped and triple oxygen isotopes"
## submitted to AGU Paleoceanography and Paleoclimatology
## Sarah A. Katz*, Naomi E. Levin, Mark B. Abbott, Donald T. Rodbell, Benjamin H. Passey, Nicole M. DeLuca, Darren J. Larsen, Arielle Woods
## Correspondence: skatzees@umich.edu (SAK)
## Updated: Dec 12, 2023

## INSTALL AND ATTACH PACKAGES
# install.packages("ggplot2")
# install.packages("ggpubr")
# install.packages("rgl")
library(ggplot2)
library(ggpubr)
library(rgl)

## PLOT PATH
plot.path <- "~/Desktop/"             ## user may update plot path

## STEADY STATE ISOTOPE LAKE BALANCE EQUATIONS (SECTIONS 3.1 AND TABLE 2)
## Rw = (aeq*Ri*(adiff(a-h) + h*(1-F)) + aeq*h*Xe*Rv*F) / (Xe + aeq*(1-Xe)*(adiff*(1-h) + h(1-F))
## Open basin lake where evaporated water contributes to atmospheric vapor (Benson and White, 1994; Passey and Ji, 2019, Eq. 6)

## SET UP
## CONSTANTS
    R18smow = 0.0020052   ## Baertschi, 1976; IAEA Reference sheet
    R17smow = 0.0003799   ## Li et al., 1988; IAEA Reference sheet
    theta.eq = 0.529      ## Barkan and Luz, 2005
    theta.diff = 0.5185   ## Barkan and Luz, 2007
    theta.ref = 0.528
    diffratio18 = 1/0.9723  ## Merlivat, 1978

## USER DEFINED VARIABLES
    temp = 14       ## Lake temperature in degrees Celsius
    Phi = 0.5       ## The relative proportion of diffusive (0; molecular diffusion) to turbulent (1; non-fractionating) transport of water vapor during evaporation
    h = c(rep(seq(0.1, 0.9, by=0.1), 11))   ## relative humidity at lake surface
    Xe = c(rep(0.0, 9),      ## Volumetric proportion of evaporation (E) to inputs (I). Xe = E/I. Ranges from open basin lakes (0), to closed basin lake (1.0)
           rep(0.1, 9), 
           rep(0.2, 9),
           rep(0.3, 9),
           rep(0.4, 9),
           rep(0.5, 9),
           rep(0.6, 9),
           rep(0.7, 9),
           rep(0.8, 9),
           rep(0.9, 9),
           rep(1.0, 9))
    
    f = 0.9   ## Fraction of atmospheric vapor derived from distal sources versus the lake itself
              ## Where 1 and 0 represent exclusively distal and lake-derived sources of atmospheric vapor, respectively. 

    ## INPUT WATER  
    ## Based on amount-weighted mean annual precipitation at Junin, Peru (Katz et al., 2023, EPSL). All in units of per mil.
        D17Oi = 0.031
        dp18Oi = -14.1
        dp17Oi = D17Oi + (dp18Oi * theta.ref) 

## Calculate R values for input waters
    Ri18 = exp(dp18Oi/1000)*R18smow
    Ri17 = exp(dp17Oi/1000)*R17smow

## Isotopic ratio of water that atmospheric vapor is in equilibrium with (estimated from input water values)
    d18Ov = -15     ## water that vapor is in equilibrium with
    D17Ov = 0.035   ## water that vapor is in equilibrium with

    dp18Ov = log(d18Ov/1000+1)*1000                   ## convert vapor d18O to d'18O
    dp17Ov = D17Ov + theta.ref*dp18Ov                 ## calculate vapor d17O to d'17O

## Temperature dependent equilibrium fractionation factor between vapor and liquid water
    aeq18 = exp((1137/((temp + 273.15)^2)) - (0.4156/(temp+273.15)) - 0.0020667)    ## Majoube 1971
    aeq17 = exp(theta.eq * log(aeq18))   

## Calculate R values for vapor
    Rv18 = (exp(dp18Ov/1000)* R18smow)/aeq18
    Rv17 = (exp(dp17Ov/1000)*R17smow)/aeq17

## Diffusion vs. pure turbulence (i.e. no fractionation). When Phi = 1, 100% diffusive fractionation; when Phi = 0, 0% diffusive fractionation (all turbulent)
    adiff18 = Phi*diffratio18 + (1-Phi)
    adiff17 = exp(theta.diff*log(adiff18))


## SET UP FOR LOOP
    dat = data.frame(cbind(Xe, h))
    dp18Ow <- vector()            ## Create empty vectors to hold products from for loop
    d18Ow <- vector()
    dp17Ow <- vector()
    Dp17Ow <- vector()

## OPEN FOR LOOP
    for (i in 1:nrow(dat)){

##  LAKE WATER CALCULATIONS
## Calculate R values for lake waters
# Rw18 = ((aeq18*Ri18*(adiff18*(1-h)+h*(1-f)))+(aeq18*h*Xe*Rv18*f))/
#   (Xe+aeq18*(1-Xe)*(adiff18*(1-h)+h*(1-f)))
# 
# Rw17 = ((aeq17*Ri17*(adiff17*(1-h)+h*(1-f)))+(aeq17*h*Xe*Rv17*f))/
#   (Xe+aeq17*(1-Xe)*(adiff17*(1-h)+h*(1-f)))

    Rw18 = ((aeq18*Ri18*(adiff18*(1-dat[i,2])+dat[i,2]*(1-f)))+(aeq18*dat[i,2]*dat[i,1]*Rv18*f))/
      (dat[i,1]+aeq18*(1-dat[i,1])*(adiff18*(1-dat[i,2])+dat[i,2]*(1-f)))
    
    Rw17 = ((aeq17*Ri17*(adiff17*(1-dat[i,2])+dat[i,2]*(1-f)))+(aeq17*dat[i,2]*dat[i,1]*Rv17*f))/
      (dat[i,1]+aeq17*(1-dat[i,1])*(adiff17*(1-dat[i,2])+dat[i,2]*(1-f)))

## Calculate delta (d) and delta prime (dp) values for lake waters in units of per mil and D'17O (Dp) in units of per meg.
    dp18Ow. = (log(Rw18/R18smow))*1000
    d18Ow. = ((Rw18/R18smow)-1)*1000
    dp17Ow. = (log(Rw17/R17smow))*1000
    
    Dp17Ow. = (dp17Ow. - (theta.ref*dp18Ow.))*1000

## Fill empty vectors with newly calculated values
    dp18Ow[i] = dp18Ow.
    d18Ow[i] = d18Ow.
    dp17Ow[i] = dp17Ow.
    Dp17Ow[i] = Dp17Ow.

}

    dat = data.frame(cbind(Xe, h, Phi, f, D17Oi, dp18Oi, temp, dp18Ow, d18Ow, dp17Ow, Dp17Ow))
    
    dat

###########
## PLOTS ##
###########

Xe0 <- subset(dat, dat[,1] == 0)
Xe0.1 <- subset(dat, dat[,1] == 0.1)
Xe0.2 <- subset(dat, dat[,1] == 0.2)
Xe0.3 <- subset(dat, dat[,1] == 0.3)
Xe0.4 <- subset(dat, dat[,1] == 0.4)
Xe0.5 <- subset(dat, dat[,1] == 0.5)
Xe0.6 <- subset(dat, dat[,1] == 0.6)
Xe0.7 <- subset(dat, dat[,1] == 0.7)
Xe0.8 <- subset(dat, dat[,1] == 0.8)
Xe0.9 <- subset(dat, dat[,1] == 0.9)
Xe1 <- subset(dat, dat[,1] == 1)


Fig6 <- ggplot()+
  theme_bw()+
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())+
  geom_line(aes(x=Xe0$h, y=Xe0$Dp17Ow-Xe0$D17Oi*1000))+
  geom_line(aes(x=Xe0.1$h, y=Xe0.1$Dp17Ow-Xe0$D17Oi*1000))+
  geom_line(aes(x=Xe0.2$h, y=Xe0.2$Dp17Ow-Xe0$D17Oi*1000))+
  geom_line(aes(x=Xe0.3$h, y=Xe0.3$Dp17Ow-Xe0$D17Oi*1000))+
  geom_line(aes(x=Xe0.4$h, y=Xe0.4$Dp17Ow-Xe0$D17Oi*1000))+
  geom_line(aes(x=Xe0.5$h, y=Xe0.5$Dp17Ow-Xe0$D17Oi*1000))+
  geom_line(aes(x=Xe0.6$h, y=Xe0.6$Dp17Ow-Xe0$D17Oi*1000))+
  geom_line(aes(x=Xe0.7$h, y=Xe0.7$Dp17Ow-Xe0$D17Oi*1000))+
  geom_line(aes(x=Xe0.8$h, y=Xe0.8$Dp17Ow-Xe0$D17Oi*1000))+
  geom_line(aes(x=Xe0.9$h, y=Xe0.9$Dp17Ow-Xe0$D17Oi*1000))+
  geom_line(aes(x=Xe1$h, y=Xe1$Dp17Ow-Xe0$D17Oi*1000))+
  
  # Late Holocene
  geom_segment(aes(x=-Inf, xend=Inf, y=-8, yend=-8), color="#44cf6c", lwd=1)+
  geom_segment(aes(x=-Inf, xend=Inf, y=-2, yend=-2), color="dodgerblue3", lwd=1)+
  geom_segment(aes(x=-Inf, xend=Inf, y=-24, yend=-24), color="red", lwd=1)+
  
  # Early and Mid Holocene
  # geom_segment(aes(x=-Inf, xend=Inf, y=-18, yend=-18), color="#44cf6c", lwd=1, lty=1)+
  # geom_segment(aes(x=-Inf, xend=Inf, y=-15, yend=-15), color="dodgerblue3", lwd=1, lty=1)+
  # geom_segment(aes(x=-Inf, xend=Inf, y=-42, yend=-42), color="red", lwd=1, lty=1)+
  # geom_rect(aes(xmin=0.5, xmax=0.9, ymax=-43+6, ymin=-43-6), fill="red", alpha=.2)+
  
  labs(x="Relative humidity", y=expression(Delta*"\u02B9"^"17"*"O"[rlw]*" - "*Delta*"\u02B9"^"17"*"O"[I]*" (per meg)"), size=5)+ 
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14),
        panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA))+
  geom_text(aes(x=0.05, y=subset(dat, dat[,2] == .1)$Dp17Ow-D17Oi*1000), alpha=1, label=c(seq(0, 1, by=0.1)),  hjust = 0, color="orange")+
  scale_y_continuous(limits = c(-120, 2), expand = c(0, 0), labels = scales::number_format(accuracy = 1))+
  scale_x_continuous(limits = c(0.05,0.9), n.breaks=8, labels = scales::number_format(accuracy = .1))

Fig6
# ggsave(filename="Fig6.pdf", plot = Fig6, path=plot.path, device=cairo_pdf, height=6, width=5 )






