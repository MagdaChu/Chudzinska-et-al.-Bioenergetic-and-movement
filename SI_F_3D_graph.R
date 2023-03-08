rm(list = ls())
require(plotly)
require(patchwork)
library(plot3D)


## load data
data <- read.csv("C:/Users/mec21/OneDrive - University of St Andrews/SMRU/DEB_SMRUC/DEB paper code/All results/SI_F.csv")
data$Dist_effect <- as.numeric(as.character(unique(data$Dist_effect)))

x <- unique(data$p_dist)
y <- unique(data$Dist_effect)

ax <- list(
  nticks = 6,
  range = c(0,1),
  title = "Exposure"
)

ay <- list(
  nticks = 6,
  range = c(1,12),
  title = "Dist. effect"
)

# Harbour porpoise
#calf birth
hpz <- matrix(data$HP_BirthRate, ncol=length(x), nrow=length(y))
hpcb <- list(x,y,hpz)
names(hpcb) <- c("x","y","z")

az1 <- list(
  nticks = 6,
  range = c(-20,0),
  title = "Birth rate"
)

# calf surv
  hpz <- matrix(data$HP_CalfSurv, ncol=length(x), nrow=length(y))
  hpcs <- list(x,y,hpz)
  names(hpcs) <- c("x","y","z")
  
  az2 <- list(
    nticks = 6,
    range = c(0,50),
    title = "Calf sr"
  )
  
# adult surv
hpz <- matrix(data$HP_AdSurv, ncol=length(x), nrow=length(y))
hpas <- list(x,y,hpz)
names(hpas) <- c("x","y","z")

az3 <- list(
  nticks = 6,
  range = c(0,70),
  title = "Adult sr"
)  

# subplot and define scene

ghpcb <- plot_ly(x=hpcb$x, y=hpcb$y,  z = hpcb$z,scene='scene1') %>% 
  add_surface(showscale=FALSE)  

ghpcs <- plot_ly(x=hpcs$x, y=hpcs$y,  z = hpcs$z,scene='scene2') %>% 
  add_surface(showscale=FALSE) 

ghpas <- plot_ly(x=hpas$x, y=hpas$y,  z = hpas$z,scene='scene3') %>% 
  add_surface(showscale=FALSE) 

figHP <- subplot(ghpcb, ghpcs, ghpas) 

figHP %>% layout(title = "Harbour porpoise",
               scene = list(domain=list(x=c(0,0.5),y=c(0.5,1)),
                            xaxis=ax, yaxis=ay, zaxis=az1,
                            aspectmode='cube'),
               scene2 = list(domain=list(x=c(0.5,1),y=c(0.5,1)),
                             xaxis=ax, yaxis=ay, zaxis=az2,
                             aspectmode='cube'),
               scene3 = list(domain=list(x=c(0,0.5),y=c(0,0.5)),
                             xaxis=ax, yaxis=ay, zaxis=az3,
                             aspectmode='cube'))


# grey seals

#pup birth
gsz <- matrix(data$GS_BirthRate, ncol=length(x), nrow=length(y))
gscb <- list(x,y,gsz)
names(gscb) <- c("x","y","z")

az4 <- list(
  nticks = 6,
  range = c(-30,0),
  title = "Birth rate"
)

# pup surv
gsz <- matrix(data$GS_CalfSurv, ncol=length(x), nrow=length(y))
gscs <- list(x,y,gsz)
names(gscs) <- c("x","y","z")

az5 <- list(
  nticks = 6,
  range = c(0,15),
  title = "Pup sr"
)

# adult surv
gsz <- matrix(data$GS_AdSurv, ncol=length(x), nrow=length(y))
gsas <- list(x,y,gsz)
names(gsas) <- c("x","y","z")

az6 <- list(
  nticks = 6,
  range = c(0,50),
  title = "Adult sr"
)


# subplot and define scene

ggscb <- plot_ly(x=gscb$x, y=gscb$y,  z = gscb$z ,scene='scene1') %>% 
  add_surface(showscale=FALSE)  

ggscs <- plot_ly(x=gscs$x, y=gscs$y,  z = gscs$z,scene='scene2') %>% 
  add_surface(showscale=FALSE) 

ggsas <- plot_ly(x=gsas$x, y=gsas$y,  z = gsas$z,scene='scene3') %>% 
  add_surface(showscale=FALSE) 

figGS <- subplot(ggscb, ggscs, ggsas) 

figGS %>% layout(title = "Grey seal",
                 scene = list(domain=list(x=c(0,0.5),y=c(0.5,1)),
                              xaxis=ax, yaxis=ay, zaxis=az4,
                              aspectmode='cube'),
                 scene2 = list(domain=list(x=c(0.5,1),y=c(0.5,1)),
                               xaxis=ax, yaxis=ay, zaxis=az5,
                               aspectmode='cube'),
                 scene3 = list(domain=list(x=c(0,0.5),y=c(0,0.5)),
                               xaxis=ax, yaxis=ay, zaxis=az6,
                               aspectmode='cube'))






