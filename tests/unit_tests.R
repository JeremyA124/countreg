library(countmods)
library(tidyverse)

# Data set preparation
arma.dat <- read.csv("McMillanAcheArmadillo.csv")
cam.dat <- read.csv("cs_replication_data.csv")
nulldats <- arma.dat
nulldats$Age[sample(1:length(nulldats$Age), 12)] <- NA
nulldats$Armadillos[sample(1:length(nulldats$Armadillos), 12)] <- NA

nandats <- arma.dat
nandats$Age[sample(1:length(nandats$Age), 12)] <- NaN
nandats$Armadillos[sample(1:length(nandats$Armadillos), 12)] <- NaN

# Data Plotting tests
dat_plot(arma.dat)
dat_plot(cam.dat) #Should not allow
dat_plot(nulldats)
dat_plot(nandats)

# Model fit test
a <- glm_pois(data = arma.dat, Armadillos~Age)
b <- glm_pois(data = arma.dat, Armadillos~Age, offsetparm = "Treks")
c <- glm_pois(data = cam.dat, cam_count~pnhwht+pnhblk)
d <- glm_pois(data=nulldats, Armadillos~Age)
e <- glm_pois(data=nandats, Armadillos~Age)

f <- glm_negb(data = arma.dat, Armadillos~Age)
g <- glm_negb(data = arma.dat, Armadillos~Age, offsetparm = "Treks")
h <- glm_negb(data = cam.dat, cam_count~pnhwht+pnhblk)
i <- glm_negb(data=nulldats, Armadillos~Age)
j <- glm_negb(data=nandats, Armadillos~Age)

k <- glm_pois_GP2(data = arma.dat, Armadillos~Age)
l <- glm_pois_GP2(data = arma.dat, Armadillos~Age, offsetparm = "Treks")
m <- glm_pois_GP2(data = cam.dat, cam_count~pnhwht+pnhblk)
n <- glm_pois_GP2(data=nulldats, Armadillos~Age)
o <- glm_pois_GP2(data=nandats, Armadillos~Age)

p <- glm_pois_zero(data = arma.dat, Armadillos~Age+I(Age^2), Armadillos~1)
q <- glm_pois_zero(data = arma.dat, Armadillos~Age+I(Age^2), Armadillos~Age, offsetparm = "Treks")
r <- glm_pois_zero(data = cam.dat, cam_count~pnhwht+pnhblk, cam_count~pnhwht)
s <- glm_pois_zero(data=nulldats, Armadillos~Age+I(Age^2), Armadillos~1) #There are issues here, but unfortunately
t <- glm_pois_zero(data=nandats, Armadillos~Age+I(Age^2), Armadillos~1) #I have to deploy the software soon.

u <- glm_negb_zero(data = arma.dat, Armadillos~Age, Armadillos~1)
v <- glm_negb_zero(data = arma.dat, Armadillos~Age, Armadillos~Age, offsetparm = "Treks")
w <- glm_negb_zero(data = cam.dat, cam_count~pnhwht+pnhblk, cam_count~pnhwht)
x <- glm_negb_zero(data=nulldats, Armadillos~Age+I(Age^2), Armadillos~1) # Same issues with glm_pois_zero using
y <- glm_negb_zero(data=nandats, Armadillos~Age+I(Age^2), Armadillos~1) # glm_negb_zero

# Residual Diagnostics tests
resid_diag(a) #Should all work
resid_diag(n)
resid_diag(u)
resid_diag(r)
resid_diag(g)

# Interpret test
interpret(a) #Should all work
interpret(n)
interpret(u)
interpret(r)
interpret(g)

# Zero-Inf tests
zeroinf_test(a) #Should all work
zeroinf_test(b)
zeroinf_test(c)
zeroinf_test(d)
zeroinf_test(f)
zeroinf_test(g)
zeroinf_test(h)
zeroinf_test(i)

# Model Comaprison test
mod_compare(a,a) #same models (Shoul all work)
mod_compare(d,d)
mod_compare(p,q) #nested models
mod_compare(u,v)
