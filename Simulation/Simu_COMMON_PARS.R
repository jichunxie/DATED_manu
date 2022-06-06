## Parameter setting

R = 100 #number of replicated datasets
seed.vec = seq(R)
L=4
alpha=0.05
N1=1000
thresh_val=5

#multi_thresh <- c(3,7,9,11)
multi_thresh <- seq(3,21,2)

#sampsize.factors for simulation
sampsize.factors = seq(4)
