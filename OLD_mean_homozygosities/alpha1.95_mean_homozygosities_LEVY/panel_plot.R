attach(mtcars)





par(mfrow=c(3,3))


dummy <-   read.table("outfile_mu0.1_rho_inverse.01")

 distance <- as.vector(dummy[,1])
 class(distance)
 mean_homozygosity <- as.vector(dummy[,2])
 #class(avg)
 upper <- dummy[,3]
 lower <- dummy[,4]


plot(distance, mean_homozygosity, ylim=range(c(lower, upper)),
    pch=19, main = "Rho_inv = 0.01, Mu = 0.1 ")
arrows(distance, lower, distance, upper, length=0.05, angle=90, code=3)
abline(a = dummy[nrow(dummy),2] +2.95*(nrow(dummy) -1), b = -2.95)

dummy <-   read.table("outfile_mu0.1_rho_inverse.1")

 distance <- as.vector(dummy[,1])
 class(distance)
 mean_homozygosity <- as.vector(dummy[,2])
 #class(avg)
 upper <- dummy[,3]
 lower <- dummy[,4]

plot(distance, mean_homozygosity, ylim=range(c(lower, upper)),
    pch=19, main = "Rho_inv = .1, Mu = 0.1 ")
arrows(distance, lower, distance, upper, length=0.05, angle=90, code=3)
abline(a = dummy[nrow(dummy),2] +2.95*(nrow(dummy) -1), b = -2.95)

dummy <-   read.table("outfile_mu0.1_rho_inverse1")

 distance <- as.vector(dummy[,1])
 class(distance)
 mean_homozygosity <- as.vector(dummy[,2])
 #class(avg)
 upper <- dummy[,3]
 lower <- dummy[,4]

plot(distance, mean_homozygosity, ylim=range(c(lower, upper)),
    pch=19, main = "Rho_inv = 1, Mu = 0.1 ")
arrows(distance, lower, distance, upper, length=0.05, angle=90, code=3)
abline(a = dummy[nrow(dummy),2] +2.95*(nrow(dummy) -1), b = -2.95)

dummy <-   read.table("outfile_mu0.01_rho_inverse.01")

 distance <- as.vector(dummy[,1])
 class(distance)
 mean_homozygosity <- as.vector(dummy[,2])
# class(avg)
 upper <- dummy[,3]
 lower <- dummy[,4]


plot(distance, mean_homozygosity, ylim=range(c(lower, upper)),
    pch=19, main = "Rho_inv = 0.01, Mu = 0.01 ")
arrows(distance, lower, distance, upper, length=0.05, angle=90, code=3)
abline(a = dummy[nrow(dummy),2] +2.95*(nrow(dummy) -1), b = -2.95)

dummy <-   read.table("outfile_mu0.01_rho_inverse.1")

 distance <- as.vector(dummy[,1])
 class(distance)
 mean_homozygosity <- as.vector(dummy[,2])
 #class(avg)
 upper <- dummy[,3]
 lower <- dummy[,4]

plot(distance, mean_homozygosity, ylim=range(c(lower, upper)),
    pch=19, main = "Rho_inv = .1, Mu = 0.01 ")
arrows(distance, lower, distance, upper, length=0.05, angle=90, code=3)
abline(a = dummy[nrow(dummy),2] +2.95*(nrow(dummy) -1), b = -2.95)

dummy <-   read.table("outfile_mu0.01_rho_inverse1")

 distance <- as.vector(dummy[,1])
 class(distance)
 mean_homozygosity <- as.vector(dummy[,2])
 #class(avg)
 upper <- dummy[,3]
 lower <- dummy[,4]

plot(distance, mean_homozygosity, ylim=range(c(lower, upper)),
    pch=19, main = "Rho_inv = 1, Mu = 0.01 ")
arrows(distance, lower, distance, upper, length=0.05, angle=90, code=3)
abline(a = dummy[nrow(dummy),2] +2.95*(nrow(dummy) -1), b = -2.95)
dummy <-   read.table("outfile_mu0.001_rho_inverse.01")

 distance <- as.vector(dummy[,1])
 class(distance)
 mean_homozygosity <- as.vector(dummy[,2])
 #class(avg)
 upper <- dummy[,3]
 lower <- dummy[,4]


plot(distance, mean_homozygosity, ylim=range(c(lower, upper)),
    pch=19, main = "Rho_inv = 0.01, Mu = 0.001 ")
arrows(distance, lower, distance, upper, length=0.05, angle=90, code=3)
abline(a = dummy[nrow(dummy),2] +2.95*(nrow(dummy) -1), b = -2.95)

dummy <-   read.table("outfile_mu0.001_rho_inverse.1")

 distance <- as.vector(dummy[,1])
 class(distance)
 mean_homozygosity <- as.vector(dummy[,2])
 #class(avg)
 upper <- dummy[,3]
 lower <- dummy[,4]

plot(distance, mean_homozygosity, ylim=range(c(lower, upper)),
    pch=19, main = "Rho_inv = .1, Mu = 0.001 ")
arrows(distance, lower, distance, upper, length=0.05, angle=90, code=3)
abline(a = dummy[nrow(dummy),2] +2.95*(nrow(dummy) -1), b = -2.95)

dummy <-   read.table("outfile_mu0.001_rho_inverse1")

 distance <- as.vector(dummy[,1])
 class(distance)
 mean_homozygosity <- as.vector(dummy[,2])
 #class(avg)
 upper <- dummy[,3]
 lower <- dummy[,4]

plot(distance, mean_homozygosity, ylim=range(c(lower, upper)),
    pch=19, main = "Rho_inv = 1, Mu = 0.001 ")
arrows(distance, lower, distance, upper, length=0.05, angle=90, code=3)
abline(a = dummy[nrow(dummy),2] +2.95*(nrow(dummy) -1), b = -2.95)


#plot(wt,mpg, main="Scatterplot of wt vs. mpg")
#plot(wt,disp, main="Scatterplot of wt vs disp")
#hist(wt, main="Histogram of wt")
#boxplot(wt, main="Boxplot of wt")