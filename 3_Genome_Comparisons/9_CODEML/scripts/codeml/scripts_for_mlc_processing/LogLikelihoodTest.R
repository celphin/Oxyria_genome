#install.packages("extRemes")
library(extRemes)
lnL <- read.table("BOTH.mlc",header=F)
new_lnL <- as.data.frame.matrix(lnL) 
Output_lnL <- cbind(new_lnL)

#This works:
#NB: The test follows D = -2*(y - x), where y stands for m1 and x m0, but the command is lr.test(x, y, alpha = 0.05, df = 1, ...).  Normally, M0 values should come first in the expression. 
#HOWEVER I found that due to negative log likelihood values, the expression needs to be turned around so that m1 comes first and m2 comes second.
for(i in 1:nrow(new_lnL)) {
filename1 <- print(new_lnL[i,1])
filename2 <- paste(filename1,".txt", sep="")
sink(filename2, append=TRUE)    
print(new_lnL[i,1]);print(lr.test(new_lnL[i,2],new_lnL[i,3], alpha = 0.05, df = 1));
sink()
}

#This is how it normally works:
#lr.test(-1245.757897, -1245.927791, alpha = 0.05, df = 1)
