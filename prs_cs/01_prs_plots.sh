R
library(data.table)
library(plotrix)
library(scales)

#Read in a results file consisting of log ORs (or betas) for quintiles of PRS risk
#Input file should have project name, the quintile, and corresponding beta and standard errors
res2 <- fread('tinnitus_OR.csv',data.table=F)
names(res2) <- c("Project","Decile","Estimate","Std. Error")

#Here we calculate plot for 3 groups. Hence we perturb the values slightly for the x axis, so the plotted dots don't overlap
res2[which(res2$Project=="UKBMVPtoR4"),]$Decile <- 1:5 + 0.075
res2[which(res2$Project=="UKBtoR3"),]$Decile <- 1:5
res2[which(res2$Project=="MVPtoUKB"),]$Decile <- 1:5 - 0.075

#Calculate OR (if applicable)
res2$OR <- exp(res2$Estimate)
res2$LCI <- exp(res2$Estimate - 1.96*res2[,"Std. Error"])
res2$UCI <- exp(res2$Estimate + 1.96*res2[,"Std. Error"])

#Pick colors for results (here 3 groups)
res2$color <- NA
caro_orange <- rgb(249,164,24,maxColorValue=255)
beauty_red <- rgb(237,33,36,maxColorValue=255)
adam_blue <- rgb(41,153,208,maxColorValue=255)

res2[which(res2$Project=="UKBMVPtoR4"),]$color <-adam_blue
res2[which(res2$Project=="UKBtoR3"),]$color <-caro_orange
res2[which(res2$Project=="MVPtoUKB"),]$color <- beauty_red


cis2 <- rbind(res2)


pdf('tinnitus_allgroups.pdf',5.5,7)
par(mar=c(5, 4, 4, 2) + 0.5)
plotCI(x=cis2$Decile,y=cis2$OR,li=cis2$LCI,ui=cis2$UCI,lwd=2,ylim=c(1,1.75),pch=19,cex.axis=1.25,xlab="PRS Quintile",ylab="Quintile Odds Ratio (95% CI)",main="",cex.lab=1.45,col=cis2$color,scol=alpha(cis2$color,.4),sfrac=0,xaxt='n',yaxt='n',cex=1.8)

legend("topleft",col=c("white",beauty_red,caro_orange,adam_blue),legend=c("Training -> Target", "MVP -> UKB", "UKB -> MVP", "Meta-analysis -> MVP R4"),bty="n",pch=19,cex=1.5)

axis(1,at=c(1:5),cex.axis=1.25)
axis(2,at=c(1,1.25,1.5,1.75), labels=c("1","1.25","1.5","1.75"),cex.axis=1.25)


dev.off()

