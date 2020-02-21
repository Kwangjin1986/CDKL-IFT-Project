# ciliary length plotting
# 2017-08-21 Kwangjin Park - Sixth Edition
# 2016-02-05 Kwangjin Park - Fifth Edition
# 2016-01-20 Kwangjin Park - Fourth Edition
# 2015-11-27 Kwangjin Park - Third Edition
# 2015-11-23 Kwangjin Park - Second Edition
# 2015-10-18 Kwangjin Park - First Edition
# rscript "C:\Users\JakePark\Documents\protein_expression_and_length\bin\ciliary_length.R" nx131.csv nx131figure.pdf pvalue-nx131.csv


if(!require("ggplot2")) {
    install.packages("ggplot2", dependencies = TRUE)
  }
if(!require("plyr")) {
    install.packages("plyr", dependencies = TRUE)
  }
  if(!require("broom")) {
    install.packages("broom", dependencies = TRUE)
  }
  if(!require("FSA")) {
    install.packages("FSA", dependencies = TRUE)
  }

library(ggplot2)
library(plyr)
library(broom)
library(FSA)

args <- commandArgs(trailingOnly = TRUE)

datafile <- args[1]
savefile <- args[2]
pvaluefile <- args[3]

#load data
cell_body <- read.csv(datafile, 
                       header=TRUE, as.is=TRUE, strip.white=TRUE)

					   
# change types of values of STRAIN column (Charachter to factors)
cell_body$STRAIN = as.factor(cell_body$STRAIN)
original_order <- unique(cell_body$STRAIN)
factor = levels(cell_body$STRAIN)
cell_body$STRAIN = factor(cell_body$STRAIN, levels(cell_body$STRAIN)[match(original_order, factor)])
					   
#plot(box and dot plots)
plotbasic <- ggplot(data=cell_body, aes(x= STRAIN, y= CILIARY_LENGTH))+
  theme(plot.title = element_text(size=12, face="bold", vjust=2), ## make the plot title larger and higher
        axis.text.x=element_text(colour="black", size = 12), ## change the x-axis values font to black
        axis.text.y=element_text(colour="black", size = 12), ## change the y-axis values font to black and make larger
        axis.title.y = element_text(size = 16, vjust = 1.3)) +
  ggtitle("ciliary length of ADL in wildtype and tm4182")+
  xlab("")+ylab("Length of ADL cilia(um)")+
  geom_jitter(alpha = 0.7, position = position_jitter(width = 0.2), size = 2, colour="gray50")+
  geom_boxplot(aes(fill=STRAIN), width = 0.6, outlier.size = 0, alpha=0.3, outlier.colour=NA)


#save figures
ggsave(plotbasic, file=savefile, h=4, w=6, units="in", dpi=300)


#check distribution of data (normal v.s. non-normal) check![Shapiro-wilk test]
Distribution_Ref <- ddply(cell_body, 'STRAIN', summarise,  Shapirotest = shapiro.test(CILIARY_LENGTH)[2])

# If the values of Shapiro-wilk is less than 0.5, the normal distribution is rejected. 
# -> calculate p-value by using "Dunn's Kruskal-Wallis Multiple Comparisons"


if(min(unlist(Distribution_Ref$Shapirotest)) < 0.05){
  
  #If strains are less then 3, P-value will be calculted by "Kruskal-Wallis test"
  if(length(unique(cell_body$STRAIN))<=2){
    P_value1 <- tidy(kruskal.test(CILIARY_LENGTH ~ STRAIN, data = cell_body))
    P_value1$test  <- "Kruskal-Wallis test" 
    write.table(P_value1, file=pvaluefile, row.names = FALSE, quote = FALSE, append = FALSE, sep = ",")
    
  
  }
  
  #If strains are more than 3, p-value will be obstained by "Dunn's Kruskal-Wallis Multiple Comparisons"    
  if(length(unique(cell_body$STRAIN))>=3){
    P_value2 <- dunnTest(CILIARY_LENGTH ~ STRAIN, 
                         data = cell_body, 
                         method="hs")
    
    P_value2$res$test <- rep("Dunn's Kruskal-Wallis test: Holm-Sidak", dim(P_value2$res)[1]) 
    write.table(P_value2$res, file=pvaluefile, row.names = FALSE, quote = FALSE, append = FALSE, sep = ",")}
}

# If the values of Shapiro-wilk is greater than 0.5, the normal distribution is rejected
# ->  calculate p-value by using "oneway anova"

if(min(unlist(Distribution_Ref$Shapirotest)) > 0.05){
  P_value3 <- tidy(TukeyHSD(aov(CILIARY_LENGTH ~ STRAIN,
                                data= cell_body), conf.level = 0.95))
  P_value3$test <- rep("TukeyHSD")
  write.table(P_value3, file=pvaluefile, row.names = FALSE, quote = FALSE, append = FALSE, sep = ",")}
  
