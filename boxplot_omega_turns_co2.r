# Modified by Kwangjin Park, Jan 9, 2020
# Turn/min code Tiffany Timbers, Oct 15, 2016

#command example) rscript boxplot_omega_turns_co2.r merged.file

# load libraries

require(plyr)
require(stringr)
require(peakPick)
#require(dplyr)
require(tidyverse)



main <- function() {

  args <- commandArgs(trailingOnly = TRUE)

  file <- args[1]

## load MWT merged file 
print("Loading data...")
#parsed.data  <- read.table("merged.file")
parsed.data  <- read.table(file)

## use function to extract column names and change time column from factor to numeric
print("Parsing data...")
parsed.data <- extract.col(parsed.data)

##remove NANs from data
turns_data <- data.frame(parsed.data$strain, parsed.data$plate, parsed.data$time, parsed.data$ID, parsed.data$morphwidth, parsed.data$midline)

#if strain names are all correct ---no need to run this part
turns_data$parsed.data.strain <- as.character(turns_data$parsed.data.strain)
turns_data$parsed.data.strain[turns_data$parsed.data.strain == "N2"] <- "n2"
turns_data$parsed.data.strain[turns_data$parsed.data.strain == "tm4812"] <- "tm4182"
turns_data$parsed.data.strain <- as.factor(turns_data$parsed.data.strain)

##remove NANs from data-CONTINUE
colnames(turns_data) <- c("strain", "plate", "time", "ID", "morphwidth", "midline")
turns_data <- turns_data[complete.cases(turns_data), ]



# detect high amplitude turns (M/m)
turns_data$aspect_ratio <- turns_data$morphwidth / turns_data$midline 


# Detect Peaks of M/m (min_peak_height = 0.25, min_peak_distance = 90)

turns_data2 <- turns_data %>% 
  group_by(plate, ID) %>% 
  nest() %>% 
  mutate(turns = map(data, detect_HA_turns)) %>% 
  filter(!is.na(turns)) %>% 
  select(turns, ID, plate) %>% 
  unnest(cols = c(turns))



# calculate mean N at each time interval (need to chnage 260 to the number that you want)
time_interval <- seq(from = 0, to = 300, by = 20)

# calculate turns per per plate per time bin
turns_by_interval <- turns_data2 %>% 
  mutate(interval = cut(time, time_interval)) 

turns_by_interval2 <- ddply(turns_by_interval,.(strain, plate, interval), summarise,turns=length(ID))
  

# calculate good number per plate per time bin
good_number_by_interval <- turns_data %>% 
  mutate(interval = cut(time, time_interval))

good_number_by_interval2<- ddply(good_number_by_interval,.(strain, plate, interval, ID), summarise, extra=length(ID))
good_number_by_interval3<- ddply(good_number_by_interval2,.(strain, plate, interval), summarise, count=length(ID))


# create dataframe to get proportion of omega-turn for each time bin 
turns_min <- left_join(good_number_by_interval3, turns_by_interval2)
turns_min$turns[is.na(turns_min$turns)] <- 0
turns_min$tpm <- (turns_min$turns / turns_min$count)

turns_min_agg <- turns_min
# create summarised data frame to plot mean turns/min for each time bin
#turns_min_agg <- ddply(turns_min,.(strain, interval), summarise, N =length(plate), 
#                       mean_tpm = mean(tpm), 
#                       sd_tpm = sd(tpm),
#                       se_tpm = sd_tpm / sqrt(N))


turns_min_agg$interval <- as.numeric(str_extract(turns_min_agg$interval, "[0-9]?[0-9]+"))
turns_min_agg$interval <- turns_min_agg$interval + 10 #set to mid of time bin

turns_min_agg_110s <- turns_min_agg[turns_min_agg$interval == 110,]
turns_min_agg_130s <- turns_min_agg[turns_min_agg$interval == 130,]

## save data as a file
print("saving summary...")
write.table(turns_min_agg_110s, file="trimmed_data_110sec.csv", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE, append=FALSE)
write.table(turns_min_agg_130s, file="trimmed_data_130sec.csv", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE, append=FALSE)

#plot(box and dot plots)
plotbasic110 <- ggplot(data=turns_min_agg_110s, aes(x= strain, y= tpm))+
  theme(plot.title = element_text(size=12, face="bold", vjust=2), ## make the plot title larger and higher
        axis.text.x=element_text(colour="black", size = 12), ## change the x-axis values font to black
        axis.text.y=element_text(colour="black", size = 12), ## change the y-axis values font to black and make larger
        axis.title.y = element_text(size = 16, vjust = 1.3)) +
  ggtitle("Omega turns")+
  xlab("")+ylab("Fraction high amplitude turns")+
  geom_jitter(alpha = 0.7, position = position_jitter(width = 0.2), size = 2, colour="gray50")+
  geom_boxplot(aes(fill=strain), width = 0.6, outlier.size = 0, alpha=0.3, outlier.colour=NA)

plotbasic130 <- ggplot(data=turns_min_agg_130s, aes(x= strain, y= tpm))+
  theme(plot.title = element_text(size=12, face="bold", vjust=2), ## make the plot title larger and higher
        axis.text.x=element_text(colour="black", size = 12), ## change the x-axis values font to black
        axis.text.y=element_text(colour="black", size = 12), ## change the y-axis values font to black and make larger
        axis.title.y = element_text(size = 16, vjust = 1.3)) +
  ggtitle("Omega turns")+
  xlab("")+ylab("Fraction high amplitude turns")+
  geom_jitter(alpha = 0.7, position = position_jitter(width = 0.2), size = 2, colour="gray50")+
  geom_boxplot(aes(fill=strain), width = 0.6, outlier.size = 0, alpha=0.3, outlier.colour=NA)



#save figures
print("saving figure...")
ggsave(plotbasic110, file="omega_turn_110sec.pdf", h=4, w=6, units="in", dpi=300)
ggsave(plotbasic130, file="omega_turn_130sec.pdf", h=4, w=6, units="in", dpi=300)
}



#################################################################################################
extract.col <- function(data){
  
  ## split up column V1 into date, plate, time and strain 
  
  date <- str_extract(data$V1, "[0-9]{8}")
  
  plate <- str_extract(data$V1, "[0-9]{8}_[0-9]{6}")
  
  time <- str_extract(data$V1, ":[0-9]+[.][0-9]+")
  
  time <- sub(":", "", time)
  
  strain <- str_extract(data$V1,"[A-Za-z]+[-]?[0-9]+")
  
  
  
  ## combine new columns with merged file
  
  new.data <- cbind(date, plate, strain, time, data[,2:dim(data)[2]])  
  
  
  
  ##rename columns  
  
  colnames(new.data) <- c("date", "plate", "strain", "time", "frame", "ID", "persistance", "area", "speed", "angularspeed", "length", "rellength", "width", "relwidth", "aspect", "relaspect", "midline", "morphwidth", "kink", "bias", "pathlen", "curve", "dir", "loc_x", "loc_y", "vel_x", "vel_y", "orient", "crab", "NA", "NA1", "NA2", "NA3")
  
  
  
  ##replace time column (factor) with time as numeric
  
  new.data$time  <- as.numeric(levels(new.data$time))[new.data$time]
  
  
  
  return(new.data)
  
    
}
#################################################################################################


#################################################################################################
#function "dection of HA turns"  -> move later to the bottom of this script
detect_HA_turns <- function(x, min_peak_height = 0.25, min_peak_distance = 90) {
  
  if (max(x$aspect_ratio) > min_peak_height) { 
    x_ts <- as.matrix(select(x, time, aspect_ratio))
    # neighlim = (min_peak_distance * 25) because neighlim is integer val, and camera
    # records at ~ 25 frames/sec
    peaks_bool <- peakpick(x_ts, neighlim = (min_peak_distance * 25), peak.npos = 50)
  } else {
    return (NA)
  }
  
  turns <- filter(x, peaks_bool[,2], aspect_ratio > min_peak_height)
  return(turns)
}
#################################################################################################
main()