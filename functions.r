#Important equations:
#  RQ = CO2 eliminated/O2 consumed
#  EE = 0.06 * (3.941 * VO2 + 1.106 * VCO2)
#  
#  from 14.4 in Leighton book
#  LabDiet 5015 = (26.101/100).71+(19.752/100).83+.54148 = .8907387
#  LabDiet low fat 5015 = (22.8/100).71+(6.6/100).83+.706 = .92266


#import data function
bring_in_data <- function(data_file)
{
  data <- paste(data_file,sep="")
  raw <- read_csv(data, skip_empty_rows=TRUE,
                  col_types = cols(Animal = col_integer(), 
                                   deltaCO2 = col_double(), 
                                   deltaH2O = col_double(),
                                   H2Oml = col_double(),
                                   Deg_C = col_double(),
                                   VCO2 = col_double()))
  
  '%!in%' <- function(x)!('%in%'(x))
  
    raw <- raw %>% 
    mutate(EE = 0.06*(3.941*VO2 + 1.106*VCO2)) %>% 
    mutate(RQ = VCO2/VO2) %>%
    rename(c("StartTime" = "time", "StartDate" = "date")) %>%
    #mutate(time = format(time, "%H:%M:%S", tz="EST")) %>%
    mutate(date, date = as.POSIXlt(date, format = "%d-%b-%y", tz="EST")) %>%
    unite("DateTime", date:time, remove = FALSE, sep =  " ") %>%
    #mutate(DateTime = as.POSIXlt(DateTime), tz="EST")  %>%
    #mutate(DateTime, DateTime = as.POSIXlt(ymd_hms(DateTime), tz="EST"))
    mutate(weight = 
             ifelse(Animal == 0, mean(mouse_metadata_rep1$weight[mouse_metadata_rep1$cage == 0], na.rm = TRUE), 
             ifelse(Animal == 1, mean(mouse_metadata_rep1$weight[mouse_metadata_rep1$cage == 1], na.rm = TRUE),
             ifelse(Animal == 2, mean(mouse_metadata_rep1$weight[mouse_metadata_rep1$cage == 2], na.rm = TRUE),
             ifelse(Animal == 3, mean(mouse_metadata_rep1$weight[mouse_metadata_rep1$cage == 3], na.rm = TRUE),
             ifelse(Animal == 4, mean(mouse_metadata_rep1$weight[mouse_metadata_rep1$cage == 4], na.rm = TRUE),
             ifelse(Animal == 5, mean(mouse_metadata_rep1$weight[mouse_metadata_rep1$cage == 5], na.rm = TRUE),
             ifelse(Animal == 6, mean(mouse_metadata_rep1$weight[mouse_metadata_rep1$cage == 6], na.rm = TRUE), NA)))))))) %>% 
      mutate(sex = 
               ifelse(Animal == 0, na.omit(mouse_metadata_rep1$sex[mouse_metadata_rep1$cage == 0])[1], 
               ifelse(Animal == 1, na.omit(mouse_metadata_rep1$sex[mouse_metadata_rep1$cage == 1])[1],
               ifelse(Animal == 2, na.omit(mouse_metadata_rep1$sex[mouse_metadata_rep1$cage == 2])[1],
               ifelse(Animal == 3, na.omit(mouse_metadata_rep1$sex[mouse_metadata_rep1$cage == 3])[1],
               ifelse(Animal == 4, na.omit(mouse_metadata_rep1$sex[mouse_metadata_rep1$cage == 4])[1],
               ifelse(Animal == 5, na.omit(mouse_metadata_rep1$sex[mouse_metadata_rep1$cage == 5])[1],
               ifelse(Animal == 6, na.omit(mouse_metadata_rep1$sex[mouse_metadata_rep1$cage == 6])[1], NA)))))))) %>% 
    mutate(Animal_ID = 
             ifelse(Animal == 0, mouse_metadata_rep1$animal_id[mouse_metadata_rep1$cage == 0][1], 
             ifelse(Animal == 1, mouse_metadata_rep1$animal_id[mouse_metadata_rep1$cage == 1][1],
             ifelse(Animal == 2, mouse_metadata_rep1$animal_id[mouse_metadata_rep1$cage == 2][1],
             ifelse(Animal == 3, mouse_metadata_rep1$animal_id[mouse_metadata_rep1$cage == 3][1],
             ifelse(Animal == 4, mouse_metadata_rep1$animal_id[mouse_metadata_rep1$cage == 4][1],
             ifelse(Animal == 5, mouse_metadata_rep1$animal_id[mouse_metadata_rep1$cage == 5][1],
             ifelse(Animal == 6, mouse_metadata_rep1$animal_id[mouse_metadata_rep1$cage == 6][1], NA)))))))) 
    #mutate(H2Omg_edit = 
             #ifelse(hour(StartTime) == 8, H2Omg,
             #ifelse(hour(StartTime) == 7, H2Omg,
             #ifelse(hour(StartTime) == 9, H2Omg,
             #ifelse(hour(StartTime) == 10, H2Omg,
             #ifelse(hour(StartTime) == 19, H2Omg,
             #ifelse(hour(StartTime) == 20, H2Omg,
             #ifelse(hour(StartTime) == 21, H2Omg,
             #ifelse(hour(StartTime) == 22, H2Omg,
             #ifelse(hour(StartTime) %!in% c(7,8,9,10,20,21,22,19), H2Omg, NA)))))))))) %>%
    #mutate_at("H2Omg_edit", as.numeric)
  
  target <- c(0,1,2,3,4,5,6,7)
  cages <- raw %>% filter(Animal %in% target)
  
  
  #start_time <- ymd_hms(subset[[5]][1])
  #begin_experiment <- start_time + dhours(2)
  #end_time <- begin_experiment + dhours(72)
  #filtered <- subset %>% filter(raw$DateTime >= begin_experiment & raw$DateTime <= end_time)
  
  
  return(cages)
}




save <- function(name, fig_len, fig_width,plot_name)
{
  library(Cairo)
  
  Cairo::Cairo(
    fig_len, #length
    fig_width, #width
    file = paste(name, ".tiff", sep = ""),
    type = "tiff", #tiff
    bg = "transparent", #white or transparent depending on your requirement 
    dpi = 600,
    units = "cm" #you can change to pixels etc 
  )
  plot(plot_name) #p is your graph object 
  dev.off()
}

## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}