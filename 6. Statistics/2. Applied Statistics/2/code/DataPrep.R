#Files
# setwd("")

# We load the files
Camlob1 <-read.delim("campy_pre2002.txt")
Camlob2 <-read.csv("campy_2002-2005.csv", sep=",")
Camlob3 <-read.csv("campy_2005-.csv", sep=",")

# 1) pre2002: Remove those with SEKTION=="res"
Camlob1<-Camlob1[!(Camlob1$SEKTION == "res"),]
# 2) pre2002: Only keep those with AKTVNR==5133
Camlob1<-Camlob1[(Camlob1$AKTVNR == 5133),]
# 3) All files: Valid CHR numbers are 10000 and above
Camlob1<-Camlob1[(Camlob1$CHR_NR >= 10000),]
Camlob2<-Camlob2[(Camlob2$Chrnr >= 10000),]
Camlob3<-Camlob3[(Camlob3$Chrnr >= 10000),]

#4) Convert dates to common format.
# Camlob1$PRV_DATO <- as.Date(Camlob1$PRV_DATO)

# The first format transformation only works for some months !?!?!?
Camlob1$PRV_DATO <- as.Date(strptime(Camlob1$PRV_DATO, format = "%d%b%Y:%H:%M:%S"))
Camlob2$Prvdato <- as.Date(Camlob2$Prvdato, format = "%m/%d/%Y")
Camlob3$Provedato <- as.Date(Camlob3$Provedato, format = "%m/%d/%Y")

# as.Date(Camlob1$PRV_DATO[4])
# as.Date(Camlob3$Provedato[4])
# as.Date(Camlob2$Prvdato[3])
# typeof(Camlob1$PRV_DATO[3])
# typeof(Camlob2$Prvdato[3])
# typeof(Camlob3$Provedato[4])

# 5)get same order of columns to keep and then rename.
## We will change numbers to lower case in 1 and put the names of 3.
colnames(Camlob1)[which(names(Camlob1) == "JNR")] <- "Jnr"
colnames(Camlob1)[which(names(Camlob1) == "DYRNR")] <- "Dyrnr"
colnames(Camlob1)[which(names(Camlob1) == "MATR")] <- "Materialeart"
colnames(Camlob1)[which(names(Camlob1) == "EPINR")] <- "Epinr"
colnames(Camlob1)[which(names(Camlob1) == "PRV_DATO")] <- "Provedato"
colnames(Camlob1)[which(names(Camlob1) == "MATR")] <- "Materialeart"
colnames(Camlob1)[which(names(Camlob1) == "CHR_NR")] <- "Chrnr"
## "BAKTFUND" Blank if no campy and subspecies if positive
colnames(Camlob1)[which(names(Camlob1) == "BAKTFUND")] <- "Resultat"

colnames(Camlob2)[which(names(Camlob2) == "Epi.nr")] <- "Epinr"
colnames(Camlob2)[which(names(Camlob2) == "Prvdato")] <- "Provedato"

colnames(Camlob3)[which(names(Camlob3) == "Tolkning")] <- "Resultat"
## So the common variables are:
common_columns = c("Jnr","Chrnr", "Epinr", "Provedato", "Materialeart", "Resultat")

# 6) Some tests are recorded in two files with different JNR!?! 
# (Due to transitions between databases ...) 
# (This step can be skipped at first and handled if time permits)
# 7) Merge the data using "rbind"
# rbind will only join if both dataframes have the same columns (not ordered necesarily)
# merge would do it better though
Camlob1 <- Camlob1[common_columns]
Camlob2 <- Camlob2[common_columns]
Camlob3 <- Camlob3[common_columns]

## Join the data sets
Camlob <- rbind(Camlob1,Camlob2)
Camlob <- rbind(Camlob,Camlob3)
head(Camlob)

# 8) Remove records with chrnr<=10000 and those with NA as epinr.
# Hint: Use "!is.na(epinr)"
Camlob <- Camlob[!(Camlob$Chrnr <= 10000),]
Camlob <- Camlob[!is.na(Camlob$Epinr),]

# 9) Reduce the levels of resultat to only "POS" or "NEG"
# We have to transform the POSITIV to POS, and NEGATIV to NEG
Camlob[(Camlob$Resultat == "NEGATIV"),]$Resultat <- "NEG"
Camlob[(Camlob$Resultat == "POSITIV"),]$Resultat <- "POS"
## We fill the empty with NEG and the rest with POC
Camlob[(Camlob$Resultat == ""),]$Resultat <- "NEG"
Camlob[(Camlob$Resultat != "NEG" & Camlob$Resultat != "POS"),]$Resultat <- "POS"
Camlob$Resultat <- factor(Camlob$Resultat)  # This way we remove the others

# 10) Remove records with duplicated jnr (Keep first record)
# Hint: Use "duplicated"
Camlob <-Camlob[!(duplicated(Camlob$Jnr)),]

# 11) Only keep records with "matr" in c("Kloaksvaber","Svaberprøve","766","772")
Camlob = Camlob[(Camlob$Materialeart == "766" | (Camlob$Materialeart == "772") |
                   (Camlob$Materialeart == "Kloaksvaber")|(Camlob$Materialeart == "Svaberprøve")),]

# 12) Add week number since week one 1998 for each record
# Camlob$Provedato
weeks = as.integer(strftime(Camlob$Provedato,format="%W"))
years = as.integer(strftime(Camlob$Provedato,format="%Y"))
Camlob$weeks <- weeks
Camlob$years <- years

unique(years)
unique(weeks)
# 13) Only keep those with positive week number
Camlob <-Camlob[Camlob$weeks >= 0,]
# 14) Some chrnr should be removed due to various reasons.
# Skip!

# 15) It may be decided only to include data from farms that 
# have delivered more than 10 flocks, as those with less may have 
# a bias. This could also be included in the analysis ...

# 16) Use "split" to split the data by week
splitted_Camlob = split(Camlob, list(Camlob$years, Camlob$weeks ))
splitted_Camlob = split(Camlob, list(years))
splitted_Camlob = split(splitted_Camlob, list(weeks))
sapply(Camlob, function(x) apply(x, Resultat, mean))

get_ratio <- function(x){
  Total = length(x)
  Npos = sum(x == "POS")
  ratio = as.double(Npos)/Total
}

jk = by(Camlob$Resultat, list(years,weeks), get_ratio)
plot(jk)

splitted_Camlob = split(Camlob, weeks)
# 17) Summarize number of flocks slaughtered and number of positive flocks per week.
summary(splitted_Camlob$'32'$Resultat == "POS")

## Apply summary to all weeks.
tapply(Camlob$Resultat, Camlob$weeks, summary)
# 18) Save your data file!
write.csv(Camlob, file = "Camlob.csv")

# Get the ratio and do some plotting !!
split_names = names(splitted_Camlob)

ratio = c()
for (sn in split_names){
  Total = length(splitted_Camlob[[sn]]$Resultat)
  print(Total)
  Npos = sum(splitted_Camlob[[sn]]$Resultat == "POS")
  print(Npos)
  print(as.double(Npos)/Total)
  ratio = c(ratio, as.double(Npos)/Total)
}

plot(Camlob$Resultat)
plot(ratio)

