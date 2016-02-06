library(doParallel)
library(foreach)
library(iterpc)
library(iterators)

# Functions definitions ---------------------------------------------------

GetNumberOfNa <- function(column) {
  return (sum(is.na(column)))
}

GetNumberOfNonNa <- function(column) {
  return (sum(!is.na(column)))
}

GetNumberOfCompleteCases <- function(data, selected.colnames) {
  return (sum(complete.cases(data[selected.colnames])))
}

GetOnlyAddedColumns <- function(vec, colnames.to.include) {
  added.colnames <- setdiff(vec, colnames.to.include)
  # return (added.colnames)
  if (any(grep("vkorc", added.colnames, ignore.case = T))
      & (any(grep("cyp", added.colnames, ignore.case = T)))) {
    return(added.colnames)
  }
}

IterpcIterator <- function(iterpc.object, iteration.length) {
  # one's own function of nextElement() because iterpc
  # returns NULL on finished iteration on subsequent getnext() invocation
  # but not 'StopIteration'
  nextEl <- function() {
    if (iteration.length > 0)
      iteration.length <<- iteration.length - 1
    else
      stop('StopIteration')
    
    getnext(iterpc.object)
  }
  obj <- list(nextElem=nextEl)
  class(obj) <- c('irep', 'abstractiter', 'iter')
  obj
}


# Plot number of NA's -----------------------------------------------------

iwpc <- read.csv("~/yandexDisk/DIPLOMA/data/ipwc_data.csv", na.strings = c("NA", ""))

n.of.na <- apply(iwpc, 2, GetNumberOfNa)
n.of.non.na <- apply(iwpc, 2, GetNumberOfNonNa)
fully.available.columns <- n.of.non.na[n.of.non.na == nrow(iwpc)]
n.of.non.na <- n.of.non.na[n.of.non.na != nrow(iwpc)]

# filter out fully available columns
n.of.non.na <- n.of.non.na[! (names(n.of.non.na) %in% names(fully.available.columns)) ]
forbidden.to.include <- n.of.non.na[n.of.non.na < 1700]
# filter out too les avaiable columns
n.of.non.na <- n.of.non.na[!(names(n.of.non.na) %in% names(forbidden.to.include))]

par(mar=c(18.1,4.1,1.1,1.1), cex=0.8)
plot(n.of.non.na, type='b', xaxt="n", xlab="")
axis(1, at=1:length(n.of.non.na), labels=names(n.of.non.na), las=2)
# abline(1700, 0)
grid()

par(mar=c(7.1,32.1,1.1,1.1), cex=0.7)
plot(n.of.non.na, 1:length(n.of.non.na), type='b', yaxt="n", ylab="", xlab="Number of non NA's in columns")
axis(2, at=1:length(n.of.non.na), labels=names(n.of.non.na), las=1)
abline(1700, 0)
grid()





colnames.to.include <- c("Age",  "Height..cm.", "Weight..kg.",
                         "Therapeutic.Dose.of.Warfarin",
                         "INR.on.Reported.Therapeutic.Dose.of.Warfarin",
                         "Race..Reported.", 
                         "Amiodarone..Cordarone.")

rest.columns.to.choose <- setdiff(names(n.of.non.na), colnames.to.include)


# register parallel backend
registerDoParallel(cores = 3)



results <- vector(mode="list")

for (len in 2:7) {
  combinations <- iterpc(length(rest.columns.to.choose), len)
  it <- IterpcIterator(combinations, getlength(combinations))
  results[[len]] <- foreach(i=it, .combine="cbind") %dopar% {
    colnames <- c(rest.columns.to.choose[i], colnames.to.include)
    total.cases <- sum(complete.cases(iwpc[, colnames]))
    if (total.cases == 1732) {
      colnames
    }
  }
}



grepped.colnames  <- vector(mode = "list")
for (num in 2:6) {
  grepped.colnames[[num]] <- apply(results[[num]], 2, 
                                   function(x) GetOnlyAddedColumns(x, colnames.to.include))
}

grepped.colnames <- sapply(grepped.colnames, function(x) Filter(Negate(is.null), x) )

apply(results[[2]], 2, function(x) GetOnlyAddedColumns(x, colnames.to.include))

experment.colnames <- c(grepped.colnames[[6]]$result.1261179, colnames.to.include, 
                        "Amiodarone..Cordarone.")
sum(complete.cases(iwpc[, experment.colnames]))


[1] "VKORC1.genotype...1639.G.A..3673...chr16.31015190..rs9923231..C.T"   
[2] "VKORC1.QC.genotype...1639.G.A..3673...chr16.31015190..rs9923231..C.T"
[3] "VKORC1.genotype..497T.G..5808...chr16.31013055..rs2884737..A.C"      
[4] "VKORC1.QC.genotype..497T.G..5808...chr16.31013055..rs2884737..A.C"   
[5] "VKORC1.genotype..1173.C.T.6484...chr16.31012379..rs9934438..A.G"     
[6] "VKORC1.QC.genotype..1173.C.T.6484...chr16.31012379..rs9934438..A.G"  
[7] "VKORC1.genotype..1542G.C..6853...chr16.31012010..rs8050894..C.G"     
[8] "VKORC1.QC.genotype..1542G.C..6853...chr16.31012010..rs8050894..C.G"  
[9] "VKORC1.genotype..3730.G.A..9041...chr16.31009822..rs7294...A.G"      
[10] "VKORC1.QC.genotype..3730.G.A..9041...chr16.31009822..rs7294...A.G"   
[11] "VKORC1.genotype..2255C.T..7566...chr16.31011297..rs2359612..A.G"     
[12] "VKORC1.QC.genotype..2255C.T..7566...chr16.31011297..rs2359612..A.G"  
[13] "VKORC1.genotype...4451.C.A..861...Chr16.31018002..rs17880887..A.C"   
[14] "VKORC1.QC.genotype...4451.C.A..861...Chr16.31018002..rs17880887..A.C"
[15] "VKORC1..1639.consensus"                                              
[16] "VKORC1.497.consensus"                                                
[17] "VKORC1.1173.consensus"                                               
[18] "VKORC1.1542.consensus"                                               
[19] "VKORC1.3730.consensus"                                               
[20] "VKORC1.2255.consensus"                                               
[21] "VKORC1..4451.consensus"



GetNumberOfCompleteCases(iwpc, rest.columns.to.choose[grep("vkorc.*conse", rest.columns.to.choose, ignore.case = T)])


save.image("~/yandexDisk/DIPLOMA/data/find_out_na.workspace")




