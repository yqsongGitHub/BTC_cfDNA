#### calculate AFV ####
dist2d <- function(a,b,c) {
  v1 <- b - c
  v2 <- a - b
  m <- cbind(v1,v2)
  d <- abs(det(m))/sqrt(sum(v1*v1))
  return(d)
} 
AFV <- function(x,y){
  dd0 <- c()
  NN <- length(x)
  for(i in 1:NN){
    pp0 <- c(as.numeric(x[i]),as.numeric(y[i]))
    dd1 <- dist2d(pp0,c(0,0),c(1,1)) 
    dd0 <- c(dd0,dd1)
  }
  return(median(dd0,na.rm = T))
}

# one sample 
dat <- read.csv("./example/Cholangiocarcinoma/CJY_cfDNA1.csv")
AFV(dat[,3],dat[,4])

# loop
wd0 <- "./example/Cholangiocarcinoma/"
suffix <- "cfDNA1"
files <- list.files(wd0)
files <- files[grepl(suffix,files)]
samples <- gsub("_cfDNA.*","",files)
df0 <- data.frame()
for(i in samples){
  print(i)
  dat <- read.csv(paste0(wd0,i,"_",suffix,".csv"))
  AFV_value <- AFV(dat[,3],dat[,4])
  df <- data.frame("sample"=paste0(i,"_",suffix),
                   "AFV"=AFV_value)
  df0 <- rbind(df0,df)
}
AFV_re <- df0
