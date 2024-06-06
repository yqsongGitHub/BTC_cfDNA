## Cell-Free DNA in Plasma Reveals Genomic Similarity Between Biliary Tract Inflammatory Lesion and Biliary Tract Cancer

### Data analysis scripts overview:
All the main and supplementary figures and statistics for the manuscript are generated in these scripts.

Script                       | Brief description                                       
-----------------------------| ------------------------------------
library.R                    |  The required R packages              
function.R                   | The required functions to be loaded  
cal_AFV.R                    | Example of Calculating AFV        
fig2.R                       | Fig. 2a, Fig. 2b, Fig. 2c, Fig. 2d, Fig. 2e          
fig3.R                       | Fig. 3a, Fig. 3b, Fig. 3c, Fig. 3d, Fig. 3e           
fig4.R                       | Fig. 4b, Fig. 4c, Fig. 4d, Fig. 4e, Fig. 4f, Fig. 4g 
fig5.R                       | Fig. 5a, Fig. 5b, Fig. 5c 
figSup.R                     | Fig. S3a, Fig. S3b, Fig. S4a, Fig. S4b, Fig. S5a, Fig. S5b 

### Running AFV overview:
We employed the allele fraction variance (AFV) metric as a measure of similarity between cfDNA samples and WBC samples. We aimed to identify somatic variants in cfDNA tumors and blood tissue variants using a minimum depth cut-off of 200. We treated allele frequency values of genomic sites as independent and identically distributed random variables with a common distribution function that indicates the fraction of tumor-derived cfDNAs. To correct for sequencing errors and germline variants, we used the WBC sample as a control when calculating the variant allele frequency (VAF). For each patient, we calculated AFV using variant allele frequencies. We generated a scatter plot for meaningful genomic sites, with the VAF of a tumor sample on the Y axis and the VAF of a paired normal or blood sample on the X axis. Next, we created a diagonal line on which the points have the same VAF between both samples.

<p align="center">
<img align="center" src="https://github.com/BTC_AFV/images/Schematic.png">
</p>   

Example data for testing is available at [BTC_AFV/example](https://github.com/BTC_AFV/images/).

#### One sample
```
dat <- read.csv("./example/Cholangiocarcinoma/CJY_cfDNA1.csv")
AFV(dat[,3],dat[,4])

```

#### Multiple samples
```
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

```
--------

## Questions & Answers
Please submit a github issue with any questions or if you experience any issues/bugs.
- Create a new [issue](https://github.com/BTC_AFV/issues).

--------
## Warning   
No tested bugs needed warning. 

--------
### Citation
Cite as: Liu, R., Song, Y., Hua, R. et al. Cell-Free DNA in Plasma Reveals Genomic Similarity Between Biliary Tract Inflammatory Lesion and Biliary Tract Cancer. Phenomics (2024). https://doi.org/10.1007/s43657-024-00160-2
