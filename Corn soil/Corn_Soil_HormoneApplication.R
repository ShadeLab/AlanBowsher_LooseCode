# Load Packages, Import Data, Create & Clean Physeq Objects ----------------------------------------------

library("phyloseq")
packageVersion("phyloseq")
library("ggplot2")
packageVersion("ggplot2")
theme_set(theme_bw())
library("reshape2")
packageVersion("reshape2")
library("vegan")
packageVersion("vegan")
sessionInfo()
set.seed(1)

#Import OTU table ('#' in cell A1 has been removed)
OTU_table <- read.table("Clean_OTU_table.txt", sep="\t", row.names=1, header=TRUE)
#Convert to matrix, then to phyloseq format.
OTUs_DNARNA <- data.matrix(OTU_table)
OTUs_DNARNA = otu_table(OTUs_DNARNA, taxa_are_rows = TRUE)

#Import taxonomy (has been formatted in Excel to be tab-delimited)
TAX <- read.table("Clean_otus_tax_assignments.txt", sep="\t", row.names=1, stringsAsFactors = FALSE, header=FALSE, na.strings="")
colnames(TAX) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
#convert to matrix, then phyloseq format.
TAX <- as.matrix(TAX)
TAX = tax_table(TAX)

#Import sample data
METADATA <- read.csv("Metadata.csv", sep=",", row.names=1, header=TRUE, na.strings="")
#convert to phyloseq format.
METADATA<- sample_data(METADATA)

#Create phyloseq object:
physeqDNARNA = phyloseq(OTUs_DNARNA, TAX, METADATA)
physeqDNARNA

#Check for singletons and remove if any.
any(taxa_sums(physeqDNARNA) <= 1) #FALSE

#Remove taxa with no phylum designation (not necessary here, every OTU was ID'd at phylum level)
#physeqDNARNA <- subset_taxa(physeqDNARNA, Phylum!= "")

#Remove taxa with Domain designated as "Eukaryota" (Chloroplast and mitochondria have already been removed on HPCC).
physeqDNARNA <- subset_taxa(physeqDNARNA, Domain!= "Eukaryota") #12695 taxa remaining
physeqDNARNA


# Rarefy and relativize data (Including both DNA and RNA) ----------------------------------------------

#Graph number of reads per OTU, and number of reads per sample)
##readsumsdf = data.frame(nreads = sort(taxa_sums(physeqDNARNA), TRUE), sorted = 1:ntaxa(physeqDNARNA),type = "OTUs")
##readsumsdf = rbind(readsumsdf, data.frame(nreads = sort(sample_sums(physeqDNARNA), TRUE), sorted = 1:nsamples(physeqDNARNA), type = "Samples"))
##title = "Total number of reads"
##readsumsdf
##p = ggplot(readsumsdf, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity")
##p + ggtitle(title) + scale_y_log10() + facet_wrap(~type, 1, scales = "free")

#Make histogram of read counts 
##sample_sum_df <- data.frame(sum = sample_sums(physeqDNARNA))
##ggplot(sample_sum_df, aes(x = sum)) +    
  ##geom_histogram(color = "black", fill = "indianred", binwidth = 2500) + 
  ##ggtitle("Distribution of sample sequencing depth of otus") + 
  ##xlab("Read counts") + 
  ##theme(axis.title.y = element_blank())

#Plot richness vs total reads per sample (if strong positive relationship, tells us that we need to rarefy)
##richness <- estimate_richness(physeqDNARNA, measures = "Observed")
##samplesums <- sample_sums(physeqDNARNA)
##richness_sums <- merge(richness,samplesums, by = "row.names")
##richness_sums
##ggplot(richness_sums, aes(x=y, y = Observed)) + 
  ##geom_point() + 
  ##ggtitle("otu Observed Richness by Total Reads") + 
  ##labs(x = "Total Reads")+
  ##geom_text(label = "") 

#sample_sums(physeqDNARNA)
#sorted <- sort(sample_sums(physeqDNARNA))
#sorted <- as.data.frame(sorted)

# mean, max, min, median of sample read counts
#min(sample_sums(physeqDNARNA)) #22556
#mean(sample_sums(physeqDNARNA)) #100371.2
#max(sample_sums(physeqDNARNA)) #161904
#median(sample_sums(physeqDNARNA)) # 107148.5

# Plot a rarefaction curve: shows #of species detected at given sequencing depth.
#the vertical line is the number of reads to which we are rarefying.
#the horizontal lines indicate the number of species in each sample at raremax.
#DNARNA.mat=t(as(otu_table(physeqDNARNA), "matrix"))
#rarecurve(DNARNA.mat, step = 1000, sample = 22556,label = FALSE)

#Actually do the rarefying here.
set.seed(1)
DNARNA.rfy <- rarefy_even_depth(physeqDNARNA, sample.size = 22556, replace = FALSE)
#1873 OTUs were removed: no longer present in any sample after subsampling.


#make physeq object containing only RNA
OTUs <- otu_table(DNARNA.rfy)
OTUs_RNA.rfy <- OTUs[ , !c(TRUE,FALSE) ] 
RNA.rfy = phyloseq(OTUs_RNA.rfy, TAX, METADATA)
RNA.rfy

#make physeq object containing only DNA
OTUs <- otu_table(DNARNA.rfy)
OTUs_DNA.rfy <- OTUs[ , !c(FALSE,TRUE) ] 
DNA.rfy = phyloseq(OTUs_DNA.rfy, TAX, METADATA)
DNA.rfy

#Convert Rarified data to Relative Abundance Data
DNARNA.rfy.ra = transform_sample_counts(DNARNA.rfy, function(x) x / sum(x) )
#Double-check it worked: sample sums should equal 1 
sample_sums(otu_table(DNARNA.rfy.ra)) 



# Analysis Code for All 8 Treatments (and both DNA & RNA) -----------------

#braycurtis
bc <- phyloseq::distance(DNARNA.rfy.ra, method="bray")

#pcoa plot
set.seed(1)
# Ordinate
DNARNA_pcoa <- ordinate(physeq = DNARNA.rfy.ra, method = "PCoA", distance = "bc")
# Plot
DNARNA_pcoa <- plot_ordination(physeq = DNARNA.rfy.ra,ordination = DNARNA_pcoa,color = "Treatment",shape = "Type",title = "Variation among rhizosphere soil communities (16s)") + geom_point(size = 3)+ scale_colour_brewer(type="qual", palette="Set1") 
#Change size of graph
DNARNA_pcoa + coord_fixed() + scale_color_manual(name="Treatment", breaks=c("Pre-dry", "Post-dry", "Post-water", "Control", "ABA", "IAA", "JA", "SA"), values=c("red","blue","yellow","orange","magenta4","hotpink","chocolate4","green4"))
ggsave("8trts_DNARNA_pcoa.pdf",height=150,units="mm") 


#Permanova test using ADONIS function
# make a data frame from the sample_data
DNARNA.rfy.df <- data.frame(sample_data(DNARNA.rfy.ra))
DNARNA.rfy.df
# ADONIS test
set.seed(1)
adonis(bc ~ Type*Treatment, data = DNARNA.rfy.df)
set.seed(1)
adonis(bc ~ Treatment*Type, data = DNARNA.rfy.df)

#Homogeneity of dispersions test
set.seed(1)
result <- betadisper(bc, DNARNA.rfy.df$Type, type=c("centroid"))
set.seed(1)
permutest(result)

set.seed(1)
result <- betadisper(bc, DNARNA.rfy.df$Treatment, type=c("centroid"))
set.seed(1)
permutest(result)



# Analysis Code for Only 5 Treatments (and both DNA & RNA) ----------------

#Subset to only the four hormones plus control:
DNARNA.rfy.ra.5trts = subset_samples(DNARNA.rfy.ra, Treatment != "Pre-dry" & Treatment != "Post-dry" & Treatment != "Post-water")
sample_data(DNARNA.rfy.ra.5trts)


#braycurtis
bc <- phyloseq::distance(DNARNA.rfy.ra.5trts, method="bray")



#pcoa plot: normalized weighted Unifrac
set.seed(1)
# Ordinate
DNARNA_pcoa <- ordinate(physeq = DNARNA.rfy.ra.5trts, method = "PCoA", distance = "bc")
# Plot
DNARNA_pcoa <- plot_ordination(physeq = DNARNA.rfy.ra.5trts,ordination = DNARNA_pcoa,color = "Treatment",shape = "Type",title = "Variation among rhizosphere soil communities (16s)") + geom_point(size = 3)+ scale_colour_brewer(type="qual", palette="Set1") 
#Change size of graph
DNARNA_pcoa + coord_fixed()+ scale_color_manual(name="Treatment", breaks=c("Control", "ABA", "IAA", "JA", "SA"), values=c("red","blue","yellow","orange","green4"))
ggsave("5trts_DNARNA_pcoa.pdf",width=178,units="mm") 

#Permanova test using ADONIS function
# make a data frame from the sample_data
DNARNA.rfy.df <- data.frame(sample_data(DNARNA.rfy.ra.5trts))
DNARNA.rfy.df
# ADONIS test
set.seed(1)
adonis(bc ~ Type*Treatment, data = DNARNA.rfy.df)

#Homogeneity of dispersions test
set.seed(1)
result <- betadisper(bc, DNARNA.rfy.df$Type, type=c("centroid"))
set.seed(1)
permutest(result)

set.seed(1)
result <- betadisper(bc, DNARNA.rfy.df$Treatment, type=c("centroid"))
set.seed(1)
permutest(result)
TukeyHSD(result, which = "group", ordered = FALSE,conf.level = 0.95)


##Some diversity analyses I explored are below.

#Calculate richness
rich_DNARNA.rfy <- (estimate_richness(DNARNA.rfy, measures = c("Observed")))
Shan_DNARNA.rfy <- (estimate_richness(DNARNA.rfy, measures = c("Shannon")))
div.DNARNA.rfy <- merge(rich_DNARNA.rfy, Shan_DNARNA.rfy, by='row.names')
div.DNARNA.rfy

#calculate pielou's evenness
div.DNARNA.rfy$Pielou=div.DNARNA.rfy$Shannon/log(div.DNARNA.rfy$Observed)
div.DNARNA.rfy
colnames(div.DNARNA.rfy) <- c("ID", "Richness","Shannon Div", "PielouEvenness")
row.names(div.DNARNA.rfy) <- sample_names(DNARNA.rfy)
div.DNARNA.rfy

#Add metadata
metaD <- data.frame(sample_data(DNARNA.rfy))
div.DNARNA.rfy <- merge(div.DNARNA.rfy,metaD, by="row.names")
row.names(div.DNARNA.rfy) <- sample_names(DNARNA.rfy)
div.DNARNA.rfy

#Do ANOVA to compare richness of treatments
options(contrasts=c(unordered='contr.sum', ordered='contr.poly'))
Richness_fit <- aov(Richness ~ Type*Treatment, data=div.DNARNA.rfy, na.action=na.exclude)
drop1(Richness_fit,~.,test="F") # type III SS and F Tests
#TESTING ASSUMPTIONS
Richness_resids <- residuals(Richness_fit)
Richness_preds <- predict(Richness_fit)
plot(Richness_resids ~ Richness_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)

#plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
#test for homogeneity of variance
#Test for normality
shapiro.test(Richness_resids)
qqnorm(Richness_resids)
qqline(Richness_resids)

#Do ANOVA to compare Evenness of treatments
options(contrasts=c(unordered='contr.sum', ordered='contr.poly'))
Evenness_fit <- aov(PielouEvenness ~ Type*Treatment, data=div.DNARNA.rfy, na.action=na.exclude)
drop1(Evenness_fit,~.,test="F") # type III SS and F Tests
#TESTING ASSUMPTIONS
Evenness_resids <- residuals(Evenness_fit)
Evenness_preds <- predict(Evenness_fit)
plot(Evenness_resids ~ Evenness_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
#plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
#test for homogeneity of variance
#Test for normality
shapiro.test(Evenness_resids)
qqnorm(Evenness_resids)
qqline(Evenness_resids)

#Prepare data for figure
DNARNAdiv.long=melt(div.DNARNA.rfy, id.vars=c("Row.names", "ID", "Shannon Div", "Type", "Treatment"))
DNARNAdiv.long
figure <- ggplot(data=DNARNAdiv.long, aes(x=Treatment, y=value))+
#need to include 'outlier.shape=NA' otherwide the outliers get plotted twice.
geom_boxplot(outlier.shape=NA, aes(fill=Type)) + 
facet_grid(variable~., scales="free_y")+
scale_size(guide=FALSE)+
scale_x_discrete(name="Treatment")+
scale_y_continuous(name="Diversity value")+
theme_bw(base_size=12)

figure





# Analysis Code for All 8 Treatments (and only RNA) -----------------------

#Subset the DNA/RNA relative abundance dataset for ONLY RNA samples.
DNARNA.rfy.ra.OTUs <- otu_table(DNARNA.rfy.ra)
OTUs_RNA_DNARNA.rfy.ra <- DNARNA.rfy.ra.OTUs[ , !c(TRUE,FALSE) ] 
test <-as(otu_table(OTUs_RNA_DNARNA.rfy.ra), "matrix")
test <- as.data.frame(test)

#Convert to matrix, then to phyloseq format.
testRNA <- data.matrix(test)
testRNA = otu_table(testRNA, taxa_are_rows = TRUE)

#Create phyloseq object containing RNA only
RNA.rfy.ra = phyloseq(testRNA, TAX, METADATA)
RNA.rfy.ra

RNA.rfy.ra
sample_sums(RNA.rfy.ra)

#braycurtis
bc <- phyloseq::distance(RNA.rfy.ra, method="bray")


#pcoa plot: normalized weighted Unifrac
set.seed(1)
# Ordinate
RNA_pcoa <- ordinate(physeq = RNA.rfy.ra, method = "PCoA", distance = "bc")
# Plot 
RNA_pcoa <- plot_ordination(physeq = RNA.rfy.ra,ordination = RNA_pcoa,color = "Treatment",title = "Treatment effects on rhizosphere soil communities (16s RNA)") + geom_point(size = 3)  
#Change size of graph
RNA_pcoa + coord_fixed()+ scale_color_manual(name="Treatment", breaks=c("Pre-dry", "Post-dry", "Post-water", "Control", "ABA", "IAA", "JA", "SA"), values=c("red","blue","yellow","orange","magenta4","hotpink","chocolate4","green4"))
ggsave("8trts_RNA_pcoa.pdf",width=178,units="mm") 

#Permanova test using ADONIS function
# make a data frame from the sample_data
RNA.rfy.df <- data.frame(sample_data(RNA.rfy.ra))
RNA.rfy.df
# Adonis test
#(Do their means/centroids differ?)
#make sure the degrees of freedom are correct!
set.seed(1)
adonis(bc ~ Treatment, data = RNA.rfy.df)

#Homogeneity of dispersions test...TESTS ASSUMPTIONS OF ADONIS TEST.
#(do their variances/dispersions differ?)
set.seed(1)
result <- betadisper(bc, RNA.rfy.df$Treatment, type=c("centroid"))
set.seed(1)
permutest(result)






# Analysis Code for Only 5 Treatments (and only RNA) ----------------------

#Subset to only the four hormones plus control:
RNA.rfy.ra = subset_samples(RNA.rfy.ra, Treatment != "Pre-dry" & Treatment != "Post-dry" & Treatment != "Post-water")
RNA.rfy.ra
sample_data(RNA.rfy.ra)

#braycurtis
bc <- phyloseq::distance(RNA.rfy.ra, method="bray")


#pcoa plot: normalized weighted Unifrac
set.seed(1)
# Ordinate
RNA_pcoa <- ordinate(physeq = RNA.rfy.ra, method = "PCoA", distance = "bc")
# Plot 
RNA_pcoa <- plot_ordination(physeq = RNA.rfy.ra,ordination = RNA_pcoa,color = "Treatment",title = "Phytohormone effects on rhizosphere soil communities (16s RNA)") + geom_point(size = 3)  
#Change size of graph
RNA_pcoa + coord_fixed()+ scale_color_manual(name="Treatment", breaks=c("Control", "ABA", "IAA", "JA", "SA"), values=c("red","blue","yellow","orange","green4"))
ggsave("5trts_RNA_pcoa.pdf",width=178,units="mm") 

#Permanova test using ADONIS function
# make a data frame from the sample_data
RNA.rfy.df <- data.frame(sample_data(RNA.rfy.ra))
RNA.rfy.df
# Adonis test
set.seed(1)
adonis(bc ~ Treatment, data = RNA.rfy.df)

#Homogeneity of dispersions test.
set.seed(1)
result <- betadisper(bc, RNA.rfy.df$Treatment, type=c("centroid"))
set.seed(1)
permutest(result)



####
#Subset to only ABA plus control:
RNA.rfy.ra.ABA = subset_samples(RNA.rfy.ra, Treatment == "ABA" | Treatment == "Control")
sample_data(RNA.rfy.ra.ABA)

#braycurtis
bc <- phyloseq::distance(RNA.rfy.ra.ABA, method="bray")


#Permanova test using ADONIS function
# make a data frame from the sample_data
RNA.rfy.ra.df <- data.frame(sample_data(RNA.rfy.ra.ABA))
RNA.rfy.ra.df
# ADONIS test
set.seed(1)
adonis(bc ~ Treatment, data = RNA.rfy.ra.df)

#Homogeneity of dispersions test
set.seed(1)
result <- betadisper(bc, RNA.rfy.ra.df$Treatment, type=c("centroid"))
set.seed(1)
permutest(result)

####
#Subset to only IAA plus control:
RNA.rfy.ra.IAA = subset_samples(RNA.rfy.ra, Treatment == "IAA" | Treatment == "Control")
sample_data(RNA.rfy.ra.IAA)


#braycurtis
bc <- phyloseq::distance(RNA.rfy.ra.IAA, method="bray")

#Permanova test using ADONIS function
# make a data frame from the sample_data
RNA.rfy.ra.df <- data.frame(sample_data(RNA.rfy.ra.IAA))
RNA.rfy.ra.df
# ADONIS test
set.seed(1)
adonis(bc ~ Treatment, data = RNA.rfy.ra.df)

#Homogeneity of dispersions test
set.seed(1)
result <- betadisper(bc, RNA.rfy.ra.df$Treatment, type=c("centroid"))
set.seed(1)
permutest(result)

####
#Subset to only JA plus control:
RNA.rfy.ra.JA = subset_samples(RNA.rfy.ra, Treatment == "JA" | Treatment == "Control")
sample_data(RNA.rfy.ra.JA)


#braycurtis
bc <- phyloseq::distance(RNA.rfy.ra.JA, method="bray")

#Permanova test using ADONIS function
# make a data frame from the sample_data
RNA.rfy.ra.df <- data.frame(sample_data(RNA.rfy.ra.JA))
RNA.rfy.ra.df
# ADONIS test
set.seed(1)
adonis(bc ~ Treatment, data = RNA.rfy.ra.df)

#Homogeneity of dispersions test
set.seed(1)
result <- betadisper(bc, RNA.rfy.ra.df$Treatment, type=c("centroid"))
set.seed(1)
permutest(result)

####
#Subset to only SA plus control:
RNA.rfy.ra.SA = subset_samples(RNA.rfy.ra, Treatment == "SA" | Treatment == "Control")
sample_data(RNA.rfy.ra.SA)

#braycurtis
bc <- phyloseq::distance(RNA.rfy.ra.SA, method="bray")

#Permanova test using ADONIS function
# make a data frame from the sample_data
RNA.rfy.ra.df <- data.frame(sample_data(RNA.rfy.ra.SA))
RNA.rfy.ra.df
# ADONIS test
set.seed(1)
adonis(bc ~ Treatment, data = RNA.rfy.ra.df)

#Homogeneity of dispersions test
set.seed(1)
result <- betadisper(bc, RNA.rfy.ra.df$Treatment, type=c("centroid"))
set.seed(1)
permutest(result)




# Compare RNA rel abund of each OTU between hormone and control -----------

#Prepare input files for t-tests.
RNA.rfy.ra.tab=(as(otu_table(RNA.rfy.ra), "matrix"))
RNA.rfy.ra.map=(sample_data(RNA.rfy.ra))
#assign which treatment belongs to which samples
trt=t(RNA.rfy.ra.map[,"Treatment"])
trt

###ABA vs Control
#Prepare output file into which the t-test results will go.
otuname=row.names(RNA.rfy.ra.tab)
ttest.out.ABA <- data.frame(matrix(nrow = nrow(RNA.rfy.ra.tab), ncol = 4))
colnames(ttest.out.ABA) <- c("OTU", "t-statistic","df", "p.value")
rownames(ttest.out.ABA) <- row.names(RNA.rfy.ra.tab)
#Perform the t-tests.
for(i in 1:length(otuname)){
  #subset the data to test one OTU at a time
  Hormone=RNA.rfy.ra.tab[row.names(RNA.rfy.ra.tab)==otuname[i],trt=="ABA"]
  Control=RNA.rfy.ra.tab[row.names(RNA.rfy.ra.tab)==otuname[i],trt=="Control"]
  #perform the test
  test=t.test(Hormone, Control, paired=FALSE, var.equal = FALSE)
  ttest.out.ABA[i,1] <- row.names(RNA.rfy.ra.tab)[i]
  ttest.out.ABA[i,2] <-(test$statistic)
  ttest.out.ABA[i,3] <-(test$parameter)
  ttest.out.ABA[i,4] <-(test$p.value)
}
#Calculate adjusted p-values and add to table. Sort by adjusted p-values.
ttest.out.ABA[,5] <-p.adjust(ttest.out.ABA$p.value, method = "fdr")
colnames(ttest.out.ABA) <- c("OTU", "t-statistic","df", "p.value", "adj.p.value")
ttest.out.ABA <- ttest.out.ABA[order(ttest.out.ABA$adj.p.value),]
ttest.out.ABA
SigDiffOTUs.ABA <- row.names(subset(ttest.out.ABA, adj.p.value<=0.05))

###IAA vs Control
#Prepare output file into which the t-test results will go.
otuname=row.names(RNA.rfy.ra.tab)
ttest.out.IAA <- data.frame(matrix(nrow = nrow(RNA.rfy.ra.tab), ncol = 4))
colnames(ttest.out.IAA) <- c("OTU", "t-statistic","df", "p.value")
rownames(ttest.out.IAA) <- row.names(RNA.rfy.ra.tab)
#Perform the t-tests.
for(i in 1:length(otuname)){
  #subset the data to test one OTU at a time
  Hormone=RNA.rfy.ra.tab[row.names(RNA.rfy.ra.tab)==otuname[i],trt=="IAA"]
  Control=RNA.rfy.ra.tab[row.names(RNA.rfy.ra.tab)==otuname[i],trt=="Control"]
  #perform the test
  test=t.test(Hormone, Control, paired=FALSE, var.equal = FALSE)
  ttest.out.IAA[i,1] <- row.names(RNA.rfy.ra.tab)[i]
  ttest.out.IAA[i,2] <-(test$statistic)
  ttest.out.IAA[i,3] <-(test$parameter)
  ttest.out.IAA[i,4] <-(test$p.value)
}
#Calculate adjusted p-values and add to table. Sort by adjusted p-values.
ttest.out.IAA[,5] <-p.adjust(ttest.out.IAA$p.value, method = "fdr")
colnames(ttest.out.IAA) <- c("OTU", "t-statistic","df", "p.value", "adj.p.value")
ttest.out.IAA <- ttest.out.IAA[order(ttest.out.IAA$adj.p.value),]
ttest.out.IAA
SigDiffOTUs.IAA <- row.names(subset(ttest.out.IAA, adj.p.value<=0.05))
SigDiffOTUs.IAA

###JA vs Control
#Prepare output file into which the t-test results will go.
otuname=row.names(RNA.rfy.ra.tab)
ttest.out.JA <- data.frame(matrix(nrow = nrow(RNA.rfy.ra.tab), ncol = 4))
colnames(ttest.out.JA) <- c("OTU", "t-statistic","df", "p.value")
rownames(ttest.out.JA) <- row.names(RNA.rfy.ra.tab)
#Perform the t-tests.
for(i in 1:length(otuname)){
  #subset the data to test one OTU at a time
  Hormone=RNA.rfy.ra.tab[row.names(RNA.rfy.ra.tab)==otuname[i],trt=="JA"]
  Control=RNA.rfy.ra.tab[row.names(RNA.rfy.ra.tab)==otuname[i],trt=="Control"]
  #perform the test
  test=t.test(Hormone, Control, paired=FALSE, var.equal = FALSE)
  ttest.out.JA[i,1] <- row.names(RNA.rfy.ra.tab)[i]
  ttest.out.JA[i,2] <-(test$statistic)
  ttest.out.JA[i,3] <-(test$parameter)
  ttest.out.JA[i,4] <-(test$p.value)
}
#Calculate adjusted p-values and add to table. Sort by adjusted p-values.
ttest.out.JA[,5] <-p.adjust(ttest.out.JA$p.value, method = "fdr")
colnames(ttest.out.JA) <- c("OTU", "t-statistic","df", "p.value", "adj.p.value")
ttest.out.JA <- ttest.out.JA[order(ttest.out.JA$adj.p.value),]
ttest.out.JA
SigDiffOTUs.JA <- row.names(subset(ttest.out.JA, adj.p.value<=0.05))
SigDiffOTUs.JA

###SA vs Control
#Prepare output file into which the t-test results will go.
otuname=row.names(RNA.rfy.ra.tab)
ttest.out.SA <- data.frame(matrix(nrow = nrow(RNA.rfy.ra.tab), ncol = 4))
colnames(ttest.out.SA) <- c("OTU", "t-statistic","df", "p.value")
rownames(ttest.out.SA) <- row.names(RNA.rfy.ra.tab)
#Perform the t-tests.
for(i in 1:length(otuname)){
  #subset the data to test one OTU at a time
  Hormone=RNA.rfy.ra.tab[row.names(RNA.rfy.ra.tab)==otuname[i],trt=="SA"]
  Control=RNA.rfy.ra.tab[row.names(RNA.rfy.ra.tab)==otuname[i],trt=="Control"]
  #perform the test
  test=t.test(Hormone, Control, paired=FALSE, var.equal = FALSE)
  ttest.out.SA[i,1] <- row.names(RNA.rfy.ra.tab)[i]
  ttest.out.SA[i,2] <-(test$statistic)
  ttest.out.SA[i,3] <-(test$parameter)
  ttest.out.SA[i,4] <-(test$p.value)
}
#Calculate adjusted p-values and add to table. Sort by adjusted p-values.
ttest.out.SA[,5] <-p.adjust(ttest.out.SA$p.value, method = "fdr")
colnames(ttest.out.SA) <- c("OTU", "t-statistic","df", "p.value", "adj.p.value")
ttest.out.SA <- ttest.out.SA[order(ttest.out.SA$adj.p.value),]
ttest.out.SA
SigDiffOTUs.SA <- row.names(subset(ttest.out.SA, adj.p.value<=0.05))

#Now prepare relative abundance data to plot with t-test results.

#Get relative abundance values for only specific treatments.
OTUs.IAA.tab=RNA.rfy.ra.tab[,trt=="IAA"]
OTUs.JA.tab=RNA.rfy.ra.tab[,trt=="JA"]
OTUs.SA.tab=RNA.rfy.ra.tab[,trt=="SA"]
OTUs.Control.tab=RNA.rfy.ra.tab[,trt=="Control"]

#Calculate  mean taxa level rel. abundance within each treatment
OTUs.IAA.mean=apply(OTUs.IAA.tab,1,mean)
OTUs.JA.mean=apply(OTUs.JA.tab,1,mean)
OTUs.SA.mean=apply(OTUs.SA.tab,1,mean)
OTUs.Control.mean=apply(OTUs.Control.tab,1,mean)

#Calculate difference between hormone and control

OTUs.IAA.diff=as.data.frame(OTUs.IAA.mean-OTUs.Control.mean)
#Keep only OTUs that significantly differed from control
SigDiffOTUs.IAA.diff <-subset(OTUs.IAA.diff, row.names(OTUs.IAA.diff)%in%SigDiffOTUs.IAA)
SigDiffOTUs.IAA.diff=as.data.frame(SigDiffOTUs.IAA.diff)
SigDiffOTUs.IAA.diff

OTUs.JA.diff=as.data.frame(OTUs.JA.mean-OTUs.Control.mean)
#Keep only OTUs that significantly differed from control
SigDiffOTUs.JA.diff <-subset(OTUs.JA.diff, row.names(OTUs.JA.diff)%in%SigDiffOTUs.JA)
SigDiffOTUs.JA.diff=as.data.frame(SigDiffOTUs.JA.diff)
SigDiffOTUs.JA.diff

OTUs.SA.diff=as.data.frame(OTUs.SA.mean-OTUs.Control.mean)
#Keep only OTUs that significantly differed from control
SigDiffOTUs.SA.diff <-subset(OTUs.SA.diff, row.names(OTUs.SA.diff)%in%SigDiffOTUs.SA)
SigDiffOTUs.SA.diff=as.data.frame(SigDiffOTUs.SA.diff)
SigDiffOTUs.SA.diff

#Plotting relative abundances of significant OTUs.

colnames(SigDiffOTUs.IAA.diff)=c("IAA.Control")
SigDiffOTUs.IAA.diff$OTU <- rownames(SigDiffOTUs.IAA.diff)
SigDiffOTUs.IAA.diff
ggplot(SigDiffOTUs.IAA.diff, aes(x=OTU, y=as.numeric(IAA.Control))) +
  geom_bar(stat='identity',fill="yellow") +
  coord_flip()+
  labs(x="",y="Change in Rel Abund (IAA-Control)")+
  geom_hline(yintercept=0, aes(size=0.25))
  #ylim(-0.015,0.015)
ggsave("OTU.IAA.diff.pdf",units="mm") 

colnames(SigDiffOTUs.JA.diff)=c("JA.Control")
SigDiffOTUs.JA.diff$OTU <- rownames(SigDiffOTUs.JA.diff)
SigDiffOTUs.JA.diff
ggplot(SigDiffOTUs.JA.diff, aes(x=OTU, y=JA.Control)) +
  geom_bar(stat='identity',fill="orange") +
  coord_flip()+
  labs(x="",y="Change in Rel Abund (JA-Control)")+
  geom_hline(yintercept=0, aes(size=0.25))
  #ylim(-0.13,0.13)
ggsave("OTU.JA.diff.pdf",units="mm") 


colnames(SigDiffOTUs.SA.diff)=c("SA.Control")
SigDiffOTUs.SA.diff$OTU <- rownames(SigDiffOTUs.SA.diff)
SigDiffOTUs.SA.diff
ggplot(SigDiffOTUs.SA.diff, aes(x=OTU, y=as.numeric(SA.Control))) +
  geom_bar(stat='identity',fill="green4") +
  coord_flip()+
  labs(x="",y="Change in Rel Abund (SA-Control)")+
  geom_hline(yintercept=0, aes(size=0.25))
ggsave("OTU.SA.diff.pdf",units="mm") 





# Compare RNA rel abund of each PHYLUM between hormone and control--------

#For each sample, find mean relative abundance by a Given Taxonomy
RNA.rfy.ra.dat <- otu_table(RNA.rfy.ra)
RNA.rfy.tax.dat <- tax_table(RNA.rfy.ra)

#Here, choose the TAXONOMIC HEIRARCHY AT WHICH YOU WANT TO DO THE TTESTS.
RNA.rfy.tax.data <- subset(RNA.rfy.tax.dat, select = c("Phylum"))

#Now build the data.frame holding relative abundance and taxonomy data.
RNA.tax.ra.data <- merge(RNA.rfy.ra.dat, RNA.rfy.tax.data, by="row.names")
row.names(RNA.tax.ra.data) <- RNA.tax.ra.data$Row.names
RNA.tax.ra.data <- subset(RNA.tax.ra.data, select = -c(Row.names))

#Sum within each row by the Taxonomy of Interest.
library("dplyr")
RNA.tax.summary <-RNA.tax.ra.data %>% group_by(Phylum) %>% summarise_all(funs(sum))
RNA.tax.summary <- as.data.frame(RNA.tax.summary)
row.names(RNA.tax.summary) <- RNA.tax.summary$Phylum
RNA.tax.summary <- subset(RNA.tax.summary, select = -c(Phylum))

#Prepare input files for t-tests.
RNA.rfy.ra.map=(sample_data(RNA.rfy.ra))
#assign which treatment belongs to which samples
trt=t(RNA.rfy.ra.map[,"Treatment"])
trt

#ABA vs Control
#Prepare output file into which the t-test results will go.
TAXAname=row.names(RNA.tax.summary)
ttest.out.ABA <- data.frame(matrix(nrow = nrow(RNA.tax.summary), ncol = 4))
colnames(ttest.out.ABA) <- c("OTU", "t-statistic","df", "p.value")
rownames(ttest.out.ABA) <- row.names(RNA.tax.summary)
#Perform the t-tests.
for(i in 1:length(TAXAname)){
  #subset the data to test one OTU at a time
  Hormone=RNA.tax.summary[row.names(RNA.tax.summary)==TAXAname[i],trt=="ABA"]
  Control=RNA.tax.summary[row.names(RNA.tax.summary)==TAXAname[i],trt=="Control"]
  #perform the test
  test=t.test(Hormone, Control, paired=FALSE, var.equal = FALSE)
  ttest.out.ABA[i,1] <- row.names(RNA.tax.summary)[i]
  ttest.out.ABA[i,2] <-(test$statistic)
  ttest.out.ABA[i,3] <-(test$parameter)
  ttest.out.ABA[i,4] <-(test$p.value)
}
#Calculate adjusted p-values and add to table. Sort by adjusted p-values.
ttest.out.ABA[,5] <-p.adjust(ttest.out.ABA$p.value, method = "fdr")
colnames(ttest.out.ABA) <- c("OTU", "t-statistic","df", "p.value", "adj.p.value")
ttest.out.ABA <- ttest.out.ABA[order(ttest.out.ABA$adj.p.value),]
ttest.out.ABA
SigDiffTAXA.ABA <- row.names(subset(ttest.out.ABA, adj.p.value<=0.05))

#IAA vs Control
#Prepare output file into which the t-test results will go.
TAXAname=row.names(RNA.tax.summary)
ttest.out.IAA <- data.frame(matrix(nrow = nrow(RNA.tax.summary), ncol = 4))
colnames(ttest.out.IAA) <- c("OTU", "t-statistic","df", "p.value")
rownames(ttest.out.IAA) <- row.names(RNA.tax.summary)
#Perform the t-tests.
for(i in 1:length(TAXAname)){
  #subset the data to test one OTU at a time
  Hormone=RNA.tax.summary[row.names(RNA.tax.summary)==TAXAname[i],trt=="IAA"]
  Control=RNA.tax.summary[row.names(RNA.tax.summary)==TAXAname[i],trt=="Control"]
  #perform the test
  test=t.test(Hormone, Control, paired=FALSE, var.equal = FALSE)
  ttest.out.IAA[i,1] <- row.names(RNA.tax.summary)[i]
  ttest.out.IAA[i,2] <-(test$statistic)
  ttest.out.IAA[i,3] <-(test$parameter)
  ttest.out.IAA[i,4] <-(test$p.value)
}
#Calculate adjusted p-values and add to table. Sort by adjusted p-values.
ttest.out.IAA[,5] <-p.adjust(ttest.out.IAA$p.value, method = "fdr")
colnames(ttest.out.IAA) <- c("OTU", "t-statistic","df", "p.value", "adj.p.value")
ttest.out.IAA <- ttest.out.IAA[order(ttest.out.IAA$adj.p.value),]
ttest.out.IAA
SigDiffTAXA.IAA <- row.names(subset(ttest.out.IAA, adj.p.value<=0.05))

#JA vs Control
#Prepare output file into which the t-test results will go.
TAXAname=row.names(RNA.tax.summary)
ttest.out.JA <- data.frame(matrix(nrow = nrow(RNA.tax.summary), ncol = 4))
colnames(ttest.out.JA) <- c("OTU", "t-statistic","df", "p.value")
rownames(ttest.out.JA) <- row.names(RNA.tax.summary)
#Perform the t-tests.
for(i in 1:length(TAXAname)){
  #subset the data to test one OTU at a time
  Hormone=RNA.tax.summary[row.names(RNA.tax.summary)==TAXAname[i],trt=="JA"]
  Control=RNA.tax.summary[row.names(RNA.tax.summary)==TAXAname[i],trt=="Control"]
  #perform the test
  test=t.test(Hormone, Control, paired=FALSE, var.equal = FALSE)
  ttest.out.JA[i,1] <- row.names(RNA.tax.summary)[i]
  ttest.out.JA[i,2] <-(test$statistic)
  ttest.out.JA[i,3] <-(test$parameter)
  ttest.out.JA[i,4] <-(test$p.value)
}
#Calculate adjusted p-values and add to table. Sort by adjusted p-values.
ttest.out.JA[,5] <-p.adjust(ttest.out.JA$p.value, method = "fdr")
colnames(ttest.out.JA) <- c("OTU", "t-statistic","df", "p.value", "adj.p.value")
ttest.out.JA <- ttest.out.JA[order(ttest.out.JA$adj.p.value),]
ttest.out.JA
SigDiffTAXA.JA <- row.names(subset(ttest.out.JA, adj.p.value<=0.05))

#SA vs Control
#Prepare output file into which the t-test results will go.
TAXAname=row.names(RNA.tax.summary)
ttest.out.SA <- data.frame(matrix(nrow = nrow(RNA.tax.summary), ncol = 4))
colnames(ttest.out.SA) <- c("OTU", "t-statistic","df", "p.value")
rownames(ttest.out.SA) <- row.names(RNA.tax.summary)
#Perform the t-tests.
for(i in 1:length(TAXAname)){
  #subset the data to test one OTU at a time
  Hormone=RNA.tax.summary[row.names(RNA.tax.summary)==TAXAname[i],trt=="SA"]
  Control=RNA.tax.summary[row.names(RNA.tax.summary)==TAXAname[i],trt=="Control"]
  #perform the test
  test=t.test(Hormone, Control, paired=FALSE, var.equal = FALSE)
  ttest.out.SA[i,1] <- row.names(RNA.tax.summary)[i]
  ttest.out.SA[i,2] <-(test$statistic)
  ttest.out.SA[i,3] <-(test$parameter)
  ttest.out.SA[i,4] <-(test$p.value)
}
#Calculate adjusted p-values and add to table. Sort by adjusted p-values.
ttest.out.SA[,5] <-p.adjust(ttest.out.SA$p.value, method = "fdr")
colnames(ttest.out.SA) <- c("OTU", "t-statistic","df", "p.value", "adj.p.value")
ttest.out.SA <- ttest.out.SA[order(ttest.out.SA$adj.p.value),]
ttest.out.SA
SigDiffTAXA.SA <- row.names(subset(ttest.out.SA, adj.p.value<=0.0501))

#Now prepare relative abundance data to plot with t-test results.

#Make separate datasets for each hormone and control.
TAXA.ABA.tab=RNA.tax.summary[,trt=="ABA"]
TAXA.IAA.tab=RNA.tax.summary[,trt=="IAA"]
TAXA.JA.tab=RNA.tax.summary[,trt=="JA"]
TAXA.SA.tab=RNA.tax.summary[,trt=="SA"]
TAXA.Control.tab=RNA.tax.summary[,trt=="Control"]

#Calculate  mean taxa level rel. abundance within each treatment
TAXA.ABA.mean=apply(TAXA.ABA.tab,1,mean)
TAXA.IAA.mean=apply(TAXA.IAA.tab,1,mean)
TAXA.JA.mean=apply(TAXA.JA.tab,1,mean)
TAXA.SA.mean=apply(TAXA.SA.tab,1,mean)
TAXA.Control.mean=apply(TAXA.Control.tab,1,mean)

#Calculate difference between hormone and control
TAXA.ABA.diff=as.data.frame(TAXA.ABA.mean-TAXA.Control.mean)
#Keep only TAXA that significantly differed from control
TAXA.ABA.diff <-subset(TAXA.ABA.diff, row.names(TAXA.ABA.diff)%in%SigDiffTAXA.ABA)
TAXA.ABA.diff=as.data.frame(TAXA.ABA.diff)
TAXA.ABA.diff

#Calculate difference between hormone and control
TAXA.IAA.diff=as.data.frame(TAXA.IAA.mean-TAXA.Control.mean)
#Keep only TAXA that significantly differed from control
TAXA.IAA.diff <-subset(TAXA.IAA.diff, row.names(TAXA.IAA.diff)%in%SigDiffTAXA.IAA)
TAXA.IAA.diff=as.data.frame(TAXA.IAA.diff)
TAXA.IAA.diff

#Calculate difference between hormone and control
TAXA.JA.diff=as.data.frame(TAXA.JA.mean-TAXA.Control.mean)
#Keep only TAXA that significantly differed from control
TAXA.JA.diff <-subset(TAXA.JA.diff, row.names(TAXA.JA.diff)%in%SigDiffTAXA.JA)
TAXA.JA.diff=as.data.frame(TAXA.JA.diff)
TAXA.JA.diff

#Calculate difference between hormone and control
TAXA.SA.diff=as.data.frame(TAXA.SA.mean-TAXA.Control.mean)
#Keep only TAXA that significantly differed from control
TAXA.SA.diff <-subset(TAXA.SA.diff, row.names(TAXA.SA.diff)%in%SigDiffTAXA.SA)
TAXA.SA.diff=as.data.frame(TAXA.SA.diff)
TAXA.SA.diff

#Now plot the relative abundance changes for phyla with the largest changes.

####ABA
colnames(TAXA.ABA.diff)=c("ABA.Control")
TAXA.ABA.diff$Phylum <- rownames(TAXA.ABA.diff)
TAXA.ABA.diff=TAXA.ABA.diff[order(abs((TAXA.ABA.diff$ABA.Control)),decreasing=TRUE),]
TAXA.ABA.diff=as.data.frame(TAXA.ABA.diff)
TAXA.ABA.diff
TAXA.ABA.diff.top10=head(TAXA.ABA.diff,10)
ggplot(TAXA.ABA.diff.top10, aes(x=reorder(Phylum, abs(ABA.Control)), y=ABA.Control)) +
  geom_bar(stat='identity',fill="red") +
  coord_flip()+
  labs(x="",y="Change in Rel Abund (ABA-Control)")+
  geom_hline(yintercept=0, aes(size=0.25))
ggsave("Phylum.ABA.pdf",units="mm") 

####IAA
colnames(TAXA.IAA.diff)=c("IAA.Control")
TAXA.IAA.diff$Phylum <- rownames(TAXA.IAA.diff)
TAXA.IAA.diff=TAXA.IAA.diff[order(abs((TAXA.IAA.diff$IAA.Control)),decreasing=TRUE),]
TAXA.IAA.diff=as.data.frame(TAXA.IAA.diff)
TAXA.IAA.diff.top10=head(TAXA.IAA.diff,10)
ggplot(TAXA.IAA.diff.top10, aes(x=reorder(Phylum, abs(IAA.Control)), y=IAA.Control)) +
  geom_bar(stat='identity',fill="yellow") +
  coord_flip()+
  labs(x="",y="Change in Rel Abund (IAA-Control)")+
  geom_hline(yintercept=0, aes(size=0.25))
ggsave("Phylum.IAA.pdf",units="mm") 

####JA
colnames(TAXA.JA.diff)=c("JA.Control")
TAXA.JA.diff$Phylum <- rownames(TAXA.JA.diff)
TAXA.JA.diff=TAXA.JA.diff[order(abs((TAXA.JA.diff$JA.Control)),decreasing=TRUE),]
TAXA.JA.diff=as.data.frame(TAXA.JA.diff)
TAXA.JA.diff
TAXA.JA.diff.top10=head(TAXA.JA.diff,11)
ggplot(TAXA.JA.diff.top10, aes(x=reorder(Phylum, abs(JA.Control)), y=JA.Control)) +
  geom_bar(stat='identity',fill="orange") +
  coord_flip()+
  labs(x="",y="Change in Rel Abund (JA-Control)")+
  geom_hline(yintercept=0, aes(size=0.25))
ggsave("Phylum.JA.pdf",units="mm") 

####SA
colnames(TAXA.SA.diff)=c("SA.Control")
TAXA.SA.diff$Phylum <- rownames(TAXA.SA.diff)
TAXA.SA.diff=TAXA.SA.diff[order(abs((TAXA.SA.diff$SA.Control)),decreasing=TRUE),]
TAXA.SA.diff=as.data.frame(TAXA.SA.diff)
TAXA.SA.diff
TAXA.SA.diff.top10=head(TAXA.SA.diff,10)
ggplot(TAXA.SA.diff.top10, aes(x=reorder(Phylum, abs(SA.Control)), y=SA.Control)) +
  geom_bar(stat='identity',fill="green4") +
  coord_flip()+
  labs(x="",y="Change in Rel Abund (SA-Control)")+
  geom_hline(yintercept=0, aes(size=0.25))
ggsave("Phylum.SA.pdf",units="mm") 



















# Compare RNA rel abund of each CLASS/FAMILY/ORDER between hormone and control ---------

#Optional: Subset to include only Proteobacteria
#RNA.rfy.ra <- subset_taxa(RNA.rfy.ra, Phylum=="Proteobacteria")
#RNA.rfy.ra

#For each sample, find mean relative abundance by a Given Taxonomy
RNA.rfy.ra.dat <- otu_table(RNA.rfy.ra)
RNA.rfy.tax.dat <- tax_table(RNA.rfy.ra)

#Here, choose the TAXONOMIC HEIRARCHY AT WHICH YOU WANT TO DO THE TTESTS.
RNA.rfy.tax.data <- subset(RNA.rfy.tax.dat, select = c("Family"))
RNA.rfy.tax.data
#Now build the data.frame holding relative abundance and taxonomy data.
RNA.tax.ra.data <- merge(RNA.rfy.ra.dat, RNA.rfy.tax.data, by="row.names")
row.names(RNA.tax.ra.data) <- RNA.tax.ra.data$Row.names
RNA.tax.ra.data <- subset(RNA.tax.ra.data, select = -c(Row.names))

#Sum within each row by the Taxonomy of Interest.
library("dplyr")
RNA.tax.summary <-RNA.tax.ra.data %>% group_by(Family) %>% summarise_all(funs(sum))
RNA.tax.summary <- as.data.frame(RNA.tax.summary)
row.names(RNA.tax.summary) <- RNA.tax.summary$Family
RNA.tax.summary <- subset(RNA.tax.summary, select = -c(Family))

#Prepare input files for t-tests.
RNA.rfy.ra.map=(sample_data(RNA.rfy.ra))
#assign which treatment belongs to which samples
trt=t(RNA.rfy.ra.map[,"Treatment"])
trt

#ABA vs Control
#Prepare output file into which the t-test results will go.
TAXAname=row.names(RNA.tax.summary)
ttest.out.ABA <- data.frame(matrix(nrow = nrow(RNA.tax.summary), ncol = 4))
colnames(ttest.out.ABA) <- c("OTU", "t-statistic","df", "p.value")
rownames(ttest.out.ABA) <- row.names(RNA.tax.summary)
#Perform the t-tests.
for(i in 1:length(TAXAname)){
  #subset the data to test one OTU at a time
  Hormone=RNA.tax.summary[row.names(RNA.tax.summary)==TAXAname[i],trt=="ABA"]
  Control=RNA.tax.summary[row.names(RNA.tax.summary)==TAXAname[i],trt=="Control"]
  #perform the test
  test=t.test(Hormone, Control, paired=FALSE, var.equal = FALSE)
  ttest.out.ABA[i,1] <- row.names(RNA.tax.summary)[i]
  ttest.out.ABA[i,2] <-(test$statistic)
  ttest.out.ABA[i,3] <-(test$parameter)
  ttest.out.ABA[i,4] <-(test$p.value)
}
#Calculate adjusted p-values and add to table. Sort by adjusted p-values.
ttest.out.ABA[,5] <-p.adjust(ttest.out.ABA$p.value, method = "fdr")
colnames(ttest.out.ABA) <- c("OTU", "t-statistic","df", "p.value", "adj.p.value")
ttest.out.ABA <- ttest.out.ABA[order(ttest.out.ABA$adj.p.value),]
ttest.out.ABA
SigDiffTAXA.ABA <- row.names(subset(ttest.out.ABA, adj.p.value<=0.05))

#IAA vs Control
#Prepare output file into which the t-test results will go.
TAXAname=row.names(RNA.tax.summary)
ttest.out.IAA <- data.frame(matrix(nrow = nrow(RNA.tax.summary), ncol = 4))
colnames(ttest.out.IAA) <- c("OTU", "t-statistic","df", "p.value")
rownames(ttest.out.IAA) <- row.names(RNA.tax.summary)
#Perform the t-tests.
for(i in 1:length(TAXAname)){
  #subset the data to test one OTU at a time
  Hormone=RNA.tax.summary[row.names(RNA.tax.summary)==TAXAname[i],trt=="IAA"]
  Control=RNA.tax.summary[row.names(RNA.tax.summary)==TAXAname[i],trt=="Control"]
  #perform the test
  test=t.test(Hormone, Control, paired=FALSE, var.equal = FALSE)
  ttest.out.IAA[i,1] <- row.names(RNA.tax.summary)[i]
  ttest.out.IAA[i,2] <-(test$statistic)
  ttest.out.IAA[i,3] <-(test$parameter)
  ttest.out.IAA[i,4] <-(test$p.value)
}
#Calculate adjusted p-values and add to table. Sort by adjusted p-values.
ttest.out.IAA[,5] <-p.adjust(ttest.out.IAA$p.value, method = "fdr")
colnames(ttest.out.IAA) <- c("OTU", "t-statistic","df", "p.value", "adj.p.value")
ttest.out.IAA <- ttest.out.IAA[order(ttest.out.IAA$adj.p.value),]
ttest.out.IAA
SigDiffTAXA.IAA <- row.names(subset(ttest.out.IAA, adj.p.value<=0.05))

#JA vs Control
#Prepare output file into which the t-test results will go.
TAXAname=row.names(RNA.tax.summary)
ttest.out.JA <- data.frame(matrix(nrow = nrow(RNA.tax.summary), ncol = 4))
colnames(ttest.out.JA) <- c("OTU", "t-statistic","df", "p.value")
rownames(ttest.out.JA) <- row.names(RNA.tax.summary)
#Perform the t-tests.
for(i in 1:length(TAXAname)){
  #subset the data to test one OTU at a time
  Hormone=RNA.tax.summary[row.names(RNA.tax.summary)==TAXAname[i],trt=="JA"]
  Control=RNA.tax.summary[row.names(RNA.tax.summary)==TAXAname[i],trt=="Control"]
  #perform the test
  test=t.test(Hormone, Control, paired=FALSE, var.equal = FALSE)
  ttest.out.JA[i,1] <- row.names(RNA.tax.summary)[i]
  ttest.out.JA[i,2] <-(test$statistic)
  ttest.out.JA[i,3] <-(test$parameter)
  ttest.out.JA[i,4] <-(test$p.value)
}
#Calculate adjusted p-values and add to table. Sort by adjusted p-values.
ttest.out.JA[,5] <-p.adjust(ttest.out.JA$p.value, method = "fdr")
colnames(ttest.out.JA) <- c("OTU", "t-statistic","df", "p.value", "adj.p.value")
ttest.out.JA <- ttest.out.JA[order(ttest.out.JA$adj.p.value),]
ttest.out.JA
SigDiffTAXA.JA <- row.names(subset(ttest.out.JA, adj.p.value<=0.05))

#SA vs Control
#Prepare output file into which the t-test results will go.
TAXAname=row.names(RNA.tax.summary)
ttest.out.SA <- data.frame(matrix(nrow = nrow(RNA.tax.summary), ncol = 4))
colnames(ttest.out.SA) <- c("OTU", "t-statistic","df", "p.value")
rownames(ttest.out.SA) <- row.names(RNA.tax.summary)
#Perform the t-tests.
for(i in 1:length(TAXAname)){
  #subset the data to test one OTU at a time
  Hormone=RNA.tax.summary[row.names(RNA.tax.summary)==TAXAname[i],trt=="SA"]
  Control=RNA.tax.summary[row.names(RNA.tax.summary)==TAXAname[i],trt=="Control"]
  #perform the test
  test=t.test(Hormone, Control, paired=FALSE, var.equal = FALSE)
  ttest.out.SA[i,1] <- row.names(RNA.tax.summary)[i]
  ttest.out.SA[i,2] <-(test$statistic)
  ttest.out.SA[i,3] <-(test$parameter)
  ttest.out.SA[i,4] <-(test$p.value)
}
#Calculate adjusted p-values and add to table. Sort by adjusted p-values.
ttest.out.SA[,5] <-p.adjust(ttest.out.SA$p.value, method = "fdr")
colnames(ttest.out.SA) <- c("OTU", "t-statistic","df", "p.value", "adj.p.value")
ttest.out.SA <- ttest.out.SA[order(ttest.out.SA$adj.p.value),]
ttest.out.SA
SigDiffTAXA.SA <- row.names(subset(ttest.out.SA, adj.p.value<=0.0501))

#Now prepare relative abundance data to plot with t-test results.

#Make separate datasets for each hormone and control.
TAXA.ABA.tab=RNA.tax.summary[,trt=="ABA"]
TAXA.IAA.tab=RNA.tax.summary[,trt=="IAA"]
TAXA.JA.tab=RNA.tax.summary[,trt=="JA"]
TAXA.SA.tab=RNA.tax.summary[,trt=="SA"]
TAXA.Control.tab=RNA.tax.summary[,trt=="Control"]

#Calculate  mean taxa level rel. abundance within each treatment
TAXA.ABA.mean=apply(TAXA.ABA.tab,1,mean)
TAXA.IAA.mean=apply(TAXA.IAA.tab,1,mean)
TAXA.JA.mean=apply(TAXA.JA.tab,1,mean)
TAXA.SA.mean=apply(TAXA.SA.tab,1,mean)
TAXA.Control.mean=apply(TAXA.Control.tab,1,mean)

#Calculate difference between hormone and control
TAXA.ABA.diff=as.data.frame(TAXA.ABA.mean-TAXA.Control.mean)
#Keep only TAXA that significantly differed from control
SigDiffTAXA.ABA.diff <-subset(TAXA.ABA.diff, row.names(TAXA.ABA.diff)%in%SigDiffTAXA.ABA)
SigDiffTAXA.ABA.diff=as.data.frame(SigDiffTAXA.ABA.diff)
TAXA.ABA.diff <- SigDiffTAXA.ABA.diff

#Calculate difference between hormone and control
TAXA.IAA.diff=as.data.frame(TAXA.IAA.mean-TAXA.Control.mean)
#Keep only TAXA that significantly differed from control
SigDiffTAXA.IAA.diff <-subset(TAXA.IAA.diff, row.names(TAXA.IAA.diff)%in%SigDiffTAXA.IAA)
SigDiffTAXA.IAA.diff=as.data.frame(SigDiffTAXA.IAA.diff)
TAXA.IAA.diff <-SigDiffTAXA.IAA.diff

TAXA.JA.diff=as.data.frame(TAXA.JA.mean-TAXA.Control.mean)
#Keep only TAXA that significantly differed from control
SigDiffTAXA.JA.diff <-subset(TAXA.JA.diff, row.names(TAXA.JA.diff)%in%SigDiffTAXA.JA)
SigDiffTAXA.JA.diff=as.data.frame(SigDiffTAXA.JA.diff)
TAXA.JA.diff <-SigDiffTAXA.JA.diff

TAXA.SA.diff=as.data.frame(TAXA.SA.mean-TAXA.Control.mean)
#Keep only TAXA that significantly differed from control
SigDiffTAXA.SA.diff <-subset(TAXA.SA.diff, row.names(TAXA.SA.diff)%in%SigDiffTAXA.SA)
SigDiffTAXA.SA.diff=as.data.frame(SigDiffTAXA.SA.diff)
TAXA.SA.diff<- SigDiffTAXA.SA.diff

#Now plot the relative abundances of the Order with largest changes.

####ABA
colnames(TAXA.ABA.diff)=c("ABA.Control")
TAXA.ABA.diff$Class <- rownames(TAXA.ABA.diff)
TAXA.ABA.diff=TAXA.ABA.diff[order(abs((TAXA.ABA.diff$ABA.Control)),decreasing=TRUE),]
TAXA.ABA.diff=as.data.frame(TAXA.ABA.diff)
TAXA.ABA.diff.top10=head(TAXA.ABA.diff,10)
ggplot(TAXA.ABA.diff.top10, aes(x=Class, y=ABA.Control)) +
  geom_bar(stat='identity', fill="red") +
  coord_flip()+
  labs(x="",y="Change in Rel Abund (ABA-Control)")+
  geom_hline(yintercept=0, aes(size=0.25))
ggsave("Class.ABA.diff.pdf",units="mm") 

####TAXA.IAA
colnames(TAXA.IAA.diff)=c("TAXA.IAA.Control")
TAXA.IAA.diff$Class <- rownames(TAXA.IAA.diff)
TAXA.IAA.diff=TAXA.IAA.diff[order(abs((TAXA.IAA.diff$TAXA.IAA.Control)),decreasing=TRUE),]
TAXA.IAA.diff=as.data.frame(TAXA.IAA.diff)
TAXA.IAA.diff.top10=head(TAXA.IAA.diff,10)
ggplot(TAXA.IAA.diff.top10, aes(x=Class, y=TAXA.IAA.Control)) +
  geom_bar(stat='identity', fill="yellow") +
  coord_flip()+
  labs(x="",y="Change in Rel Abund (IAA-Control)")+
  geom_hline(yintercept=0, aes(size=0.25))+
  ylim(-0.06,0.02)
ggsave("Class.IAA.diff.pdf",units="mm") 

####TAXA.JA
colnames(TAXA.JA.diff)=c("TAXA.JA.Control")
TAXA.JA.diff$Class <- rownames(TAXA.JA.diff)
TAXA.JA.diff=TAXA.JA.diff[order(abs((TAXA.JA.diff$TAXA.JA.Control)),decreasing=TRUE),]
TAXA.JA.diff=as.data.frame(TAXA.JA.diff)
TAXA.JA.diff.top10=head(TAXA.JA.diff,10)
ggplot(TAXA.JA.diff.top10, aes(x=Class, y=TAXA.JA.Control)) +
  geom_bar(stat='identity',fill="orange") +
  coord_flip()+
  labs(x="",y="Change in Rel Abund (JA-Control)")+
  geom_hline(yintercept=0, aes(size=0.25))
ggsave("Class.JA.diff.pdf",units="mm") 

####TAXA.SA
colnames(TAXA.SA.diff)=c("TAXA.SA.Control")
TAXA.SA.diff$Class <- rownames(TAXA.SA.diff)
TAXA.SA.diff=TAXA.SA.diff[order(abs((TAXA.SA.diff$TAXA.SA.Control)),decreasing=TRUE),]
TAXA.SA.diff=as.data.frame(TAXA.SA.diff)
TAXA.SA.diff.top10=head(TAXA.SA.diff,10)
ggplot(TAXA.SA.diff.top10, aes(x=Class, y=TAXA.SA.Control))+
  geom_bar(stat='identity', fill="green4") +
  coord_flip()+
  labs(x="",y="Change in Rel Abund (SA-Control)")+
  geom_hline(yintercept=0, aes(size=0.25))
ggsave("Class.SA.diff.pdf",units="mm") 









# Analysis Code for All 8 Treatments (and only DNA) -----------------------

#Subset the DNA/RNA relative abundance dataset for ONLY DNA samples.
DNARNA.rfy.ra.OTUs <- otu_table(DNARNA.rfy.ra)
OTUs_DNA_DNARNA.rfy.ra <- DNARNA.rfy.ra.OTUs[ , !c(FALSE,TRUE) ] 
test <-as(otu_table(OTUs_DNA_DNARNA.rfy.ra), "matrix")
test <- as.data.frame(test)

#Convert to matrix, then to phyloseq format.
testDNA <- data.matrix(test)
testDNA = otu_table(testDNA, taxa_are_rows = TRUE)

#Create phyloseq object containing DNA only
DNA.rfy.ra = phyloseq(testDNA, TAX, METADATA)
DNA.rfy.ra

sample_sums(DNA.rfy.ra)


#braycurtis
bc <- phyloseq::distance(DNA.rfy.ra, method="bray")

#pcoa plot: normalized weighted Unifrac
set.seed(1)
# Ordinate
DNA_pcoa <- ordinate(physeq = DNA.rfy.ra, method = "PCoA", distance = "bc")
# Plot 
DNA_pcoa <- plot_ordination(physeq = DNA.rfy.ra,ordination = DNA_pcoa,color = "Treatment",title = "Treatment effects on rhizosphere soil communities (16s DNA)") + geom_point(size = 3)  
#Change size of graph
DNA_pcoa + coord_fixed()+ scale_color_manual(name="Treatment", breaks=c("Pre-dry", "Post-dry", "Post-water", "Control", "ABA", "IAA", "JA", "SA"), values=c("red","blue","yellow","orange","magenta4","hotpink","chocolate4","green4"))
ggsave("8trts_DNA_pcoa.pdf",width=178,units="mm") 

#Permanova test using ADONIS function
# make a data frame from the sample_data
DNA.rfy.df <- data.frame(sample_data(DNA.rfy.ra))
DNA.rfy.df
# Adonis test
set.seed(1)
adonis(bc ~ Treatment, data = DNA.rfy.df)

#Homogeneity of dispersions test
set.seed(1)
result <- betadisper(bc, DNA.rfy.df$Treatment, type=c("centroid"))
set.seed(1)
permutest(result)






# Analysis Code for Only 5 Treatments (and only DNA) ----------------------

#Subset to only the four hormones plus control:
DNA.rfy.ra = subset_samples(DNA.rfy.ra, Treatment != "Pre-dry" & Treatment != "Post-dry" & Treatment != "Post-water")
sample_data(DNA.rfy.ra)


#braycurtis
bc <- phyloseq::distance(DNA.rfy.ra, method="bray")

#pcoa plot: normalized weighted Unifrac
set.seed(1)
# Ordinate
DNA_pcoa <- ordinate(physeq = DNA.rfy.ra, method = "PCoA", distance = "bc")
# Plot 
DNA_pcoa <- plot_ordination(physeq = DNA.rfy.ra,ordination = DNA_pcoa,color = "Treatment",title = "Phytohormone effects on rhizosphere soil communities (16s DNA)") + geom_point(size = 3)  
#Change size of graph
DNA_pcoa + coord_fixed()+ scale_color_manual(name="Treatment", breaks=c("Control", "ABA", "IAA", "JA", "SA"), values=c("red","blue","yellow","orange","green4"))
ggsave("5trts_DNA_pcoa.pdf",width=178,units="mm") 

#Permanova test using ADONIS function
# make a data frame from the sample_data
DNA.rfy.df <- data.frame(sample_data(DNA.rfy.ra))
DNA.rfy.df
# Adonis test
set.seed(1)
adonis(bc ~ Treatment, data = DNA.rfy.df)

#Homogeneity of dispersions test
set.seed(1)
result <- betadisper(bc, DNA.rfy.df$Treatment, type=c("centroid"))
set.seed(1)
permutest(result)






####
#Subset to only ABA plus control:
DNA.rfy.ra.ABA = subset_samples(DNA.rfy.ra, Treatment == "ABA" | Treatment == "Control")
sample_data(DNA.rfy.ra.ABA)


#braycurtis
bc <- phyloseq::distance(DNA.rfy.ra.ABA, method="bray")

#Permanova test using ADONIS function
# make a data frame from the sample_data
DNA.rfy.ra.df <- data.frame(sample_data(DNA.rfy.ra.ABA))
DNA.rfy.ra.df
# ADONIS test
set.seed(1)
adonis(bc ~ Treatment, data = DNA.rfy.ra.df)

#Homogeneity of dispersions test
set.seed(1)
result <- betadisper(bc, DNA.rfy.ra.df$Treatment, type=c("centroid"))
set.seed(1)
permutest(result)

####
#Subset to only IAA plus control:
DNA.rfy.ra.IAA = subset_samples(DNA.rfy.ra, Treatment == "IAA" | Treatment == "Control")
sample_data(DNA.rfy.ra.IAA)


#braycurtis
bc <- phyloseq::distance(DNA.rfy.ra.IAA, method="bray")

#Permanova test using ADONIS function
# make a data frame from the sample_data
DNA.rfy.ra.df <- data.frame(sample_data(DNA.rfy.ra.IAA))
DNA.rfy.ra.df
# ADONIS test
set.seed(1)
adonis(bc ~ Treatment, data = DNA.rfy.ra.df)

#Homogeneity of dispersions test
set.seed(1)
result <- betadisper(bc, DNA.rfy.ra.df$Treatment, type=c("centroid"))
set.seed(1)
permutest(result)

####
#Subset to only JA plus control:
DNA.rfy.ra.JA = subset_samples(DNA.rfy.ra, Treatment == "JA" | Treatment == "Control")
sample_data(DNA.rfy.ra.JA)


#braycurtis
bc <- phyloseq::distance(DNA.rfy.ra.JA, method="bray")

#Permanova test using ADONIS function
# make a data frame from the sample_data
DNA.rfy.ra.df <- data.frame(sample_data(DNA.rfy.ra.JA))
DNA.rfy.ra.df
# ADONIS test
set.seed(1)
adonis(bc ~ Treatment, data = DNA.rfy.ra.df)

#Homogeneity of dispersions test
set.seed(1)
result <- betadisper(bc, DNA.rfy.ra.df$Treatment, type=c("centroid"))
set.seed(1)
permutest(result)

####
#Subset to only SA plus control:
DNA.rfy.ra.SA = subset_samples(DNA.rfy.ra, Treatment == "SA" | Treatment == "Control")
sample_data(DNA.rfy.ra.SA)

#braycurtis
bc <- phyloseq::distance(DNA.rfy.ra.SA, method="bray")

#Permanova test using ADONIS function
# make a data frame from the sample_data
DNA.rfy.ra.df <- data.frame(sample_data(DNA.rfy.ra.SA))
DNA.rfy.ra.df
# ADONIS test
set.seed(1)
adonis(bc ~ Treatment, data = DNA.rfy.ra.df)

#Homogeneity of dispersions test
set.seed(1)
result <- betadisper(bc, DNA.rfy.ra.df$Treatment, type=c("centroid"))
set.seed(1)
permutest(result)









# Subsetting for only "Active" OTUs and plotting ordinations -------------------

#Get the OTU table from the rarefied DNA&RNA dataset.
OTUs <- otu_table(DNARNA.rfy)
OTUs <- as.data.frame(OTUs)

#Make OTU table containing only 16sRNA
OTUsRNA <- OTUs[,grep('RNA',names(OTUs))]
#Make OTU table containing only 16sDNA
OTUsDNA <- OTUs[,grep('DNA',names(OTUs))]

#verify that these tables only contain the correct samples (and all have the same number of reads since rarefied)
colSums(OTUsRNA)
colSums(OTUsDNA)

#add 1 to each DNA OTU that is currently a zero, so no 0's in denominator
OTUsDNA<-replace(OTUsDNA, OTUsDNA == 0, 1)
#Make table of 16s ratios (RNA/DNA).
OTUs_16sRatio <- data.frame(OTUsRNA/OTUsDNA)

#Count number of OTUs with 16sRatio greater than a given threshold.
#Assign 'active' as 1, and 'inactive' as 0.
OTUs_16sRatio[OTUs_16sRatio<=1] = 0
OTUs_16sRatio[OTUs_16sRatio>1] = 1

#Next, subset the DNA/RNA rarefied dataset for ONLY RNA samples.
DNARNA.rfy.OTUs <- otu_table(DNARNA.rfy)
OTUs_RNA_DNARNA.rfy <- DNARNA.rfy.OTUs[,grep('RNA',names(OTUs))]

#Make that subset into dataframe
OTUs_RNA_DNARNA.rfy <-as(otu_table(OTUs_RNA_DNARNA.rfy), "matrix")
OTUs_RNA_DNARNA.rfy <- as.data.frame(OTUs_RNA_DNARNA.rfy)

#Next, multiply the rarefied RNA OTU table by the 1vs0 ActiveVsInactive table.
#This will return the rarefied#reads for all active taxa, and change to 0 all inactive taxa.
ActiveRNA <- data.frame(OTUs_RNA_DNARNA.rfy*OTUs_16sRatio)

OTUs_RNA_DNARNA.rfy["OTU8","A1RNA"]
OTUs_16sRatio["OTU8","A1RNA"]
ActiveRNA["OTU8","A1RNA"]

#Convert to matrix, then to phyloseq format.
ActiveRNA <- data.matrix(ActiveRNA)
ActiveRNA = otu_table(ActiveRNA, taxa_are_rows = TRUE)

#Create phyloseq object containing active taxa only
physeqActiveRNA.rfy = phyloseq(ActiveRNA, TAX, METADATA)
physeqActiveRNA.rfy

#Convert Rarified data to Relative Abundance Data
physeqActiveRNA.rfy.ra = transform_sample_counts(physeqActiveRNA.rfy, function(x) x / sum(x) )
#Double-check it worked: sample sums should equal 1 
sample_sums(otu_table(physeqActiveRNA.rfy.ra)) 


#braycurtis
bc <- phyloseq::distance(physeqActiveRNA.rfy.ra, method="bray")

#pcoa plot
set.seed(1)
# Ordinate
RNA_ACTIVE_pcoa <- ordinate(physeq = physeqActiveRNA.rfy.ra, method = "PCoA", distance = "bc")
# Plot 
RNA_ACTIVE_pcoa <- plot_ordination(physeq = physeqActiveRNA.rfy.ra,ordination = RNA_ACTIVE_pcoa,color = "Treatment",title = "Treatment effects on Active rhizosphere soil communities") + geom_point(size = 3)  
#Change size of graph
RNA_ACTIVE_pcoa + coord_fixed()+ scale_color_manual(name="Treatment", breaks=c("Pre-dry", "Post-dry", "Post-water", "Control", "ABA", "IAA", "JA", "SA"), values=c("red","blue","yellow","orange","magenta4","hotpink","chocolate4","green4"))
ggsave("8trts_RNA_Active_pcoa.pdf",width=178,units="mm") 

#Permanova test using ADONIS function
# make a data frame from the sample_data
physeqActiveRNA.rfy.ra.df <- data.frame(sample_data(physeqActiveRNA.rfy.ra))
physeqActiveRNA.rfy.ra.df
# ADONIS test
set.seed(1)
adonis(bc ~ Treatment, data = physeqActiveRNA.rfy.ra.df)

#Homogeneity of dispersions test
set.seed(1)
result <- betadisper(bc, physeqActiveRNA.rfy.ra.df$Treatment, type=c("centroid"))
set.seed(1)
permutest(result)


####
#Subset to only ABA plus control:
physeqActiveRNA.rfy.ra.ABA = subset_samples(physeqActiveRNA.rfy.ra, Treatment == "ABA" | Treatment == "Control")
sample_data(physeqActiveRNA.rfy.ra.ABA)



#braycurtis
bc <- phyloseq::distance(physeqActiveRNA.rfy.ra.ABA, method="bray")

#Permanova test using ADONIS function
# make a data frame from the sample_data
physeqActiveRNA.rfy.ra.df <- data.frame(sample_data(physeqActiveRNA.rfy.ra.ABA))
physeqActiveRNA.rfy.ra.df
# ADONIS test
set.seed(1)
adonis(bc ~ Treatment, data = physeqActiveRNA.rfy.ra.df)

#Homogeneity of dispersions test
set.seed(1)
result <- betadisper(bc, physeqActiveRNA.rfy.ra.df$Treatment, type=c("centroid"))
set.seed(1)
permutest(result)

####
#Subset to only IAA plus control:
physeqActiveRNA.rfy.ra.IAA = subset_samples(physeqActiveRNA.rfy.ra, Treatment == "IAA" | Treatment == "Control")
sample_data(physeqActiveRNA.rfy.ra.IAA)



#braycurtis
bc <- phyloseq::distance(physeqActiveRNA.rfy.ra.IAA, method="bray")

#Permanova test using ADONIS function
# make a data frame from the sample_data
physeqActiveRNA.rfy.ra.df <- data.frame(sample_data(physeqActiveRNA.rfy.ra.IAA))
physeqActiveRNA.rfy.ra.df
# ADONIS test
set.seed(1)
adonis(bc ~ Treatment, data = physeqActiveRNA.rfy.ra.df)

#Homogeneity of dispersions test
set.seed(1)
result <- betadisper(bc, physeqActiveRNA.rfy.ra.df$Treatment, type=c("centroid"))
set.seed(1)
permutest(result)

####
#Subset to only JA plus control:
physeqActiveRNA.rfy.ra.JA = subset_samples(physeqActiveRNA.rfy.ra, Treatment == "JA" | Treatment == "Control")
sample_data(physeqActiveRNA.rfy.ra.JA)


#braycurtis
bc <- phyloseq::distance(physeqActiveRNA.rfy.ra.JA, method="bray")

#Permanova test using ADONIS function
# make a data frame from the sample_data
physeqActiveRNA.rfy.ra.df <- data.frame(sample_data(physeqActiveRNA.rfy.ra.JA))
physeqActiveRNA.rfy.ra.df
# ADONIS test
set.seed(1)
adonis(bc ~ Treatment, data = physeqActiveRNA.rfy.ra.df)

#Homogeneity of dispersions test
set.seed(1)
result <- betadisper(bc, physeqActiveRNA.rfy.ra.df$Treatment, type=c("centroid"))
set.seed(1)
permutest(result)

####
#Subset to only SA plus control:
physeqActiveRNA.rfy.ra.SA = subset_samples(physeqActiveRNA.rfy.ra, Treatment == "SA" | Treatment == "Control")
sample_data(physeqActiveRNA.rfy.ra.SA)



#braycurtis
bc <- phyloseq::distance(physeqActiveRNA.rfy.ra.SA, method="bray")

#Permanova test using ADONIS function
# make a data frame from the sample_data
physeqActiveRNA.rfy.ra.df <- data.frame(sample_data(physeqActiveRNA.rfy.ra.SA))
physeqActiveRNA.rfy.ra.df
# ADONIS test
set.seed(1)
adonis(bc ~ Treatment, data = physeqActiveRNA.rfy.ra.df)

#Homogeneity of dispersions test
set.seed(1)
result <- betadisper(bc, physeqActiveRNA.rfy.ra.df$Treatment, type=c("centroid"))
set.seed(1)
permutest(result)









#Subset to only the four hormones plus control:
physeqActiveRNA.rfy.ra = subset_samples(physeqActiveRNA.rfy.ra, Treatment != "Pre-dry" & Treatment != "Post-dry" & Treatment != "Post-water")
physeqActiveRNA.rfy.ra
sample_data(physeqActiveRNA.rfy.ra)


#braycurtis
bc <- phyloseq::distance(physeqActiveRNA.rfy.ra, method="bray")

#pcoa plot: normalized weighted Unifrac
set.seed(1)
# Ordinate
RNA_pcoa <- ordinate(physeq = physeqActiveRNA.rfy.ra, method = "PCoA", distance = "bc")
# Plot 
RNA_pcoa <- plot_ordination(physeq = physeqActiveRNA.rfy.ra,ordination = RNA_pcoa,color = "Treatment",title = "Treatment effects on Active rhizosphere soil communities") + geom_point(size = 3)  
#Change size of graph
RNA_pcoa + coord_fixed()+
  scale_color_manual(name="Treatment", breaks=c("Control", "ABA", "IAA", "JA", "SA"), values=c("red","blue","yellow","orange","green4"))
ggsave("5trts_RNA_Active_pcoa.pdf",width=178,units="mm") 

#Permanova test using ADONIS function
# make a data frame from the sample_data
physeqActiveRNA.rfy.ra.df <- data.frame(sample_data(physeqActiveRNA.rfy.ra))
physeqActiveRNA.rfy.ra.df
# ADONIS test
set.seed(1)
adonis(bc ~ Treatment, data = physeqActiveRNA.rfy.ra.df)

#Homogeneity of dispersions test
set.seed(1)
result <- betadisper(bc, physeqActiveRNA.rfy.ra.df$Treatment, type=c("centroid"))
set.seed(1)
permutest(result)
TukeyHSD(result, which = "group", ordered = FALSE,conf.level = 0.95)









# T-tests comparing OTU RNA/DNA ratios across TREATMENTS ----------------------------------------


#Get the OTU table from the rarefied DNA&RNA dataset.
OTUs <- otu_table(DNARNA.rfy)
OTUs <- as.data.frame(OTUs)

#Make OTU table containing only 16sRNA
OTUsRNA <- OTUs[,grep('RNA',names(OTUs))]
#Make OTU table containing only 16sDNA
OTUsDNA <- OTUs[,grep('DNA',names(OTUs))]

#verify that these tables only contain the correct samples (and all have the same number of reads since rarefied)
colSums(OTUsRNA)
colSums(OTUsDNA)

#add 1 to each DNA OTU that is currently a zero, so no 0's in denominator
OTUsDNA<-replace(OTUsDNA, OTUsDNA == 0, 1)
#Make table of 16s ratios (RNA/DNA).
OTUs_16sRatio <- data.frame(OTUsRNA/OTUsDNA)


#Prepare input files for t-tests.
DNARNA.rfy.map=(sample_data(DNARNA.rfy))
#retain only the RNA rows.
DNARNA.rfy.map <- DNARNA.rfy.map[grep('RNA',row.names(DNARNA.rfy.map)),]
#assign which treatment belongs to which samples
trt=t(DNARNA.rfy.map[,"Treatment"])
trt

###ABA vs Control
#Prepare output file into which the t-test results will go.
otuname=row.names(OTUs_16sRatio)
ttest.out.ABA <- data.frame(matrix(nrow = nrow(OTUs_16sRatio), ncol = 4))
colnames(ttest.out.ABA) <- c("OTU", "t-statistic","df", "p.value")
rownames(ttest.out.ABA) <- row.names(OTUs_16sRatio)
#Perform the t-tests.
for(i in 1:length(otuname)){
  #subset the data to test one OTU at a time
  Hormone=OTUs_16sRatio[row.names(OTUs_16sRatio)==otuname[i],trt=="ABA"]
  Control=OTUs_16sRatio[row.names(OTUs_16sRatio)==otuname[i],trt=="Control"]
  #perform the test
  test=t.test(Hormone, Control, paired=FALSE, var.equal = FALSE)
  ttest.out.ABA[i,1] <- row.names(OTUs_16sRatio)[i]
  ttest.out.ABA[i,2] <-(test$statistic)
  ttest.out.ABA[i,3] <-(test$parameter)
  ttest.out.ABA[i,4] <-(test$p.value)
}
#Calculate adjusted p-values and add to table. Sort by adjusted p-values.
ttest.out.ABA[,5] <-p.adjust(ttest.out.ABA$p.value, method = "fdr")
colnames(ttest.out.ABA) <- c("OTU", "t-statistic","df", "p.value", "adj.p.value")
ttest.out.ABA <- ttest.out.ABA[order(ttest.out.ABA$adj.p.value),]
ttest.out.ABA
SigDiffOTUs.ABA <- row.names(subset(ttest.out.ABA, adj.p.value<=0.05))

###IAA vs Control
#Prepare output file into which the t-test results will go.
otuname=row.names(OTUs_16sRatio)
ttest.out.IAA <- data.frame(matrix(nrow = nrow(OTUs_16sRatio), ncol = 4))
colnames(ttest.out.IAA) <- c("OTU", "t-statistic","df", "p.value")
rownames(ttest.out.IAA) <- row.names(OTUs_16sRatio)
#Perform the t-tests.
for(i in 1:length(otuname)){
  #subset the data to test one OTU at a time
  Hormone=OTUs_16sRatio[row.names(OTUs_16sRatio)==otuname[i],trt=="IAA"]
  Control=OTUs_16sRatio[row.names(OTUs_16sRatio)==otuname[i],trt=="Control"]
  #perform the test
  test=t.test(Hormone, Control, paired=FALSE, var.equal = FALSE)
  ttest.out.IAA[i,1] <- row.names(OTUs_16sRatio)[i]
  ttest.out.IAA[i,2] <-(test$statistic)
  ttest.out.IAA[i,3] <-(test$parameter)
  ttest.out.IAA[i,4] <-(test$p.value)
}
#Calculate adjusted p-values and add to table. Sort by adjusted p-values.
ttest.out.IAA[,5] <-p.adjust(ttest.out.IAA$p.value, method = "fdr")
colnames(ttest.out.IAA) <- c("OTU", "t-statistic","df", "p.value", "adj.p.value")
ttest.out.IAA <- ttest.out.IAA[order(ttest.out.IAA$adj.p.value),]
ttest.out.IAA
SigDiffOTUs.IAA <- row.names(subset(ttest.out.IAA, adj.p.value<=0.05))


###JA vs Control
#Prepare output file into which the t-test results will go.
otuname=row.names(OTUs_16sRatio)
ttest.out.JA <- data.frame(matrix(nrow = nrow(OTUs_16sRatio), ncol = 4))
colnames(ttest.out.JA) <- c("OTU", "t-statistic","df", "p.value")
rownames(ttest.out.JA) <- row.names(OTUs_16sRatio)
#Perform the t-tests.
for(i in 1:length(otuname)){
  #subset the data to test one OTU at a time
  Hormone=OTUs_16sRatio[row.names(OTUs_16sRatio)==otuname[i],trt=="JA"]
  Control=OTUs_16sRatio[row.names(OTUs_16sRatio)==otuname[i],trt=="Control"]
  #perform the test
  test=t.test(Hormone, Control, paired=FALSE, var.equal = FALSE)
  ttest.out.JA[i,1] <- row.names(OTUs_16sRatio)[i]
  ttest.out.JA[i,2] <-(test$statistic)
  ttest.out.JA[i,3] <-(test$parameter)
  ttest.out.JA[i,4] <-(test$p.value)
}
#Calculate adjusted p-values and add to table. Sort by adjusted p-values.
ttest.out.JA[,5] <-p.adjust(ttest.out.JA$p.value, method = "fdr")
colnames(ttest.out.JA) <- c("OTU", "t-statistic","df", "p.value", "adj.p.value")
ttest.out.JA <- ttest.out.JA[order(ttest.out.JA$adj.p.value),]
ttest.out.JA
SigDiffOTUs.JA <- row.names(subset(ttest.out.JA, adj.p.value<=0.05))
SigDiffOTUs.JA

###SA vs Control
#Prepare output file into which the t-test results will go.
otuname=row.names(OTUs_16sRatio)
ttest.out.SA <- data.frame(matrix(nrow = nrow(OTUs_16sRatio), ncol = 4))
colnames(ttest.out.SA) <- c("OTU", "t-statistic","df", "p.value")
rownames(ttest.out.SA) <- row.names(OTUs_16sRatio)
#Perform the t-tests.
for(i in 1:length(otuname)){
  #subset the data to test one OTU at a time
  Hormone=OTUs_16sRatio[row.names(OTUs_16sRatio)==otuname[i],trt=="SA"]
  Control=OTUs_16sRatio[row.names(OTUs_16sRatio)==otuname[i],trt=="Control"]
  #perform the test
  test=t.test(Hormone, Control, paired=FALSE, var.equal = FALSE)
  ttest.out.SA[i,1] <- row.names(OTUs_16sRatio)[i]
  ttest.out.SA[i,2] <-(test$statistic)
  ttest.out.SA[i,3] <-(test$parameter)
  ttest.out.SA[i,4] <-(test$p.value)
}
#Calculate adjusted p-values and add to table. Sort by adjusted p-values.
ttest.out.SA[,5] <-p.adjust(ttest.out.SA$p.value, method = "fdr")
colnames(ttest.out.SA) <- c("OTU", "t-statistic","df", "p.value", "adj.p.value")
ttest.out.SA <- ttest.out.SA[order(ttest.out.SA$adj.p.value),]
ttest.out.SA
SigDiffOTUs.SA <- row.names(subset(ttest.out.SA, adj.p.value<=0.05))









# T-tests comparing Class RNA/DNA ratios across TREATMENTS ----------------


RNA.rfy.dat <- otu_table(RNA.rfy)
RNA.rfy.tax.dat <- tax_table(RNA.rfy)

#Here, choose the TAXONOMIC HEIRARCHY AT WHICH YOU WANT TO DO THE TTESTS.
RNA.rfy.tax.data <- subset(RNA.rfy.tax.dat, select = c("Family"))
RNA.rfy.tax.data
#Now build the data.frame holding relative abundance and taxonomy data.
RNA.tax.data <- merge(RNA.rfy.dat, RNA.rfy.tax.data, by="row.names")
row.names(RNA.tax.data) <- RNA.tax.data$Row.names
RNA.tax.data <- subset(RNA.tax.data, select = -c(Row.names))

#Sum within each row by the Taxonomy of Interest.
library("dplyr")
RNA.tax.summary <-RNA.tax.data %>% group_by(Family) %>% summarise_all(funs(sum))
RNA.tax.summary <- as.data.frame(RNA.tax.summary)
row.names(RNA.tax.summary) <- RNA.tax.summary$Family
RNA.tax.summary <- subset(RNA.tax.summary, select = -c(Family))


###now do the same thing for DNA

DNA.rfy.dat <- otu_table(DNA.rfy)
DNA.rfy.tax.dat <- tax_table(DNA.rfy)

#Here, choose the TAXONOMIC HEIRARCHY AT WHICH YOU WANT TO DO THE TTESTS.
DNA.rfy.tax.data <- subset(DNA.rfy.tax.dat, select = c("Family"))
DNA.rfy.tax.data
#Now build the data.frame holding relative abundance and taxonomy data.
DNA.tax.data <- merge(DNA.rfy.dat, DNA.rfy.tax.data, by="row.names")
row.names(DNA.tax.data) <- DNA.tax.data$Row.names
DNA.tax.data <- subset(DNA.tax.data, select = -c(Row.names))

#Sum within each row by the Taxonomy of Interest.
library("dplyr")
DNA.tax.summary <-DNA.tax.data %>% group_by(Family) %>% summarise_all(funs(sum))
DNA.tax.summary <- as.data.frame(DNA.tax.summary)
row.names(DNA.tax.summary) <- DNA.tax.summary$Family
DNA.tax.summary <- subset(DNA.tax.summary, select = -c(Family))



#verify that these tables only contain the correct samples (and all have the same number of reads since rarefied)
colSums(RNA.tax.summary)
colSums(DNA.tax.summary)

#add 1 to each DNA OTU that is currently a zero, so no 0's in denominator
DNA.tax.summary<-replace(DNA.tax.summary, DNA.tax.summary == 0, 1)
#Make table of 16s ratios (RNA/DNA).
Class_16sRatio <- data.frame(RNA.tax.summary/DNA.tax.summary)


#Prepare input files for t-tests.
RNA.rfy.map=(sample_data(RNA.rfy))
#assign which treatment belongs to which samples
trt=t(RNA.rfy.map[,"Treatment"])
trt

###ABA vs Control
#Prepare output file into which the t-test results will go.
classname=row.names(Class_16sRatio)
ttest.out.ABA <- data.frame(matrix(nrow = nrow(Class_16sRatio), ncol = 4))
colnames(ttest.out.ABA) <- c("Class", "t-statistic","df", "p.value")
rownames(ttest.out.ABA) <- row.names(Class_16sRatio)
#Perform the t-tests.
for(i in 1:length(classname)){
  #subset the data to test one OTU at a time
  Hormone=Class_16sRatio[row.names(Class_16sRatio)==classname[i],trt=="ABA"]
  Control=Class_16sRatio[row.names(Class_16sRatio)==classname[i],trt=="Control"]
  #perform the test
  test=t.test(Hormone, Control, paired=FALSE, var.equal = FALSE)
  ttest.out.ABA[i,1] <- row.names(Class_16sRatio)[i]
  ttest.out.ABA[i,2] <-(test$statistic)
  ttest.out.ABA[i,3] <-(test$parameter)
  ttest.out.ABA[i,4] <-(test$p.value)
}
#Calculate adjusted p-values and add to table. Sort by adjusted p-values.
ttest.out.ABA[,5] <-p.adjust(ttest.out.ABA$p.value, method = "fdr")
colnames(ttest.out.ABA) <- c("Class", "t-statistic","df", "p.value", "adj.p.value")
ttest.out.ABA <- ttest.out.ABA[order(ttest.out.ABA$adj.p.value),]
ttest.out.ABA
SigDiffClass.ABA <- row.names(subset(ttest.out.ABA, adj.p.value<=0.05))

###IAA vs Control
#Prepare output file into which the t-test results will go.
classname=row.names(Class_16sRatio)
ttest.out.IAA <- data.frame(matrix(nrow = nrow(Class_16sRatio), ncol = 4))
colnames(ttest.out.IAA) <- c("Class", "t-statistic","df", "p.value")
rownames(ttest.out.IAA) <- row.names(Class_16sRatio)
#Perform the t-tests.
for(i in 1:length(classname)){
  #subset the data to test one OTU at a time
  Hormone=Class_16sRatio[row.names(Class_16sRatio)==classname[i],trt=="IAA"]
  Control=Class_16sRatio[row.names(Class_16sRatio)==classname[i],trt=="Control"]
  #perform the test
  test=t.test(Hormone, Control, paired=FALSE, var.equal = FALSE)
  ttest.out.IAA[i,1] <- row.names(Class_16sRatio)[i]
  ttest.out.IAA[i,2] <-(test$statistic)
  ttest.out.IAA[i,3] <-(test$parameter)
  ttest.out.IAA[i,4] <-(test$p.value)
}
#Calculate adjusted p-values and add to table. Sort by adjusted p-values.
ttest.out.IAA[,5] <-p.adjust(ttest.out.IAA$p.value, method = "fdr")
colnames(ttest.out.IAA) <- c("Class", "t-statistic","df", "p.value", "adj.p.value")
ttest.out.IAA <- ttest.out.IAA[order(ttest.out.IAA$adj.p.value),]
ttest.out.IAA
SigDiffClass.IAA <- row.names(subset(ttest.out.IAA, adj.p.value<=0.05))


###JA vs Control
#Prepare output file into which the t-test results will go.
classname=row.names(Class_16sRatio)
ttest.out.JA <- data.frame(matrix(nrow = nrow(Class_16sRatio), ncol = 4))
colnames(ttest.out.JA) <- c("Class", "t-statistic","df", "p.value")
rownames(ttest.out.JA) <- row.names(Class_16sRatio)
#Perform the t-tests.
for(i in 1:length(classname)){
  #subset the data to test one OTU at a time
  Hormone=Class_16sRatio[row.names(Class_16sRatio)==classname[i],trt=="JA"]
  Control=Class_16sRatio[row.names(Class_16sRatio)==classname[i],trt=="Control"]
  #perform the test
  test=t.test(Hormone, Control, paired=FALSE, var.equal = FALSE)
  ttest.out.JA[i,1] <- row.names(Class_16sRatio)[i]
  ttest.out.JA[i,2] <-(test$statistic)
  ttest.out.JA[i,3] <-(test$parameter)
  ttest.out.JA[i,4] <-(test$p.value)
}
#Calculate adjusted p-values and add to table. Sort by adjusted p-values.
ttest.out.JA[,5] <-p.adjust(ttest.out.JA$p.value, method = "fdr")
colnames(ttest.out.JA) <- c("Class", "t-statistic","df", "p.value", "adj.p.value")
ttest.out.JA <- ttest.out.JA[order(ttest.out.JA$adj.p.value),]
ttest.out.JA
SigDiffClass.JA <- row.names(subset(ttest.out.JA, adj.p.value<=0.05))
SigDiffClass.JA

###SA vs Control
#Prepare output file into which the t-test results will go.
classname=row.names(Class_16sRatio)
ttest.out.SA <- data.frame(matrix(nrow = nrow(Class_16sRatio), ncol = 4))
colnames(ttest.out.SA) <- c("Class", "t-statistic","df", "p.value")
rownames(ttest.out.SA) <- row.names(Class_16sRatio)
#Perform the t-tests.
for(i in 1:length(classname)){
  #subset the data to test one OTU at a time
  Hormone=Class_16sRatio[row.names(Class_16sRatio)==classname[i],trt=="SA"]
  Control=Class_16sRatio[row.names(Class_16sRatio)==classname[i],trt=="Control"]
  #perform the test
  test=t.test(Hormone, Control, paired=FALSE, var.equal = FALSE)
  ttest.out.SA[i,1] <- row.names(Class_16sRatio)[i]
  ttest.out.SA[i,2] <-(test$statistic)
  ttest.out.SA[i,3] <-(test$parameter)
  ttest.out.SA[i,4] <-(test$p.value)
}
#Calculate adjusted p-values and add to table. Sort by adjusted p-values.
ttest.out.SA[,5] <-p.adjust(ttest.out.SA$p.value, method = "fdr")
colnames(ttest.out.SA) <- c("Class", "t-statistic","df", "p.value", "adj.p.value")
ttest.out.SA <- ttest.out.SA[order(ttest.out.SA$adj.p.value),]
ttest.out.SA
SigDiffClass.SA <- row.names(subset(ttest.out.SA, adj.p.value<=0.05))

SigDiffClass.SA


# Rank Abundance Curves for RNA versus DNA --------------------------------

### Phantom Taxa
Corn.rfy <- otu_table(DNARNA.rfy)
Corn.rfy <- as.data.frame(Corn.rfy)
#Make OTU table containing only 16sRNA
Corn_OTUsRNA <- Corn.rfy[,grep('RNA',names(Corn.rfy))]
#Make OTU table containing only 16sDNA
Corn_OTUsDNA <- Corn.rfy[,grep('DNA',names(Corn.rfy))]

#in RNA table, replace RNA>=1 as 1, and keep all zeros as zero.
Corn_OTUsRNA[Corn_OTUsRNA>=1] = 1
Corn_OTUsRNA[Corn_OTUsRNA<1] = 0
#In DNA table, replace all 0's with 1, and everything else with 0.
#First, make things that are currently >=1 into a higher number.
Corn_OTUsDNA[Corn_OTUsDNA>=1] = 2
#then, make everything that's currently 0 into a 1.
Corn_OTUsDNA[Corn_OTUsDNA<1] = 1
#then, make everything greater than 1 into 0.
Corn_OTUsDNA[Corn_OTUsDNA>1] = 0
#sum the two tables
Corn_RNAnoDNA <- data.frame(Corn_OTUsRNA+Corn_OTUsDNA)

#NOW, everthing with a 2 is a phantom taxa.
#replace everything with a 2 with "Phantom"...all else with "NotPhantom".
Corn_RNAnoDNA[Corn_RNAnoDNA=="2"] <- "Phantom"
Corn_RNAnoDNA[Corn_RNAnoDNA=="1"] <- "NotPhantom"
Corn_RNAnoDNA[Corn_RNAnoDNA=="0"] <- "NotPhantom"



#make figure: rank abundance curve for RNA

#make corn physeq object containing only RNA
Corn_OTUs <- otu_table(DNARNA.rfy)
Corn_OTUs <- as.data.frame(Corn_OTUs)
#Make OTU table containing only 16sRNA
Corn_OTUsRNA <- Corn_OTUs[,grep('RNA',names(Corn_OTUs))]

#convert to matrix, then phyloseq format
Corn_OTUsRNA <- data.matrix(Corn_OTUsRNA)
Corn_OTUsRNA = otu_table(Corn_OTUsRNA, taxa_are_rows = TRUE)
Corn_rare_RNA = phyloseq(Corn_OTUsRNA, TAX, METADATA)

# this converts taxa counts in each sample to a percentage
phyloTemp_CornRNA = transform_sample_counts(Corn_rare_RNA, function(x) 1e+02 * x/sum(x))
clusterData_CornRNA = psmelt(phyloTemp_CornRNA)
clusterData_CornRNA = filter(clusterData_CornRNA,Abundance > 0)
# this is where the mean is calculated and the taxa to display is chosen
clusterAgg_CornRNA = aggregate(Abundance ~ OTU,data=clusterData_CornRNA,mean)
# filtering and picking the number to display
clusterAgg_CornRNA = clusterAgg_CornRNA[order(-clusterAgg_CornRNA$Abundance),][1:100,]


#combine with data on Phantom taxa if you want to label them in the Figure
#ClusterAgg_Fig <- merge(clusterAgg_CornRNA,Corn_RNAnoDNA,by.x="OTU",by.y="row.names",all=FALSE)
#rename columns
#names(ClusterAgg_Fig) <- c("OTU", "Abundance","Phantom")

Corn_RNAFig<-ggplot(clusterAgg_CornRNA,aes(x=reorder(OTU,-Abundance),y=Abundance)) +
  geom_point() + 
  theme(axis.ticks = element_blank(), axis.text.x = element_blank()) +
  ylab("Relative Abundance")+
  xlab("OTU Rank")+
  ggtitle("Corn RNA Dataset")+
  coord_cartesian(ylim=c(0,15))
Corn_RNAFig



#make figure: rank abundance curve for DNA

#make corn physeq object containing only DNA
Corn_OTUs <- otu_table(DNARNA.rfy)
Corn_OTUs <- as.data.frame(Corn_OTUs)
#Make OTU table containing only 16sRNA
Corn_OTUsDNA <- Corn_OTUs[,grep('DNA',names(Corn_OTUs))]
#convert to matrix, then phyloseq format
Corn_OTUsDNA <- data.matrix(Corn_OTUsDNA)
Corn_OTUsDNA = otu_table(Corn_OTUsDNA, taxa_are_rows = TRUE)
Corn_rare_DNA = phyloseq(Corn_OTUsDNA, TAX, METADATA)

# this converts taxa counts in each sample to a percentage
phyloTemp_CornDNA = transform_sample_counts(Corn_rare_DNA, function(x) 1e+02 * x/sum(x))
clusterData_CornDNA = psmelt(phyloTemp_CornDNA)
clusterData_CornDNA = filter(clusterData_CornDNA,Abundance > 0)
# this is where the mean is calculated and the taxa to display is chosen
clusterAgg_CornDNA = aggregate(Abundance ~ OTU,data=clusterData_CornDNA,mean)
# filtering and picking the number to display
clusterAgg_CornDNA = clusterAgg_CornDNA[order(-clusterAgg_CornDNA$Abundance),][1:100,]

Corn_DNAFig<-ggplot(clusterAgg_CornDNA,aes(x=reorder(OTU,-Abundance),y=Abundance)) +
  geom_point() + 
  theme(axis.ticks = element_blank(), axis.text.x = element_blank()) +
  ylab("Relative Abundance")+
  xlab("OTU Rank")+
  ggtitle("Corn DNA Dataset")+
  coord_cartesian(ylim=c(0,15))
Corn_DNAFig


