library("phyloseq")
packageVersion("phyloseq")
library("ggplot2")
packageVersion("ggplot2")
theme_set(theme_bw())
library("reshape2")
packageVersion("reshape2")
set.seed(1)

#Import OTU table ('#' in cell A1 has been removed)
Bean_OTU_table <- read.table("Clean_OTU_table.txt", sep="\t", row.names=1, header=TRUE)
#Convert to matrix, then to phyloseq format.
Bean_OTU_table <- data.matrix(Bean_OTU_table)
Bean_OTU_table = otu_table(Bean_OTU_table, taxa_are_rows = TRUE)
#Import taxonomy (has been formatted in Excel to be tab-delimited)
Bean_TAX <- read.table("Clean_otus_tax_assignments.txt", sep="\t", row.names=1, stringsAsFactors = FALSE, header=FALSE, na.strings="")
colnames(Bean_TAX) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
#convert to matrix, then phyloseq format.
Bean_TAX <- as.matrix(Bean_TAX)
Bean_TAX = tax_table(Bean_TAX)
#Import sample data
Bean_METADATA <- read.csv("bean_seed_rhizo_map.txt", sep="\t", row.names=1, header=TRUE, na.strings="")
#convert to phyloseq format.
Bean_METADATA<- sample_data(Bean_METADATA)
#Create phyloseq object:
Bean_physeq = phyloseq(Bean_OTU_table, Bean_TAX, Bean_METADATA)
Bean_physeq
#Check for singletons and remove if any.
any(taxa_sums(Bean_physeq) <= 1) #FALSE
#Remove taxa with Domain designated as "Eukaryota" (Chloroplast and mitochondria have already been removed on HPCC).
Bean_physeq <- subset_taxa(Bean_physeq, Domain!= "Eukaryota") #8148 taxa remaining
Bean_physeq
sample_sums(Bean_physeq)
sorted <- sort(sample_sums(Bean_physeq))
sorted <- as.data.frame(sorted)
sorted
#Actually do the rarefying here.
set.seed(1)
Bean.rfy <- rarefy_even_depth(Bean_physeq, sample.size = 37815, replace = FALSE)
#52 OTUs were removed: no longer present in any sample after subsampling.

#Convert Rarified data to Relative Abundance Data
bean.rfy.ra = transform_sample_counts(Bean.rfy, function(x) x / sum(x) )
#Double-check it worked: sample sums should equal 1 
sample_sums(otu_table(bean.rfy.ra))




# Ordination and PERMANOVA for RNA and DNA --------------------------------

#braycurtis
set.seed(1)
bc  <-distance(bean.rfy.ra, method="bray")
#pcoa plot
# Ordinate
pcoa.bc <- ordinate(physeq = bean.rfy.ra, method = "PCoA", distance = "bc")
# Plot
pcoa_bc <- plot_ordination(physeq = bean.rfy.ra,ordination = pcoa.bc,shape = "Type",color = "Treatment") + geom_point(size = 3)+ scale_colour_brewer(type="qual", palette="Set1") 
#Change size of graph
pcoa_bc + 
  coord_fixed()+ 
  scale_color_manual(name="Treatment", breaks=c("Control", "Drought", "Nutrients"), values=c("blue","red","green4"))+
  scale_shape_discrete(name="Type", breaks=c("DNA", "RNA"), labels=c("16s rRNA Gene","16s rRNA"))+
  xlab("PCoA1: 25.1% var. explained")+
  ylab("PCoA2: 18.2% var. explained")+
  #The next line moves the "type" legend above the "treatment" legend.
  guides(colour = guide_legend(order = 2), 
         shape = guide_legend(order = 1))

ggsave("pcoa_all_bc.pdf",width=178,units="mm") 

#Permanova test using ADONIS function
# make a data frame from the sample_data
physeq.rfy.ra.df <- data.frame(sample_data(bean.rfy.ra))
physeq.rfy.ra.df
# ADONIS test
library("vegan")
set.seed(1)
adonis(bc ~ Type*Treatment, data = physeq.rfy.ra.df)
set.seed(1)
adonis(bc ~ Treatment*Type, data = physeq.rfy.ra.df)

#Homogeneity of dispersions test
set.seed(1)
result <- betadisper(bc, physeq.rfy.ra.df$Type, type=c("centroid"))
set.seed(1)
permutest(result)
set.seed(1)
result <- betadisper(bc, physeq.rfy.ra.df$Treatment, type=c("centroid"))
set.seed(1)
permutest(result)




# Ordination and PERMANOVA for DNA ----------------------------------------

#Subset the  relative abundance dataset for ONLY DNA samples.
DNA.rfy.ra.OTUs <- otu_table(bean.rfy.ra)
OTUs_DNA.rfy.ra <- DNA.rfy.ra.OTUs[ , !c(TRUE,FALSE) ] 
test <-as(otu_table(OTUs_DNA.rfy.ra), "matrix")
test <- as.data.frame(test)

#Convert to matrix, then to phyloseq format.
testDNA <- data.matrix(test)
testDNA = otu_table(testDNA, taxa_are_rows = TRUE)

#Import taxonomy (has been formatted in Excel to be tab-delimited)
Bean_TAX <- read.table("Clean_otus_tax_assignments.txt", sep="\t", row.names=1, stringsAsFactors = FALSE, header=FALSE, na.strings="")
colnames(Bean_TAX) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
#convert to matrix, then phyloseq format.
Bean_TAX <- as.matrix(Bean_TAX)
Bean_TAX = tax_table(Bean_TAX)
#Import sample data
Bean_METADATA <- read.csv("bean_seed_rhizo_map.txt", sep="\t", row.names=1, header=TRUE, na.strings="")
#convert to phyloseq format.
Bean_METADATA<- sample_data(Bean_METADATA)

#Create phyloseq object containing DNA only
DNA.rfy.ra = phyloseq(testDNA, Bean_TAX, Bean_METADATA)
DNA.rfy.ra

sample_sums(DNA.rfy.ra)


#braycurtis
set.seed(1)
bc  <-distance(DNA.rfy.ra, method="bray")
#pcoa plot
# Ordinate
pcoa.bc <- ordinate(physeq = DNA.rfy.ra, method = "PCoA", distance = "bc")
# Plot
pcoa_bc <- plot_ordination(physeq = DNA.rfy.ra,ordination = pcoa.bc,color = "Treatment") + geom_point(size = 3)+ scale_colour_brewer(type="qual", palette="Set1") 
#Change size of graph
pcoa_bc + 
  coord_fixed()+ 
  scale_color_manual(name="Treatment", breaks=c("Control", "Drought", "Nutrients"), values=c("blue","red","green4"))+
  #scale_shape_discrete(name="Type", breaks=c("DNA", "RNA"), labels=c("16s rRNA Gene","16s rRNA"))+
  xlab("PCoA1: 16.2% var. explained")+
  ylab("PCoA2: 13.3% var. explained")+
  #The next line moves the "type" legend above the "treatment" legend.
  guides(colour = guide_legend(order = 2), 
         shape = guide_legend(order = 1))

ggsave("pcoa_DNA_bc.pdf",width=178,units="mm") 



#Permanova test using ADONIS function
# make a data frame from the sample_data
DNA.rfy.ra.df <- data.frame(sample_data(DNA.rfy.ra))
DNA.rfy.ra.df
# ADONIS test
library("vegan")
set.seed(1)
adonis(bc ~ Treatment, data = DNA.rfy.ra.df)

#Homogeneity of dispersions test
set.seed(1)
result <- betadisper(bc, DNA.rfy.ra.df$Treatment, type=c("centroid"))
set.seed(1)
permutest(result)


####
#Subset to only drought plus control:
DNA.rfy.ra.drctrl = subset_samples(DNA.rfy.ra, Treatment == "Drought" | Treatment == "Control")
sample_data(DNA.rfy.ra.drctrl)
#normalized Weighted Unifrac: considers relative abundance of taxa shared between samples.
bc <-distance(DNA.rfy.ra.drctrl, method="bray")

#Permanova test using ADONIS function
# make a data frame from the sample_data
DNA.rfy.ra.df <- data.frame(sample_data(DNA.rfy.ra.drctrl))
DNA.rfy.ra.df
# ADONIS test
set.seed(1)
adonis(bc ~ Treatment, data = DNA.rfy.ra.df)

#Homogeneity of dispersions test
set.seed(1)
result <- betadisper(bc, DNA.rfy.ra.df$Treatment, type=c("centroid"))
set.seed(1)
permutest(result)



#Subset to only nutrients plus control:
DNA.rfy.ra.nutctrl = subset_samples(DNA.rfy.ra, Treatment == "Nutrients" | Treatment == "Control")
sample_data(DNA.rfy.ra.nutctrl)
#normalized Weighted Unifrac: considers relative abundance of taxa shared between samples.
bc <-distance(DNA.rfy.ra.nutctrl, method="bray")

#Permanova test using ADONIS function
# make a data frame from the sample_data
DNA.rfy.ra.df <- data.frame(sample_data(DNA.rfy.ra.nutctrl))
DNA.rfy.ra.df
# ADONIS test
set.seed(1)
adonis(bc ~ Treatment, data = DNA.rfy.ra.df)

#Homogeneity of dispersions test
set.seed(1)
result <- betadisper(bc, DNA.rfy.ra.df$Treatment, type=c("centroid"))
set.seed(1)
permutest(result)




#Subset to only drought plus nutrients:
DNA.rfy.ra.nutdr = subset_samples(DNA.rfy.ra, Treatment == "Nutrients" | Treatment == "Drought")
sample_data(DNA.rfy.ra.nutdr)
#normalized Weighted Unifrac: considers relative abundance of taxa shared between samples.
bc <-distance(DNA.rfy.ra.nutdr, method="bray")

#Permanova test using ADONIS function
# make a data frame from the sample_data
DNA.rfy.ra.df <- data.frame(sample_data(DNA.rfy.ra.nutdr))
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
OTUs <- otu_table(Bean.rfy)
OTUs <- as.data.frame(OTUs)

names(OTUs) = gsub(pattern = "cDNA", replacement = "RNA", x = names(OTUs))

#Make OTU table containing only 16sRNA
OTUsRNA <- OTUs[,grep('RNA',names(OTUs))]
#Make OTU table containing only 16sDNA
OTUsDNA <- OTUs[,grep('DNA',names(OTUs))]

#verify that these tables only contain the correct samples (and all have the same number of reads since rarefied)
colSums(OTUsRNA)
colSums(OTUsDNA)


#Make table of 16s ratios (RNA/DNA) and rename column headers.
OTUs_16sRatio <- data.frame(OTUsRNA/OTUsDNA)

#Correct the values in the table as follows:
#1. When RNA=0 and DNA>0, the value is '0'. We want these cases to remain as '0'.
#2. When RNA=0 and DNA=0, the value is 'NaN'. Must change these to 'NA'.
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
OTUs_16sRatio[is.nan(OTUs_16sRatio)] <- NA
#3. When RNA>0 and DNA=0, the value is 'Inf'. Must change these to 100 (simply to indicate that we'd call these OTUs 'active', regardless of what threshold is chosen as 'active')
is.infinite.data.frame <- function(x)
  do.call(cbind, lapply(x, is.infinite))
OTUs_16sRatio[is.infinite(OTUs_16sRatio)] <- 100

#Count number of OTUs with 16sRatio greater than a given threshold.
#Assign 'active' as 1, and 'inactive' as 0.
OTUs_16sRatio[OTUs_16sRatio<=1] = 0
OTUs_16sRatio[OTUs_16sRatio>1] = 1

#Next, subset the DNA/RNA rarefied dataset for ONLY RNA samples.
Bean.rfy.OTUs <- otu_table(Bean.rfy)
Bean_RNA.rfy <- Bean.rfy.OTUs[ , !c(FALSE,TRUE) ]
#Make that subset into dataframe
OTUs_RNA_DNARNA.rfy <-as(otu_table(Bean_RNA.rfy), "matrix")
OTUs_RNA_DNARNA.rfy <- as.data.frame(OTUs_RNA_DNARNA.rfy)

#Next, multiply the rarefied RNA OTU table by the 1vs0 ActiveVsInactive table.
#This will return the rarefied#reads for all active taxa, and change to 0 all inactive taxa.
ActiveRNA <- data.frame(OTUs_RNA_DNARNA.rfy*OTUs_16sRatio)

OTUs_RNA_DNARNA.rfy["OTU256","C2cDNA"]
OTUs_16sRatio["OTU256","C2RNA"]
ActiveRNA["OTU256","C2cDNA"]

#Replace all NA's with 0.
is.na.data.frame <- function(x)
  do.call(cbind, lapply(x, is.na))
ActiveRNA[is.na(ActiveRNA)] <- 0

#Convert to matrix, then to phyloseq format.
ActiveRNA <- data.matrix(ActiveRNA)
ActiveRNA = otu_table(ActiveRNA, taxa_are_rows = TRUE)

#Import taxonomy (has been formatted in Excel to be tab-delimited)
Bean_TAX <- read.table("Clean_otus_tax_assignments.txt", sep="\t", row.names=1, stringsAsFactors = FALSE, header=FALSE, na.strings="")
colnames(Bean_TAX) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
#convert to matrix, then phyloseq format.
Bean_TAX <- as.matrix(Bean_TAX)
Bean_TAX = tax_table(Bean_TAX)
#Import sample data
Bean_METADATA <- read.csv("bean_seed_rhizo_map.txt", sep="\t", row.names=1, header=TRUE, na.strings="")
#convert to phyloseq format.
Bean_METADATA<- sample_data(Bean_METADATA)

#Create phyloseq object containing active taxa only
physeqActiveRNA.rfy = phyloseq(ActiveRNA, Bean_TAX, Bean_METADATA)
physeqActiveRNA.rfy

#Convert Rarified data to Relative Abundance Data
physeqActiveRNA.rfy.ra = transform_sample_counts(physeqActiveRNA.rfy, function(x) x / sum(x) )
#Double-check it worked: sample sums should equal 1 
sample_sums(otu_table(physeqActiveRNA.rfy.ra)) 



#braycurtis
set.seed(1)
bc  <-distance(physeqActiveRNA.rfy.ra, method="bray")
#pcoa plot
# Ordinate
pcoa.bc <- ordinate(physeq = physeqActiveRNA.rfy.ra, method = "PCoA", distance = "bc")
# Plot
pcoa_bc <- plot_ordination(physeq = physeqActiveRNA.rfy.ra,ordination = pcoa.bc,color = "Treatment") + geom_point(size = 3)+ scale_colour_brewer(type="qual", palette="Set1") 
#Change size of graph
pcoa_bc + 
  coord_fixed()+ 
  scale_color_manual(name="Treatment", breaks=c("Control", "Drought", "Nutrients"), values=c("blue","red","green4"))+
  #scale_shape_discrete(name="Type", breaks=c("DNA", "RNA"), labels=c("16s rRNA Gene","16s rRNA"))+
  xlab("PCoA1: 39.1% var. explained")+
  ylab("PCoA2: 14.2% var. explained")+
  #The next line moves the "type" legend above the "treatment" legend.
  guides(colour = guide_legend(order = 2), 
         shape = guide_legend(order = 1))

ggsave("pcoa_RNA_bc.pdf",width=178,units="mm") 

#Permanova test using ADONIS function
# make a data frame from the sample_data
physeqActiveRNA.rfy.ra.df <- data.frame(sample_data(physeqActiveRNA.rfy.ra))
physeqActiveRNA.rfy.ra.df
# ADONIS test
library("vegan")
set.seed(1)
adonis(bc ~ Treatment, data = physeqActiveRNA.rfy.ra.df)

#Homogeneity of dispersions test
set.seed(1)
result <- betadisper(bc, physeqActiveRNA.rfy.ra.df$Treatment, type=c("centroid"))
set.seed(1)
permutest(result)



####

#Subset to only drought plus control:
physeqActiveRNA.rfy.ra.drctrl = subset_samples(physeqActiveRNA.rfy.ra, Treatment == "Drought" | Treatment == "Control")
sample_data(physeqActiveRNA.rfy.ra.drctrl)
#normalized Weighted Unifrac: considers relative abundance of taxa shared between samples.
bc <-distance(physeqActiveRNA.rfy.ra.drctrl, method="bray")

#Permanova test using ADONIS function
# make a data frame from the sample_data
physeqActiveRNA.rfy.ra.df <- data.frame(sample_data(physeqActiveRNA.rfy.ra.drctrl))
physeqActiveRNA.rfy.ra.df
# ADONIS test
set.seed(1)
adonis(bc ~ Treatment, data = physeqActiveRNA.rfy.ra.df)

#Homogeneity of dispersions test
set.seed(1)
result <- betadisper(bc, physeqActiveRNA.rfy.ra.df$Treatment, type=c("centroid"))
set.seed(1)
permutest(result)


###

#Subset to only nutrients plus control:
physeqActiveRNA.rfy.ra.nutctrl = subset_samples(physeqActiveRNA.rfy.ra, Treatment == "Nutrients" | Treatment == "Control")
sample_data(physeqActiveRNA.rfy.ra.nutctrl)
#normalized Weighted Unifrac: considers relative abundance of taxa shared between samples.
bc <-distance(physeqActiveRNA.rfy.ra.nutctrl, method="bray")

#Permanova test using ADONIS function
# make a data frame from the sample_data
physeqActiveRNA.rfy.ra.df <- data.frame(sample_data(physeqActiveRNA.rfy.ra.nutctrl))
physeqActiveRNA.rfy.ra.df
# ADONIS test
set.seed(1)
adonis(bc ~ Treatment, data = physeqActiveRNA.rfy.ra.df)

#Homogeneity of dispersions test
set.seed(1)
result <- betadisper(bc, physeqActiveRNA.rfy.ra.df$Treatment, type=c("centroid"))
set.seed(1)
permutest(result)

###

#Subset to only drought plus nutrients:
physeqActiveRNA.rfy.ra.nutdr = subset_samples(physeqActiveRNA.rfy.ra, Treatment == "Nutrients" | Treatment == "Drought")
sample_data(physeqActiveRNA.rfy.ra.nutdr)
#normalized Weighted Unifrac: considers relative abundance of taxa shared between samples.
bc <-distance(physeqActiveRNA.rfy.ra.nutdr, method="bray")

#Permanova test using ADONIS function
# make a data frame from the sample_data
physeqActiveRNA.rfy.ra.df <- data.frame(sample_data(physeqActiveRNA.rfy.ra.nutdr))
physeqActiveRNA.rfy.ra.df
# ADONIS test
set.seed(1)
adonis(bc ~ Treatment, data = physeqActiveRNA.rfy.ra.df)

#Homogeneity of dispersions test
set.seed(1)
result <- betadisper(bc, physeqActiveRNA.rfy.ra.df$Treatment, type=c("centroid"))
set.seed(1)
permutest(result)



#### Rank Abundance Curve for RNA versus DNA -------------------------------------

###determine which taxa are Phantom taxa

Bean_OTU.rfy <- otu_table(Bean.rfy)
Bean_OTU.rfy <- as.data.frame(Bean_OTU.rfy)
names(Bean_OTU.rfy) = gsub(pattern = "cDNA", replacement = "RNA", x = names(Bean_OTU.rfy))
#Make OTU table containing only 16sRNA
Bean_OTUsRNA <- Bean_OTU.rfy[,grep('RNA',names(Bean_OTU.rfy))]
#Make OTU table containing only 16sDNA
Bean_OTUsDNA <- Bean_OTU.rfy[,grep('DNA',names(Bean_OTU.rfy))]

#in RNA table, replace RNA>=1 as 1, and keep all zeros as zero.
Bean_OTUsRNA[Bean_OTUsRNA>=1] = 1
Bean_OTUsRNA[Bean_OTUsRNA<1] = 0
#In DNA table, replace all 0's with 1, and everything else with 0.
#First, make things that are currently >=1 into a higher number.
Bean_OTUsDNA[Bean_OTUsDNA>=1] = 2
#then, make everything that's currently 0 into a 1.
Bean_OTUsDNA[Bean_OTUsDNA<1] = 1
#then, make everything greater than 1 into 0.
Bean_OTUsDNA[Bean_OTUsDNA>1] = 0
#sum the two tables
Bean_RNAnoDNA <- data.frame(Bean_OTUsRNA+Bean_OTUsDNA)

#NOW, everthing with a 2 is a phantom taxa.
#replace everything with a 2 with "Phantom"...all else with "NotPhantom".
Bean_RNAnoDNA[Bean_RNAnoDNA=="2"] <- "Phantom"
Bean_RNAnoDNA[Bean_RNAnoDNA=="1"] <- "NotPhantom"
Bean_RNAnoDNA[Bean_RNAnoDNA=="0"] <- "NotPhantom"



#make figure: rank abundance curve for RNA

#make bean physeq object containing only RNA
Bean_OTU.rfy <- otu_table(Bean.rfy)
Bean_OTU.rfy <- as.data.frame(Bean_OTU.rfy)
#Make OTU table containing only 16sRNA
Bean_OTUsRNA <- Bean_OTU.rfy[,grep('cDNA',names(Bean_OTU.rfy))]
#convert to matrix, then phyloseq format
Bean_OTUsRNA <- data.matrix(Bean_OTUsRNA)
Bean_OTUsRNA = otu_table(Bean_OTUsRNA, taxa_are_rows = TRUE)
Bean_rare_RNA = phyloseq(Bean_OTUsRNA, Bean_TAX, Bean_METADATA)

# this converts taxa counts in each sample to a percentage
phyloTemp_BeanRNA = transform_sample_counts(Bean_rare_RNA, function(x) 1e+02 * x/sum(x))
clusterData_BeanRNA = psmelt(phyloTemp_BeanRNA)
clusterData_BeanRNA = filter(clusterData_BeanRNA, Abundance > 0)
# this is where the mean is calculated and the taxa to display is chosen
clusterAgg_BeanRNA = aggregate(Abundance ~ OTU,data=clusterData_BeanRNA,mean)
# filtering and picking the number to display
clusterAgg_BeanRNA = clusterAgg_BeanRNA[order(-clusterAgg_BeanRNA$Abundance),][1:100,]

#combine with data on Phantom taxa in case we want to label phantoms on the Figure
#ClusterAgg_Fig <- merge(clusterAgg_BeanRNA,Bean_RNAnoDNA,by.x="OTU",by.y="row.names",all=FALSE)
#rename columns
#names(ClusterAgg_Fig) <- c("OTU", "Abundance","Phantom")


Bean_RNAFig<-ggplot(clusterAgg_BeanRNA,aes(x=reorder(OTU,-Abundance),y=Abundance)) +
  geom_point() + 
  theme(axis.ticks = element_blank(), axis.text.x = element_blank()) +
  ylab("Relative Abundance")+
  xlab("OTU Rank")+
  ggtitle("Bean RNA Dataset")+
  coord_cartesian(ylim=c(0,5))
Bean_RNAFig


#make figure: rank abundance curve for DNA


#Make OTU table containing only 16sDNA
Bean_OTU.rfy <- otu_table(Bean.rfy)
Bean_OTU.rfy <- as.data.frame(Bean_OTU.rfy)
#Make OTU table containing only 16sDNA
names(Bean_OTU.rfy) = gsub(pattern = "cDNA", replacement = "RNA", x = names(Bean_OTU.rfy))
Bean_OTUsDNA <- Bean_OTU.rfy[,grep('DNA',names(Bean_OTU.rfy))]
#convert to matrix, then phyloseq format
Bean_OTUsDNA <- data.matrix(Bean_OTUsDNA)
Bean_OTUsDNA = otu_table(Bean_OTUsDNA, taxa_are_rows = TRUE)
Bean_rare_DNA = phyloseq(Bean_OTUsDNA, Bean_TAX, Bean_METADATA)


# this converts taxa counts in each sample to a percentage
phyloTemp_BeanDNA = transform_sample_counts(Bean_rare_DNA, function(x) 1e+02 * x/sum(x))
clusterData_BeanDNA = psmelt(phyloTemp_BeanDNA)
clusterData_BeanDNA = filter(clusterData_BeanDNA,Abundance > 0)
# this is where the mean is calculated and the taxa to display is chosen
clusterAgg_BeanDNA = aggregate(Abundance ~ OTU,data=clusterData_BeanDNA,mean)
# filtering and picking the number to display
clusterAgg_BeanDNA = clusterAgg_BeanDNA[order(-clusterAgg_BeanDNA$Abundance),][1:100,]

Bean_DNAFig<-ggplot(clusterAgg_BeanDNA,aes(x=reorder(OTU,-Abundance),y=Abundance)) +
  geom_point() + 
  theme(axis.ticks = element_blank(), axis.text.x = element_blank()) +
  ylab("Relative Abundance")+
  xlab("OTU Rank")+
  ggtitle("Bean DNA Dataset")+
  coord_cartesian(ylim=c(0,5))
Bean_DNAFig





