#Brad Fortunato


#install necessary packages;tidyverse for ggplot and complexheatmap

install.packages("tidyverse", repos = "http://cran.us.r-project.org")
library(tidyverse)
install.packages("BiocManager", repos = "http://cran.us.r-project.org")
library("ComplexHeatmap")

#read in aorta and kidney data, get headers from each and set those as column names

aorta <- read.table("C:/Users/b42na/Downloads/GSE158197_RawCounts_CS_200225_aorta.txt.gz", skip = 1,  header = FALSE, sep = "", dec = ".")
headers_a <- read.table("C:/Users/b42na/Downloads/GSE158197_RawCounts_CS_200225_aorta.txt.gz", header = F, nrows = 1, as.is = T)
kidney <- read.table("C:/Users/b42na/Documents/Broad Inst/GSE158197_RawCounts_BS_200803_kidney.txt", skip = 1,  header = FALSE, sep = "", dec = ".")
headers_k <- read.table("C:/Users/b42na/Documents/Broad Inst/GSE158197_RawCounts_BS_200803_kidney.txt", header = F, nrows = 1, as.is = T) 
colnames(aorta) <- headers_a
colnames(kidney) <- headers_k

#Picked sox17 gene to compare across the two WT groups; combine into one dataframe (rotating it)

data_ak_rp1 <- aorta[20,5:14] 
data_ak_rp1[nrow(data_ak_rp1) +1,] =  kidney[20,5:14]
df<-as.data.frame(t(data_ak_rp1)) 

#change names of columns

names(df) <- gsub("20", "gene_1", names(df))
names(df) <- gsub("2", "gene_2", names(df))
denom <- rownames(df)
df$names = denom

#scatter plot
p <- ggplot(df, aes(x = gene_1, y = gene_2)) + geom_point(color = "red") + geom_smooth(method = 'lm') + 
  geom_text(aes(label = names), vjust=2, size = 3)+
  labs(x = "Sox17 Expression: WT Aorta", y = "Sox17 Expression: WT Kidney", title = "Scatterplot: Expression of Sox17 in Wild Type Mus Muscula (Aorta vs. Kidney)" ) 
show(p)

#pearson correlation test
cor.test(df$gene_1, df$gene_2)

###################################################################################################

#new data frame; Mrpl15
df2 <- as.data.frame(t(aorta[32,5:14]))
names(df2) <- gsub("32", "gene_3", names(df2))

#add name column
denom <- rownames(df2)
df2$names = denom

#barplot

b <- ggplot(df2, aes(y = gene_3, x = names,label = "Mrpl15")) + geom_bar(stat = "identity", fill = "red") +   labs(x = "Wild Type Replicates (Aorta)", y = "Expression Across All Replicates", title = "Mrpl15 Expression Across all WT Replicates (Aorta)" ) 

b + theme(axis.text.x = element_text(angle = 45, hjust=1))


###################################################################################################

#Bonus: heatmap of top 50 highest expressed genes (aorta)
df3 <- as.matrix(aorta[,5:14])

#order by high to low, then top 50 cutoff
df3 = df3[order(df3[,1],decreasing=TRUE),] 
df3 <- df3[1:50,]

#heatmap
Heatmap(df3, name = "Expression")
