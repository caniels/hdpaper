####  Tom G. Caniels, 2021   
####  Analysis of NGS dataset of healthy donor B cells, specific for SARS-CoV-2 S
####  02-2021 All rights reserved

####  Library loading####
packages <- c("ggplot2","extrafont","RColorBrewer","gridExtra","patchwork","ggrepel","ggpubr","dplyr",
"plotrix","seqinr","plyr","dplyr","data.table","tidyverse","doParallel","colorspace",
"pheatmap","viridis","janitor","FactoMineR","factoextra","ggfortify","seqRFLP","Rtsne",
"cluster","networkD3","ggforce","GISTools","umap","gg.gap","plotly",  "data.table", "ggbeeswarm", "ggsci", 
"htmlwidgets", "webshot", "khroma", "scales")


# lapply(packages, install.packages, character.only = TRUE)

lapply(packages, library, character.only = TRUE)

####  Read files and filter data####
set.seed(69)
full <- read.csv(file = "F:/full_final_2021.csv", sep = ",")

full <- full[grepl("C",full$MULTI_ID) &
               full$V_sub_x_nunique == 1 &
               full$V_sub_y_nunique == 1 &
               full$mut.count.V.heavy_mean != "#N/A",]

full$V_sub_x_concat[full$V_sub_x_concat == "IGHV1-69D"] <- "IGHV1-69" 


full69 <- full[full$V_sub_x_concat == "IGHV1-69",]
fullnon69 <- full[!(full$V_sub_x_concat == "IGHV1-69"),]
full69vk <- full[full$V_sub_x_concat == "IGHV1-69" & full$V_sub_y_concat == "IGKV3-11",]
fullnon69vk <- full[!(full$V_sub_x_concat == "IGHV1-69" & full$V_sub_y_concat == "IGKV3-11"),]

names(full)[23] <- "CD27_CLR"
names(full)[24] <- "IgD_CLR"

donor <- read.csv("F:/from home/donoroverview.csv", sep = ";")

full$donor <- if(full$MULTI_ID == donor$hashtag_id) {donor$donor_id2}

full$donor <- case_when(full$MULTI_ID == "C0251" ~ "HD01",
                        full$MULTI_ID == "C0252" ~ "HD02",
                        full$MULTI_ID == "C0253" ~ "HD03",
                        full$MULTI_ID == "C0254" ~ "HD04",
                        full$MULTI_ID == "C0255" ~ "HD05",
                        full$MULTI_ID == "C0256" ~ "HD06",
                        full$MULTI_ID == "C0257" ~ "HD07",
                        full$MULTI_ID == "C0258" ~ "HD08",
                        full$MULTI_ID == "C0259" ~ "HD09",
                        full$MULTI_ID == "C0260" ~ "HD10")

full$age <- case_when(full$MULTI_ID == "C0251" ~ "41",
                      full$MULTI_ID == "C0252" ~ "54",
                      full$MULTI_ID == "C0253" ~ "42",
                      full$MULTI_ID == "C0254" ~ "61",
                      full$MULTI_ID == "C0255" ~ "65",
                      full$MULTI_ID == "C0256" ~ "39",
                      full$MULTI_ID == "C0257" ~ "28",
                      full$MULTI_ID == "C0258" ~ "58",
                      full$MULTI_ID == "C0259" ~ "59",
                      full$MULTI_ID == "C0260" ~ "51")

full$isotype <- case_when(full$IGH == "IGHA1" ~ "IgA",
                          full$IGH == "IGHA2" ~ "IgA",
                          full$IGH == "IGHG1" ~ "IgG",
                          full$IGH == "IGHG2" ~ "IgG",
                          full$IGH == "IGHG3" ~ "IgG",
                          full$IGH == "IGHD" ~ "IgD",
                          full$IGH == "IGHM" ~ "IgM")

colors2 <- c("#1b9e77", "#7570b3")
colors3 <- c("#1b9e77", "#d95f02", "#7570b3")
colors4 <- c("#1b9e77", "#d95f02", "#7570b3", "grey60")
colors4_2 <- c("#1b9e77", "#d95f02", "#7570b3","#e6ab02")
colors5 <- c("#1b9e77", "#d95f02", "#7570b3","#e6ab02", "grey60")
colors6 <- c("#332288", "#88CCEE","#DDCC77","#117733","#882255","#DDDDDD")

tolcol4 <- c("#0E84A4", "#F3CC5F", '#D88448', "#868282")
tolcol5 <- c("#0E84A4", "#F3CC5F", "#fdbf6f", "#66a61e", '#D88448', "#868282")


#### Pairing matrices ####
usage_scale <- colorRampPalette(c('#ffffd4','#fee391','#fec44f','#fe9929','#d95f0e','#993404'), bias = 2)
usage_scale2 <- colorRampPalette(c('#edf8fb','#bfd3e6','#9ebcda','#8c96c6','#8856a7','#810f7c'), bias = 2)
usage_scale3 <- colorRampPalette(c('#edf8fb','#ccece6','#99d8c9','#66c2a4','#2ca25f','#006d2c'), bias = 2)
usage_scale4 <- colorRampPalette(c("#F2F0F7", "#DADAEB", "#BCBDDC", "#9E9AC8", "#756BB1", "#54278F"), bias = 2)


full$mut.count.V.heavy_mean <- as.numeric(full$mut.count.V.heavy_mean) 
full$mut.frac.V.heavy_mean <- as.numeric(full$mut.frac.V.heavy_mean) 
names(full)[51] <- "v_gene"
names(full)[63] <- "v_gene2"

fullwo2 <- full[full$hash.ID != "Doublet",]

full$totalmuts <- as.numeric(full$mut.count.V.heavy_median) + 
  as.numeric(full$mut.count.V.light_median)

c_region <- read.csv("F:/count-c-region_2021.csv", sep = ",")
full <- merge(full, c_region, by = "barcode")

full$isotype2 <- case_when(full$cregion.name == "IGHA1" ~ "IgA",
                          full$cregion.name == "IGHA2" ~ "IgA",
                          full$cregion.name == "IGHG1" ~ "IgG",
                          full$cregion.name == "IGHG2" ~ "IgG",
                          full$cregion.name == "IGHG3" ~ "IgG",
                          full$cregion.name == "IGHD" ~ "IgD",
                          full$cregion.name == "IGHM" ~ "IgM",
                          full$cregion.name == "IGHM,IGHD" ~ "IgD/IgM")

sum(is.na(full$isotype2))


full$geno <- case_when(full$v_gene == "IGHV1-69" ~ "VH1-69",
                       full$v_gene != "IGHV1-69" ~ "non-VH1-69")

dekosky <- read.csv("F:/dekosky.txt", sep = "\t")

color42 <- c("#75bdf0","#9982f2","#ff915c","grey70")

full$v_gene2 <- gsub("\\;.*","",full$v_gene2)
full$v_gene2 <- str_remove(full$v_gene2, "[D]")


#ggplot(test, aes(x=cregion.name, fill = isotype)) + geom_bar(stat = "count")


####  Clustering analyses####
#full$totalmuts <- as.numeric(full$mut.count.V.heavy_median) + 
#  as.numeric(full$mut.count.V.light_median)

#full_pca <- full[c(23,24,94)]
#full_pca <- mutate_all(full_pca, function(x) as.numeric(as.character(x)))
#full_pca <- scale(full_pca, center = TRUE, scale = TRUE)


#k2 <- kmeans(full_pca, centers = 2, nstart = 100)
#fviz_cluster(k2, data = full_pca, geom = c("point","text"))
#kmeans_res <- factor(k2$cluster)
#full$cluster <- kmeans_res
#fviz_nbclust(full_pca, kmeans, method = "silhouette")


#### Cluster based on: IgD counts, CD27 counts, mutations HC, mutations LC ####

fullcluster <- full[c(23,24,89)]
fullcluster <- mutate_all(fullcluster, function(x) as.numeric(as.character(x)))
matrix <- dist(fullcluster, "euclidean")


silhouette <- c()
silhouette = c(silhouette, NA)
for(i in 2:10){
  pam_clusters = pam(as.matrix(matrix),
                     diss = TRUE,
                     k = i)
  silhouette = c(silhouette ,pam_clusters$silinfo$avg.width)
}


plot(1:10, silhouette,
     xlab = "Clusters",
     ylab = "Silhouette Width") + lines(1:10, silhouette)

pamfull <-  pam(matrix, k = 2)
fullcluster[pamfull$medoids, ]

rtsnedf <- Rtsne(matrix, is_distance = TRUE)
tsnedf <- rtsnedf$Y %>%
  data.frame() %>%
  setNames(c("X", "Y")) %>%
  mutate(cluster = factor(pamfull$clustering))

full$coord_X <- tsnedf$X
full$coord_Y <- tsnedf$Y
full$cluster_tsne <- tsnedf$cluster

full$pheno <- ifelse(full$cluster_tsne == "1", "naive", "memory")
sum(full$pheno == "naive")

ggplot(full, aes(x=coord_X, y=coord_Y)) +
  geom_point(aes(color = cluster_tsne)) +
  theme(panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.border = element_rect(colour = "black", fill = "transparent"), 
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12, face = "plain"),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 20, hjust = .5),
        aspect.ratio = 1)

ifelse(sum(full$pheno == "naive") + sum(full$pheno == "memory") == 1965, "OK.", "ERROR")

sum(full$pheno == "naive")
sum(full$pheno == "memory")


####  Number of cells per donor####
set.seed(69)
per_donor <- ggplot(full, aes(x = reorder(donor, desc(donor)), fill = pheno)) +
  geom_histogram(stat = "count", color = "black", alpha = 0.8) +
  scale_fill_manual(values = rev(colors2)) +
  theme(panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent", color = NA), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 20, color = "black"),
        axis.text.y = element_text(size = 20, color = "black"),
        axis.title.x = element_text(size = 20, face = "plain", color = "black"),
        axis.title.y = element_blank(),
        axis.ticks.length = unit(.2,"cm"),
        legend.text = element_text(size = 20),
        legend.position = c(0.8,0.9),
        legend.key = element_blank(),
        legend.title = element_blank(),
        plot.title = element_text(size = 20, hjust = .5),
        aspect.ratio = 2)  + coord_flip() +
  scale_y_continuous(expand = c(0,0), limits = c(0,400), breaks = seq(0,420,100)) +
  ylab("\nNumber of HC/LC pairs")

set.seed(69)
pop_pheno <- ggplot(full, aes(x= log10(as.numeric(CD27_CLR)), y = log10(as.numeric(IgD_CLR)), size = as.numeric(mut.frac.V.heavy_mean)*100)) +
  geom_point(shape = 16, aes(color = cluster_tsne), alpha = .7, position = position_jitter(seed = 42, width = 0.2)) +
  scale_size(breaks = c(0,2,4,6,8,10), range = c(2,7), name="IGHV mismatch\nto germline (%)") +
  theme(panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.border = element_rect(colour = "black", fill = "transparent"), 
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 20, color = "black"),
        axis.title = element_text(size = 20, face = "plain", color = "black"),
        legend.text = element_text(size = 20, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.key = element_blank(),
        aspect.ratio = 1) +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  ylab("anti-IgD CLR\n") + xlab("\nanti-CD27 CLR") +
  scale_color_manual(values = colors2, name = "Cluster") 

pop_pheno


donor_pheno_plot <- pop_pheno | per_donor
donor_pheno_plot


ggsave("F:/images 25-3/donor_pheno_plot.pdf", donor_pheno_plot, dpi = 300)

#### Analysis on isotype ####
fullnaive <- full[full$pheno == "naive",]
fullmem <- full[full$pheno == "memory",]
fullnaive$v_gene2 <- gsub("\\;.*","", fullnaive$v_gene2)
fullnaive$v_gene2 <- str_remove(fullnaive$v_gene2, "[D]")
fullmem$v_gene2 <- gsub("\\;.*","",fullmem$v_gene2)
fullmem$v_gene2 <- str_remove(fullmem$v_gene2, "[D]")


iso_colors <- c("#f4777f",
"#a5d5d8",
"#ffffe0",
"#ffbcaf",
"#73a2c6", "#DDDDDD")

naive_iso <- ggplot(fullnaive, aes(x = reorder(donor, desc(donor)), fill = isotype2)) +
  geom_histogram(stat = "count", color = "black", alpha = 0.8) +
  scale_fill_manual(values = iso_colors) +
  theme(panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent", color = NA), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 20, color = "black"),
        axis.text.y = element_text(size = 20, color = "black"),
        axis.title.x = element_text(size = 20, face = "plain", color = "black"),
        axis.title.y = element_blank(),
        axis.ticks.length = unit(.2,"cm"),
        legend.text = element_text(size = 20),
        legend.position = c(1.0,0.9),
        legend.key = element_blank(),
        legend.title = element_blank(),
        plot.title = element_text(size = 20, hjust = .5),
        aspect.ratio = 2)  + coord_flip() +
  scale_y_continuous(expand = c(0,0), limits = c(0,350), breaks = seq(0,420,100)) + scale_y_reverse() +
  scale_x_discrete(position = "top") +
  ylab("\nNumber of HC/LC pairs")

naive_iso


mem_iso <- ggplot(fullmem, aes(x = reorder(donor, desc(donor)), fill = isotype2)) +
  geom_histogram(stat = "count", color = "black", alpha = 0.8) +
  scale_fill_manual(values = iso_colors) +
  theme(panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent", color = NA), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 20, color = "black"),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 20, face = "plain", color = "black"),
        axis.title.y = element_blank(),
        axis.ticks.length = unit(.2,"cm"),
        legend.text = element_text(size = 20),
        legend.position = c(1,0.9),
        legend.key = element_blank(),
        legend.title = element_blank(),
        plot.title = element_text(size = 20, hjust = .5),
        aspect.ratio = 2)  + coord_flip() +
  scale_y_continuous(expand = c(0,0), limits = c(0,40), breaks = seq(0,40,10)) +
  ylab("\nNumber of HC/LC pairs") 
mem_iso


set.seed(69)
pop_iso <- ggplot(full, aes(x= log10(as.numeric(CD27_CLR)), y = log10(as.numeric(IgD_CLR)), size = as.numeric(mut.frac.V.heavy_mean)*100)) +
  geom_point(shape = 16, aes(color = isotype2), alpha = .7, position = position_jitter(seed = 42, width = 0.2)) +
  scale_size(breaks = c(0,2,4,6,8,10), range = c(2,7), name="IGHV mismatch\nto germline (%)") +
  theme(panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.border = element_rect(colour = "black", fill = "transparent"), 
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 20, color = "black"),
        axis.title = element_text(size = 20, face = "plain", color = "black"),
        legend.text = element_text(size = 20, color = "black"),
        legend.title = element_text(size = 20, color = "black"),
        legend.key = element_blank(),
        aspect.ratio = 1) +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  ylab("anti-IgD CLR\n") + xlab("\nanti-CD27 CLR") +
  scale_color_manual(values = iso_colors, name = "Isotype") +
  set.seed(69)

pop_iso

pop_iso

pops <- pop_iso | pop_pheno
pops
iso_pheno_plot <- pop_iso | naive_iso | mem_iso
iso_pheno_plot

barplot <- per_donor | naive_iso | mem_iso
barplot

ggsave("F:/images 25-3/iso_pheno_plot.pdf", iso_pheno_plot, dpi = 300)
ggsave("F:/images 25-3//pops.pdf", pops, dpi = 300)
ggsave("F:/images 25-3//barplots.pdf", barplot, dpi = 300)


####  Gene usage analyses####
kosky1 <- read.csv("F:/dekosky/kosky1.csv", sep = ";")
kosky2 <- read.csv("F:/dekosky/kosky2.csv", sep = ";")   
kosky3 <- read.csv("F:/dekosky/kosky3.csv", sep = ";")

kosky_naive <- as.data.frame(rbind(kosky1[,c(9,10)], kosky2[,c(9,10)], kosky3[,c(9,10)]))

##kosky3_vh <- aggregate(kosky3$Read_Count, by=list(kosky3$VH_Gene), FUN=sum)
##kosky3_vh$percent <- kosky3_vh$x/sum(kosky3_vh$x)*100

#### DEKOSKY NAIVE USAGE PNAS VS NGS NAIVE####

kosky1$VH_Gene[kosky1$VH_gene == "IGHV1-69-2"] <- "IGHV1-69" 
kosky2$VH_Gene[kosky2$VH_gene == "IGHV1-69-2"] <- "IGHV1-69" 
kosky3$VH_Gene[kosky3$VH_gene == "IGHV1-69-2"] <- "IGHV1-69" 

koskypair1 <- tabyl(kosky1, VH_gene, VL_gene)
kpm1 <- melt(koskypair1)
kpm1$vh <- paste0(kpm1$VH_gene, "/", kpm1$variable)
kpm1 <- kpm1[,-(1:2)]

koskypair2 <- tabyl(kosky2, VH_gene, VL_gene)
kpm2 <- melt(koskypair2)
kpm2$vh <- paste0(kpm2$VH_gene, "/", kpm2$variable)
kpm2 <- kpm2[,-(1:2)]

koskypair3 <- tabyl(kosky3, VH_gene, VL_gene)
koskypair3 <- koskypair3[-1,-2]
kpm3 <- melt(koskypair3)
kpm3$vh <- paste0(kpm3$VH_gene, "/", kpm3$variable)
kpm3 <- kpm3[,-(1:2)]

m1 <- merge(kpm1, kpm2, by = "vh", all = T)
m2 <- merge(m1, kpm3, by = "vh", all = T)
m2[is.na(m2)] <- "0"
names(m2) <- c("pair", "donor1", "donor2", "donor3")
m2$donor1 <- as.numeric(m2$donor1)
m2$donor2 <- as.numeric(m2$donor2)
m2$donor3 <- as.numeric(m2$donor3)

m2percent <- data.frame(pair = m2$pair,
                        donor1 = m2$donor1/sum(m2$donor1)*100,
                        donor2 = m2$donor2/sum(m2$donor2)*100,
                        donor3 = m2$donor3/sum(m2$donor3)*100)

m2percent$mean <- rowMeans(m2percent[,c(2:4)])
m2percent$sd <- apply(m2percent[2:4], 1, sd)

fullnaive$v_gene2 <- gsub("\\;.*","",fullnaive$v_gene2)
fullnaive$v_gene2 <- str_remove(fullnaive$v_gene2, "[D]")
fullmem$v_gene2 <- gsub("\\;.*","",fullmem$v_gene2)
fullmem$v_gene2 <- str_remove(fullmem$v_gene2, "[D]")


fullpairing <- tabyl(fullnaive, v_gene, v_gene2, donor)
fullpairing <- melt(fullpairing)
fullpairing <- fullpairing %>% spread(L1, value)
fullpairing <- as.data.frame(fullpairing)
fullpairing$pair <- paste0(fullpairing$v_gene, "/", fullpairing$variable)

fullpp <- data.frame(pair = fullpairing$pair, 
                     HD01 = fullpairing$HD01/sum(fullpairing$HD01)*100,
                     HD02 = fullpairing$HD02/sum(fullpairing$HD02)*100,
                     HD03 = fullpairing$HD03/sum(fullpairing$HD03)*100,
                     HD04 = fullpairing$HD04/sum(fullpairing$HD04)*100,
                     HD05 = fullpairing$HD05/sum(fullpairing$HD05)*100,
                     HD06 = fullpairing$HD06/sum(fullpairing$HD06)*100,
                     HD07 = fullpairing$HD07/sum(fullpairing$HD07)*100,
                     HD08 = fullpairing$HD08/sum(fullpairing$HD08)*100,
                     HD09 = fullpairing$HD09/sum(fullpairing$HD09)*100,
                     HD10 = fullpairing$HD10/sum(fullpairing$HD10)*100)

fullpp <- fullpp[,-6]

colSums(fullpp[2:10])

fullpp$mean <- rowMeans(fullpp[,c(2:10)])
fullpp$sd <- apply(fullpp[2:10], 1, sd)

koskysamples <- m2percent[,c(1,2,3,4)]
fullsamples <- fullpp[,c(1:10)]
paircompsamples <- merge(koskysamples, fullsamples, by = "pair", all = T)
pcompfilter <- paircompsamples[paircompsamples$HD01 != 0,] 
pcompfilter <- na.omit(pcompfilter)


koskypairspc <- m2percent[, c(1,5,6)]
fullpairspc <- fullpp[, c(1, 11, 12)]

paircomp <- merge(koskypairspc, fullpairspc, by = "pair", all = T)
pcfilter <- na.omit(paircomp)
write.table(pcompfilter, file = "F:/pccompfilter.csv", sep = ",")                

pcfilter$mean <- rowMeans(pcfilter[5:13])
pcfilter2 <- pcfilter %>% slice_max(mean, n=49)
write.table(pcfilter2, "clipboard", sep = "\t")  


#### DEKOSKY MEMORY SET PNAS ####
koskymem1 <- read.csv("F:/dekosky/koskymem1.csv", sep = ";")
koskymem2 <- read.csv("F:/dekosky/koskymem2.csv", sep = ";")   
koskymem1$VH_Gene[koskymem1$VH_gene == "IGHV1-69-2"] <- "IGHV1-69" 
koskymem2$VH_Gene[koskymem2$VH_gene == "IGHV1-69-2"] <- "IGHV1-69" 


koskymempair1 <- tabyl(koskymem1, VH_gene, VL_gene)
mem1 <- melt(koskymempair1)
mem1$vh <- paste0(mem1$VH_gene, "/", mem1$variable)
mem1 <- mem1[,-(1:2)]

koskymempair2 <- tabyl(koskymem2, VH_gene, VL_gene)
mem2 <- melt(koskymempair2)
mem2$vh <- paste0(mem2$VH_gene, "/", mem2$variable)
mem2 <- mem2[,-(1:2)]

memmerged <- merge(mem1, mem2, by = "vh", all = T)

memmerged[is.na(memmerged)] <- "0"
names(memmerged) <- c("pair", "donor1", "donor2")
memmerged$donor1 <- as.numeric(memmerged$donor1)
memmerged$donor2 <- as.numeric(memmerged$donor2)

memmergedpercent <- data.frame(pair = memmerged$pair,
                        donor1 = memmerged$donor1/sum(memmerged$donor1)*100,
                        donor2 = memmerged$donor2/sum(memmerged$donor2)*100)

memmergedpercent$mean <- rowMeans(memmergedpercent[,c(2:3)])
memmergedpercent$sd <- apply(memmergedpercent[2:3], 1, sd)


### scim paper
scim <- read.csv("F:/scimmunol.csv", sep = ";")
scim <- scim[-(1214:1215),]
scimtab <- tabyl(scim, VH.Germline, VL.Germline, Donor)
scimtab <- melt(scimtab)
scimtab <- scimtab %>% spread(L1, value)

scimtab$pair <- paste0(scimtab$VH.Germline, "/", scimtab$variable)
scimtab <- as.data.frame(scimtab[,-c(1:2)])
names(scimtab) <- c("P01", "P02", "P03", "P04", "P05", "P06", "P07", "P08", "pair")

colSums(scimtab[1:8])

scimtabready <- data.frame(pair = scimtab$pair, 
                           P01 = scimtab$P01/sum(scimtab$P01)*100,
                           P02 = scimtab$P02/sum(scimtab$P02)*100,
                           P03 = scimtab$P03/sum(scimtab$P03)*100,
                           P04 = scimtab$P04/sum(scimtab$P04)*100,
                           P05 = scimtab$P05/sum(scimtab$P05)*100,
                           P06 = scimtab$P06/sum(scimtab$P06)*100,
                           P07 = scimtab$P07/sum(scimtab$P07)*100,
                           P08 = scimtab$P08/sum(scimtab$P08)*100)
colSums(scimtabready[2:9])

scimtabready <- scimtabready[,-2]

memcomparison <- merge(memmergedpercent, scimtabready, by = "pair", all = T)
memcomp <- na.omit(memcomparison)
write.table(memcomp, file = "F:/memcomp_raw.csv", sep = ",")  
                     
              
#### comparison naive vs dekosky naive bar chart ####
naivecomp <- na.omit(paircomp)
names(naivecomp) <- c("pair", "mean_kosky", "sd_kosky", "mean_full", "sd_full")
naivecomp <- naivecomp %>% slice_max(mean_full, n = 10)
naivecomp2 <- melt(naivecomp)
naivecomp$kosky <- "kosky"
naivecomp$ngs <- "ngs"

t <- data.frame(pair = c(naivecomp[,"pair"], naivecomp[,"pair"]),
                set = c(naivecomp$kosky, naivecomp$ngs),
                mean = c(naivecomp$mean_kosky, naivecomp$mean_full),
                sd = c(naivecomp$sd_kosky, naivecomp$sd_full))

ggplot(t, aes(x = reorder(pair, -mean), y = mean, fill = set)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "black") +
  geom_errorbar(position = position_dodge(width = 0.9), aes(x=pair, ymin=mean, ymax=mean+sd), width=0.33) +
  theme(panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent", color = NA), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1),
        axis.title.y = element_text(size = 12, face = "plain", color = "black"),
        axis.title.x = element_blank(),
        axis.ticks = element_line(color = "black"),
        legend.text = element_text(size = 12),
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_blank(),
        legend.title = element_text(size = 12, color = "black"),
        legend.position = c(0.7, 0.5),
        plot.title = element_text(size = 12, hjust = .5),
        aspect.ratio = .5) +
  scale_y_continuous(expand = c(0,0), limits = c(0,12), breaks = seq(0,12,2))

#### heavy chain/light chain split ####

kosky1_vh <- tabyl(kosky1, VH_gene)
kosky1_vh$percent <- kosky1_vh$percent*100

kosky2_vh <- tabyl(kosky2, VH_gene)
kosky2_vh$percent <- kosky2_vh$percent*100

kosky3_vh <- tabyl(kosky3, VH_gene)
kosky3_vh$percent <- kosky3_vh$percent*100

colnames(kosky1_vh)[1] <- "v_gene"
colnames(kosky2_vh)[1] <- "v_gene"
colnames(kosky3_vh)[1] <- "v_gene"


usage <- tabyl(full, v_gene, v_gene2, donor)
usage_naive <- tabyl(fullnaive, v_gene, donor)

usage_naive_pc <- data.frame(v_gene = usage_naive$v_gene, 
                     HD01 = usage_naive$HD01/sum(usage_naive$HD01)*100,
                     HD02 = usage_naive$HD02/sum(usage_naive$HD02)*100,
                     HD03 = usage_naive$HD03/sum(usage_naive$HD03)*100,
                     HD04 = usage_naive$HD04/sum(usage_naive$HD04)*100,
                     HD05 = usage_naive$HD05/sum(usage_naive$HD05)*100,
                     HD06 = usage_naive$HD06/sum(usage_naive$HD06)*100,
                     HD07 = usage_naive$HD07/sum(usage_naive$HD07)*100,
                     HD08 = usage_naive$HD08/sum(usage_naive$HD08)*100,
                     HD09 = usage_naive$HD09/sum(usage_naive$HD09)*100,
                     HD10 = usage_naive$HD10/sum(usage_naive$HD10)*100)

usage_naive_pc <- usage_naive_pc[,-6]
colSums(usage_naive_pc[2:10])


hcusagenaive <- merge(kosky1_vh, kosky2_vh, by = "v_gene")
hcusagenaive <- merge(hcusagenaive, kosky3_vh, by = "v_gene")
hcusagenaive <- merge(hcusagenaive, usage_naive_pc, by = "v_gene")
hcusagenaive <- hcusagenaive[,-c(2,4,6)]
write.table(hcusagenaive, "clipboard", sep = "\t")


hcusagenaive$kosky_mean <- rowMeans(hcusagenaive[2:4])
hcusagenaive$kosky_sd <- apply(hcusagenaive[2:4], 1, sd)
hcusagenaive$ngs_mean <- rowMeans(hcusagenaive[5:13])
hcusagenaive$ngs_sd <- apply(hcusagenaive[5:13], 1, sd)


usage_naive10 <- hcusagenaive %>% slice_max(ngs_mean, n = 10)
usage_naive10$fold <- usage_naive10$ngs_mean/usage_naive10$kosky_mean
usage_naive10$log <- log10(usage_naive10$fold)
  
limit_heavy <- max(abs(usage_naive10$log)) * c(-1, 1)     

vh_usage <- ggplot(usage_naive10, aes(x = reorder(v_gene, -ngs_mean), y = ngs_mean, fill = log)) +
  geom_histogram(stat = "identity", color = "black", alpha = 0.8) +
  geom_errorbar(aes(x=v_gene, ymin=ngs_mean, ymax=ngs_mean+ngs_sd), width=0.33) +
  theme(panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent", color = NA), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1),
        axis.title.y = element_text(size = 12, face = "plain", color = "black"),
        axis.title.x = element_blank(),
        axis.ticks = element_line(color = "black"),
        legend.text = element_text(size = 12),
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_blank(),
        legend.title = element_text(size = 12, color = "black"),
        legend.position = c(0.7, 0.5),
        plot.title = element_text(size = 12, hjust = .5),
        aspect.ratio = .5)  + 
  scale_y_continuous(expand = c(0,0), limits = c(0,40), breaks = seq(0,40,5)) +
  scale_fill_gradient2(high = "blue", mid = "white", low = "red", 
                       breaks = c( -0.6, -0.3,0, 0.3, 0.6), 
                       midpoint = 0, limit = c(-.6, .6), 
                       labels = c("-4", "-2", "1", "2", "4")) +
  ylab("\nIGHV gene usage (%)") +
  guides(fill = guide_colorbar(barwidth = 12, barheight = 1, ticks = TRUE, frame.colour = "black", 
                               ticks.colour = "black", title = "Fold difference\nvs. naive repertoire", 
                               direction = "horizontal", title.position = "top", title.hjust = 0.5, title.vjust = 1))


write.table(usage_naive10, "F:/from home/HD project/usage_HC_naive_10.csv", sep = ";")
vh_usage

usage_naive10_supp <- usage_naive10[,c(1,5,6,7,8,9,10,11,12,13)]
usage_naive10_supp_melt <- melt(usage_naive10_supp)
usage_naive10_supp_melt$v_gene <- factor(usage_naive10_supp_melt$v_gene,levels = c("IGHV1-69", "IGHV3-23", "IGHV1-46", "IGHV3-30", "IGHV1-18", "IGHV3-11", "IGHV3-21", "IGHV4-39", "IGHV4-59", "IGHV5-51"))

vh_usage_supp <- ggplot(usage_naive10_supp_melt, aes(x = v_gene, y = value, color = variable )) +
  geom_point() +
  geom_point(data = usage_naive10_supp_melt %>% 
               group_by(v_gene) %>% 
               summarise(mean_usage = mean(value)), 
             mapping = aes(y = mean_usage), 
             size = 10, color = 'red', shape = "_") +
  theme(panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent", color = NA), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1),
        axis.title.y = element_text(size = 12, face = "plain", color = "black"),
        axis.title.x = element_blank(),
        axis.ticks = element_line(color = "black"),
        legend.text = element_text(size = 12),
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_blank(),
        legend.title = element_text(size = 12, color = "black"),
        legend.position = c(0.9, 0.7),
        plot.title = element_text(size = 12, hjust = .5),
        aspect.ratio = .5)  + 
  scale_y_continuous(expand = c(0,0), limits = c(0,40), breaks = seq(0,40,5)) +
  ylab("\nIGHV gene usage (%)")

vh_usage_supp

usagelight_naive10_supp <- usagelight_naive10[,c(1,5,6,7,8,9,10,11,12,13)]
usagelight_naive10_supp_melt <- melt(usagelight_naive10_supp)
usagelight_naive10_supp_melt$v_gene2 <- factor(usagelight_naive10_supp_melt$v_gene2,levels = c("IGKV3-11", "IGKV3-20", "IGKV4-1", "IGKV3-15", "IGKV1-39", "IGKV1-5", "IGLV2-14", "IGKV2-28", "IGLV3-25", "IGKV1-8"))


vhlight_usage_supp <- ggplot(usagelight_naive10_supp_melt, aes(x = v_gene2, y = value, color = variable )) +
  geom_point() +
  geom_point(data = usagelight_naive10_supp_melt %>% 
               group_by(v_gene2) %>% 
               summarise(mean_usage = mean(value)), 
             mapping = aes(y = mean_usage), 
             size = 10, color = 'red', shape = "_") +
  theme(panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent", color = NA), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1),
        axis.title.y = element_text(size = 12, face = "plain", color = "black"),
        axis.title.x = element_blank(),
        axis.ticks = element_line(color = "black"),
        legend.text = element_text(size = 12),
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_blank(),
        legend.title = element_text(size = 12, color = "black"),
        legend.position = c(0.9, 0.7),
        plot.title = element_text(size = 12, hjust = .5),
        aspect.ratio = .5)  + 
  scale_y_continuous(expand = c(0,0), limits = c(0,25), breaks = seq(0,40,5)) +
  ylab("\nIGHV gene usage (%)")

vhlight_usage_supp


usagenaive_supp <- vh_usage_supp | vhlight_usage_supp
usagenaive_supp

ggsave("F:/from home/HD project/Submission/Editorial/usagenaive_supp.pdf", usagenaive_supp, dpi = 300)



#### LIGHT CHAIN USAGE ####

kosky1_vl <- tabyl(kosky1, VL_gene)

kosky1_vl$percent <- kosky1_vl$percent*100

kosky2_vl <- tabyl(kosky2, VL_gene)
kosky2_vl$percent <- kosky2_vl$percent*100

kosky3_vl <- tabyl(kosky3, VL_gene)
kosky3_vl$percent <- kosky3_vl$percent*100

colnames(kosky1_vl)[1] <- "v_gene2"
colnames(kosky2_vl)[1] <- "v_gene2"
colnames(kosky3_vl)[1] <- "v_gene2"

usagelight_naive <- tabyl(fullnaive, v_gene2, donor)

usagelight_naive_pc <- data.frame(v_gene2 = usagelight_naive$v_gene2, 
                             HD01 = usagelight_naive$HD01/sum(usagelight_naive$HD01)*100,
                             HD02 = usagelight_naive$HD02/sum(usagelight_naive$HD02)*100,
                             HD03 = usagelight_naive$HD03/sum(usagelight_naive$HD03)*100,
                             HD04 = usagelight_naive$HD04/sum(usagelight_naive$HD04)*100,
                             HD05 = usagelight_naive$HD05/sum(usagelight_naive$HD05)*100,
                             HD06 = usagelight_naive$HD06/sum(usagelight_naive$HD06)*100,
                             HD07 = usagelight_naive$HD07/sum(usagelight_naive$HD07)*100,
                             HD08 = usagelight_naive$HD08/sum(usagelight_naive$HD08)*100,
                             HD09 = usagelight_naive$HD09/sum(usagelight_naive$HD09)*100,
                             HD10 = usagelight_naive$HD10/sum(usagelight_naive$HD10)*100)

usagelight_naive_pc <- usagelight_naive_pc[,-6]
colSums(usagelight_naive_pc[2:10])


hcusagelightnaive <- merge(kosky1_vl, kosky2_vl, by = "v_gene2")
hcusagelightnaive <- merge(hcusagelightnaive, kosky3_vl, by = "v_gene2")
hcusagelightnaive <- merge(hcusagelightnaive, usagelight_naive_pc, by = "v_gene2")
hcusagelightnaive <- hcusagelightnaive[,-c(2,4,6)]
write.table(hcusagelightnaive, "clipboard", sep = "\t")


hcusagelightnaive$kosky_mean <- rowMeans(hcusagelightnaive[2:4])
hcusagelightnaive$kosky_sd <- apply(hcusagelightnaive[2:4], 1, sd)
hcusagelightnaive$ngs_mean <- rowMeans(hcusagelightnaive[5:13])
hcusagelightnaive$ngs_sd <- apply(hcusagelightnaive[5:13], 1, sd)


usagelight_naive10 <- hcusagelightnaive %>% slice_max(ngs_mean, n = 10)
usagelight_naive10$fold <- usagelight_naive10$ngs_mean/usagelight_naive10$kosky_mean
usagelight_naive10$log <- log10(usagelight_naive10$fold)


limit_light <- max(abs(usagelight_naive10$log)) * c(-1, 1)     

vl_usage <- ggplot(usagelight_naive10, aes(x = reorder(v_gene2, -ngs_mean), y = ngs_mean, fill = log)) +
  geom_histogram(stat = "identity", color = "black", alpha = 0.8) +
  geom_errorbar(aes(x=v_gene2, ymin=ngs_mean, ymax=ngs_mean+ngs_sd), width=0.33) +
  theme(panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent", color = NA), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1),
        axis.title.y = element_text(size = 12, face = "plain", color = "black"),
        axis.title.x = element_blank(),
        axis.ticks = element_line(color = "black"),
        legend.text = element_text(size = 12),
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_blank(),
        legend.title = element_text(size = 12, color = "black"),
        legend.position = c(0.7, 0.5),
        plot.title = element_text(size = 12, hjust = .5),
        aspect.ratio = .5) + 
  scale_y_continuous(expand = c(0,0), limits = c(0,25), breaks = seq(0,25,5)) +
  scale_fill_gradient2(high = "blue", mid = "white", low = "red", 
                       breaks = c( -0.6, -0.3,0, 0.3, 0.6), 
                       midpoint = 0, limit = c(-.6, .6), 
                       labels = c("-4", "-2", "1", "2", "4")) +
  ylab("\nIGKV/IGLV gene usage (%)") +
  guides(fill = guide_colorbar(barwidth = 12, barheight = 1, ticks = TRUE, frame.colour = "black", 
                               ticks.colour = "black", title = "Fold difference\nvs. naive repertoire", 
                               direction = "horizontal", title.position = "top", title.hjust = 0.5, title.vjust = 1))
write.table(usagelight_naive10, "F:/from home/HD project/usage_LC_naive_10.csv", sep = ";")
vl_usage


usagenaive <- vh_usage | vl_usage
usagenaive

ggsave("F:/images 25-5/usagenaive.pdf", usagenaive, dpi = 300)


#### HEAVY CHAIN/LIGHT CHAIN USAGE MEMORY CELLS ####
usage_mem <- tabyl(fullmem, v_gene, donor)

usage_mem$sum <- rowSums(usage_mem[2:11])
usage_mem10 <- usage_mem %>% slice_max(usage_mem$sum, n = 10)
usage_mem10 <- usage_mem10[,-ncol(usage_mem10)]
usage_mem10 <- melt(usage_mem10)


usage_mem_plot <- ggplot(usage_mem10, aes(x = reorder(v_gene, -value), y = value , fill = variable )) +
                  geom_bar(stat = "identity", color = "black") +
                  theme(panel.background = element_rect(fill = "transparent"), 
                        plot.background = element_rect(fill = "transparent", color = NA), 
                        axis.line = element_line(colour = "black"),
                        axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1),
                        axis.text.y = element_text(size = 12, color = "black"),
                        axis.title.y = element_text(size = 12, face = "plain", color = "black"),
                        axis.title.x = element_blank(),
                        axis.ticks.length = unit(.2,"cm"),
                        legend.text = element_text(size = 12),
                        legend.key = element_blank(),
                        legend.title = element_blank(),
                        plot.title = element_text(size = 12, hjust = .5),
                        aspect.ratio = 1/1.5) +
                  guides(colour = guide_legend(override.aes = list(size=3))) +
                  #scale_fill_manual(values = color4) +
                  scale_y_continuous(expand = c(0,0), limits = c(0,30)) +
                  ylab("Number of sequences\n")
usage_mem_plot



usage_mem_light <- tabyl(fullmem, v_gene2, donor)

usage_mem_light$sum <- rowSums(usage_mem_light[2:11])
usage_mem_light10 <- usage_mem_light %>% slice_max(usage_mem_light$sum, n = 10)
usage_mem_light10 <- usage_mem_light10[,-ncol(usage_mem_light10)]
usage_mem_light10 <- usage_mem_light10[-c(11,12),]
usage_mem_light10 <- melt(usage_mem_light10)


usage_mem_light_plot <- ggplot(usage_mem_light10, aes(x = reorder(v_gene2, -value), y = value , fill = variable )) +
  geom_bar(stat = "identity", color = "black") +
  theme(panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent", color = NA), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 12, face = "plain", color = "black"),
        axis.title.x = element_blank(),
        axis.ticks.length = unit(.2,"cm"),
        legend.text = element_text(size = 12),
        legend.key = element_blank(),
        legend.title = element_blank(),
        plot.title = element_text(size = 12, hjust = .5),
        aspect.ratio = 1/1.5) +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  #scale_fill_manual(values = color4) +
  scale_y_continuous(expand = c(0,0), limits = c(0,30)) +
  ylab("Number of sequences\n")


usage_mem_final<- usage_mem_plot | usage_mem_light_plot
usage_mem_final

ggsave("F:/images 25-5/usagemem.pdf", usage_mem_final, dpi = 300)

###########
koskymem1_tbl <- tabyl(koskymem1, VH_gene)
koskymem1_tbl$percent <- koskymem1_tbl$percent*100

koskymem2_tbl <- tabyl(koskymem2, VH_gene)
koskymem2_tbl$percent <- koskymem2_tbl$percent*100

koskymem <- merge(koskymem1_tbl, koskymem2_tbl, by = "VH_gene")

koskymem$mean <- rowMeans(koskymem[c(3,5)])
koskymem$sd <- apply(koskymem[c(3,5)], 1, sd)
koskymem10 <- koskymem %>% slice_max(mean, n=10)

usage_mempc <- data.frame(usage_mem, pc =usage_mem$sum/sum(usage_mem$sum)*100)
usage_mempc <- usage_mempc %>% slice_max(pc, n = 10)
names(koskymem)[1] <- "v_gene"

usage_memvskos <- merge(usage_mempc, koskymem, by = "v_gene")
usage_memvskos <- usage_memvskos[,c(1,13,18)]
names(usage_memvskos) <- c("v_gene", "mean_ngs", "mean_kosky")
usage_memvskos <- usage_memvskos %>% arrange(mean_ngs)
usage_memvskos <- melt(usage_memvskos)

memvskosplot <- ggplot(usage_memvskos, aes(x=reorder(v_gene, -value), y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent", color = NA), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 12, face = "plain", color = "black"),
        axis.title.x = element_blank(),
        axis.ticks.length = unit(.2,"cm"),
        legend.text = element_text(size = 12),
        legend.key = element_blank(),
        legend.title = element_blank(),
        plot.title = element_text(size = 12, hjust = .5),
        aspect.ratio = 1/1.5) +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  #scale_fill_manual(values = color4) +
  scale_y_continuous(expand = c(0,0), limits = c(0,15)) +
  ylab("IGHV gene usage (%)")

memvskosplot


##########
koskymem1_light_tbl <- tabyl(koskymem1, VL_gene)
koskymem1_light_tbl$percent <- koskymem1_light_tbl$percent*100

koskymem2_light_tbl <- tabyl(koskymem2, VL_gene)
koskymem2_light_tbl$percent <- koskymem2_light_tbl$percent*100

koskymem_light <- merge(koskymem1_light_tbl, koskymem2_light_tbl, by = "VL_gene")
names(koskymem_light)[1] <- "v_gene2"
koskymem_light$mean <- rowMeans(koskymem_light[c(3,5)])
koskymem_light$sd <- apply(koskymem_light[c(3,5)], 1, sd)
koskymem_light10 <- koskymem_light %>% slice_max(mean, n=10)

usage_mem_light_pc <- data.frame(usage_mem_light, pc = usage_mem_light$sum/sum(usage_mem_light$sum)*100)
usage_mem_light_pc <- usage_mem_light_pc %>% slice_max(pc, n = 10)

usage_memvskos_light <- merge(usage_mem_light_pc, koskymem_light, by = "v_gene2")
usage_memvskos_light <- usage_memvskos_light[,c(1,13,18)]
names(usage_memvskos_light) <- c("v_gene", "mean_ngs", "mean_kosky")
usage_memvskos_light <- usage_memvskos_light[-c(11,12),]

usage_memvskos_light$v_gene <- factor(usage_memvskos_light$v_gene, levels = usage_memvskos_light$v_gene[order(usage_memvskos_light$mean_ngs)])

usage_memvskos_light <- melt(usage_memvskos_light)

memvskosplotlight <-ggplot(usage_memvskos_light, aes(x=reorder(v_gene, -value/3), y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent", color = NA), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 12, face = "plain", color = "black"),
        axis.title.x = element_blank(),
        axis.ticks.length = unit(.2,"cm"),
        legend.text = element_text(size = 12),
        legend.key = element_blank(),
        legend.title = element_blank(),
        plot.title = element_text(size = 12, hjust = .5),
        aspect.ratio = 1/1.5) +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  #scale_fill_manual(values = color4) +
  scale_y_continuous(expand = c(0,0), limits = c(0,15)) +
  ylab("IGHV gene usage (%)")

memvskosplot | memvskosplotlight

memory_vs_kosky_plot <- memvskosplot | memvskosplotlight

supp_mem_usage <- usage_mem_final / memory_vs_kosky_plot


ggsave("F:/images 25-5/final_usagemem_supp.pdf", supp_mem_usage, dpi = 300)

#### usage in plot ####

set.seed(69)
heavypop <- ggplot(full, aes(x= log10(as.numeric(CD27_CLR)), y = log10(as.numeric(IgD_CLR)), size = as.numeric(mut.frac.V.heavy_mean)*100)) +
  geom_jitter(shape = 16, aes(color = v_gene == "IGHV1-69"), width = 0.25, alpha = .7) +
  scale_size(breaks = c(0,2,4,6,8,10), range = c(2,7), name="IGHV mismatch\nto germline (%)") +
  theme(panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.border = element_rect(colour = "black", fill = "transparent"), 
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 12, face = "plain", color = "black"),
        legend.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 12, color = "black"),
        legend.key = element_blank(),
        aspect.ratio = 1) +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  ylab("anti-IgD CLR\n") + xlab("\nanti-CD27 CLR") +
  scale_color_manual(values = tolcol4, name = "Isotype")

set.seed(69)
lightpop <- ggplot(full, aes(x= log10(as.numeric(CD27_CLR)), y = log10(as.numeric(IgD_CLR)), size = as.numeric(mut.frac.V.heavy_mean)*100)) +
  geom_jitter(shape = 16, aes(color = v_gene2 == "IGKV3-11"), width = 0.25, alpha = .7) +
  scale_size(breaks = c(0,2,4,6,8,10), range = c(2,7), name="IGHV mismatch\nto germline (%)") +
  theme(panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.border = element_rect(colour = "black", fill = "transparent"), 
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 12, face = "plain", color = "black"),
        legend.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 12, color = "black"),
        legend.key = element_blank(),
        aspect.ratio = 1) +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  ylab("anti-IgD CLR\n") + xlab("\nanti-CD27 CLR") +
  scale_color_manual(values = colors4, name = "Isotype")

  
heavypop + lightpop


ggsave("F:/images 25-5/test_pops.pdf", dpi = 300)




### !!!! WATCH HERE FOR NAIVE ONLY !!!! ####
pairingnaive <- as.data.frame(tabyl(fullnaive, v_gene, v_gene2))
rownames(pairingnaive) <- pairingnaive$v_gene
pairingnaive <- pairingnaive[,-1]
pairingnaive$sums <- rowSums(pairingnaive)
pairingnaive <- pairingnaive %>% slice_max(sums, n=7)
pairingnaive <- pairingnaive[,-ncol(pairingnaive)]
row_pairingnaive <- rownames(pairingnaive)
col_pairingnaive <- colnames(pairingnaive)
pairingnaive2 <- t(pairingnaive)
colnames(pairingnaive2) <- row_pairingnaive
rownames(pairingnaive2) <- col_pairingnaive
pairingnaive2 <- as.data.frame(pairingnaive2)
pairingnaive2$sums <- rowSums(pairingnaive2)
pairingnaive2 <- pairingnaive2 %>% slice_max(sums, n=7)
pairingnaive2 <- pairingnaive2[,-ncol(pairingnaive2)]

pairingnaive_perc <- pairingnaive2
for(col in names(pairingnaive2)) {
  pairingnaive_perc[col] = (pairingnaive2[col] / sum(pairingnaive2[col]) * 100)
}



pheatmap(pairingnaive2, color = usage_scale2(100), cellwidth = 30, fontsize = 16,
         fontsize_col = 20, fontsize_row = 20,
         cellheight = 30, cluster_rows = FALSE, cluster_cols = FALSE, number_format = "%.0f",
         cutree_cols = 1, cutree_rows = 1, angle_col = 45, display_numbers = T, number_color = "black",
         border_color = "black")

#### !!!! WATCH HERE FOR MEMORY ONLY !!!! ####

pairingmem <- as.data.frame(tabyl(fullmem, v_gene, v_gene2))
rownames(pairingmem) <- pairingmem$v_gene
pairingmem <- pairingmem[,-1]
pairingmem$sums <- rowSums(pairingmem)
pairingmem <- pairingmem %>% slice_max(sums, n=7)
pairingmem <- pairingmem[,-ncol(pairingmem)]
row_pairingmem <- rownames(pairingmem)
col_pairingmem <- colnames(pairingmem)
pairingmem2 <- t(pairingmem)
colnames(pairingmem2) <- row_pairingmem
rownames(pairingmem2) <- col_pairingmem
pairingmem2 <- as.data.frame(pairingmem2)
pairingmem2$sums <- rowSums(pairingmem2)
pairingmem2 <- pairingmem2 %>% slice_max(sums, n=7)
pairingmem2 <- pairingmem2[,-ncol(pairingmem2)]

pairingmem_perc <- pairingmem2
for(col in names(pairingmem2)) {
  pairingmem_perc[col] = (pairingmem2[col] / sum(pairingmem2[col]) * 100)
}



pheatmap(pairingmem2, color = usage_scale2(100), cellwidth = 30, fontsize = 16,
         fontsize_col = 20, fontsize_row = 20,
         cellheight = 30, cluster_rows = FALSE, cluster_cols = FALSE, number_format = "%.0f",
         cutree_cols = 1, cutree_rows = 1, angle_col = 45, display_numbers = T, number_color = "black",
         border_color = "black")


#### ABDAB pairing ####
abdab <- read.csv('F:/covabdab24-3.csv', sep = ";")
abdab$epitope <- case_when(abdab$Protein...Epitope == "S; RBD" ~ "RBD",
                           abdab$Protein...Epitope != "S; RBD" ~ "non-RBD")
abdab$Heavy.V.Gene[abdab$Heavy.V.Gene == "IGHV1-69-2"] <- "IGHV1-69"
abdab69 <- abdab[abdab$Heavy.V.Gene == "IGHV1-69",]

abdabnonrbd <- abdab[abdab$epitope == "non-RBD",]

abdabpairing <- as.data.frame(tabyl(abdabnonrbd, Heavy.V.Gene, Light.V.Gene))
rownames(abdabpairing) <- abdabpairing$Heavy.V.Gene
abdabpairing <- abdabpairing[,-1]
abdabpairing$sums <- rowSums(abdabpairing)
abdabpairing <- abdabpairing %>% slice_max(sums, n=7)
abdabpairing <- abdabpairing[,-ncol(abdabpairing)]
row_abdabpairing <- rownames(abdabpairing)
col_abdabpairing <- colnames(abdabpairing)
abdabpairing2 <- t(abdabpairing)
colnames(abdabpairing2) <- row_abdabpairing
rownames(abdabpairing2) <- col_abdabpairing
abdabpairing2 <- as.data.frame(abdabpairing2)
abdabpairing2$sums <- rowSums(abdabpairing2)
abdabpairing2 <- abdabpairing2 %>% slice_max(sums, n=7)
abdabpairing2 <- abdabpairing2[,-ncol(abdabpairing2)]


hm_abdab <- pheatmap(abdabpairing2, color = usage_scale3(100), cellwidth = 30, fontsize = 16,
         fontsize_col = 20, fontsize_row = 20,
         cellheight = 30, cluster_rows = FALSE, cluster_cols = FALSE, number_format = "%.0f",
         cutree_cols = 1, cutree_rows = 1, angle_col = 45, display_numbers = T, number_color = "black",
         border_color = "black")

hm_abdab

#### Science Immunology paper pairing ####
scim <- read.csv("F:/scimmunol.csv", sep = ";")
scim <- scim[-(1214:1215),]

scim$epitope <- case_when(str_count(scim$Antigenic.Site, "RBD") == TRUE ~ "RBD",
                          str_count(scim$Antigenic.Site, "S2") == TRUE ~ "S2",
                          str_count(scim$Antigenic.Site, "NTD") == TRUE ~ "NTD")

scimpairing <- as.data.frame(tabyl(scim, VH.Germline, VL.Germline))
rownames(scimpairing) <- scimpairing$VH.Germline
scimpairing <- scimpairing[,-1]
scimpairing$sums <- rowSums(scimpairing)
scimpairing <- scimpairing %>% slice_max(sums, n=8)
scimpairing <- scimpairing[,-ncol(scimpairing)]
row_scimpairing <- rownames(scimpairing)
col_scimpairing <- colnames(scimpairing)
scimpairing <- scimpairing[-nrow(scimpairing),]
scimpairing2 <- t(scimpairing)
colnames(scimpairing2) <- row_scimpairing
rownames(scimpairing2) <- col_scimpairing
scimpairing2 <- as.data.frame(scimpairing2)
scimpairing2$sums <- rowSums(scimpairing2)
scimpairing2 <- scimpairing2 %>% slice_max(sums, n=8)
scimpairing2 <- scimpairing2[-nrow(scimpairing2),]
scimpairing2 <- scimpairing2[,-ncol(scimpairing2)]

hm_scim <- pheatmap(scimpairing2, color = usage_scale(100), cellwidth = 30, fontsize = 16,
         fontsize_col = 20, fontsize_row = 20,
         cellheight = 30, cluster_rows = FALSE, cluster_cols = FALSE, number_format = "%.0f",
         cutree_cols = 1, cutree_rows = 1, angle_col = 45, display_numbers = T, number_color = "black",
         border_color = "black")

hm_scim

ggplot(scim, aes(x=Antigenic.Site, fill = VL.Germline == "IGKV3-20" | VL.Germline == "IGKV3-11")) + geom_bar(stat = "count")


#### DEKOSKY pairing ####

koskypairing <- as.data.frame(tabyl(dekosky, Heavy_V.Gene, Light_V.Gene))
rownames(koskypairing) <- koskypairing$Heavy_V.Gene
koskypairing <- koskypairing[,-1]
koskypairing$sums <- rowSums(koskypairing)
koskypairing <- koskypairing %>% slice_max(sums, n=12)
koskypairing <- koskypairing[,-ncol(koskypairing)]
row_koskypairing <- rownames(koskypairing)
col_koskypairing <- colnames(koskypairing)
koskypairing <- koskypairing[-nrow(koskypairing),]
koskypairing2 <- t(koskypairing)
colnames(koskypairing2) <- rownames(koskypairing)
rownames(koskypairing2) <- colnames(koskypairing)
koskypairing2 <- as.data.frame(koskypairing2)
koskypairing2$sums <- rowSums(koskypairing2)
koskypairing2 <- koskypairing2 %>% slice_max(sums, n=11)
koskypairing2 <- koskypairing2[,-ncol(koskypairing2)]
koskypairing2 <- as.matrix(koskypairing2)


hm_kosky <- pheatmap(koskypairing2, color = usage_scale3(100), cellwidth = 30, fontsize = 16,
                    fontsize_col = 20, fontsize_row = 20,
                    cellheight = 30, cluster_rows = FALSE, cluster_cols = FALSE, number_format = "%.0f",
                    cutree_cols = 1, cutree_rows = 1, angle_col = 45, display_numbers = T, number_color = "black",
                    border_color = "black")

ggsave("F:/images 25-3/heatmaps/heatmap_ngs.pdf", hm_full, dpi = 300)
#ggsave("F:/images 25-3/heatmaps/heatmap_abdab.pdf", hm_abdab, dpi = 300)
ggsave("F:/images 25-3/heatmaps/heatmap_scim.pdf", hm_scim, dpi = 300)
ggsave("F:/images 25-3/heatmaps/heatmap_kosky.pdf", hm_kosky, dpi = 300)



#### PAIRING MEMORY AND NAIVE SEPARATELY ####

mempairing <- as.data.frame(tabyl(fullmem, v_gene, v_gene2))
rownames(mempairing) <- mempairing$v_gene
mempairing <- mempairing[,-1]
mempairing$sums <- rowSums(mempairing)
mempairing <- mempairing %>% slice_max(sums, n=12)
mempairing <- mempairing[,-ncol(mempairing)]
row_mempairing <- rownames(mempairing)
col_mempairing <- colnames(mempairing)
mempairing <- mempairing[-nrow(mempairing),]
mempairing2 <- t(mempairing)
colnames(mempairing2) <- rownames(mempairing)
rownames(mempairing2) <- colnames(mempairing)
mempairing2 <- as.data.frame(mempairing2)
mempairing2$sums <- rowSums(mempairing2)
mempairing2 <- mempairing2 %>% slice_max(sums, n=11)
mempairing2 <- mempairing2[,-ncol(mempairing2)]
mempairing2 <- as.matrix(mempairing2)


hm_mem <- pheatmap(mempairing2, color = usage_scale3(100), cellwidth = 30, fontsize = 16,
                     fontsize_col = 20, fontsize_row = 20,
                     cellheight = 30, cluster_rows = FALSE, cluster_cols = FALSE, number_format = "%.0f",
                     cutree_cols = 1, cutree_rows = 1, angle_col = 45, display_numbers = T, number_color = "black",
                     border_color = "black")



naivepairing <- as.data.frame(tabyl(fullnaive, v_gene, v_gene2))
rownames(naivepairing) <- naivepairing$v_gene
naivepairing <- naivepairing[,-1]
naivepairing$sums <- rowSums(naivepairing)
naivepairing <- naivepairing %>% slice_max(sums, n=12)
naivepairing <- naivepairing[,-ncol(naivepairing)]
row_naivepairing <- rownames(naivepairing)
col_naivepairing <- colnames(naivepairing)
naivepairing <- naivepairing[-nrow(naivepairing),]
naivepairing2 <- t(naivepairing)
colnames(naivepairing2) <- rownames(naivepairing)
rownames(naivepairing2) <- colnames(naivepairing)
naivepairing2 <- as.data.frame(naivepairing2)
naivepairing2$sums <- rowSums(naivepairing2)
naivepairing2 <- naivepairing2 %>% slice_max(sums, n=11)
naivepairing2 <- naivepairing2[,-ncol(naivepairing2)]
naivepairing2 <- as.matrix(naivepairing2)


hm_naive <- pheatmap(naivepairing2, color = usage_scale3(100), cellwidth = 30, fontsize = 16,
                   fontsize_col = 20, fontsize_row = 20,
                   cellheight = 30, cluster_rows = FALSE, cluster_cols = FALSE, number_format = "%.0f",
                   cutree_cols = 1, cutree_rows = 1, angle_col = 45, display_numbers = T, number_color = "black",
                   border_color = "black")

#### CDRH3 with kosky naive ####
colnames(kosky_naive) <- c("cdrh3_len", "x")
kosky_cdr3 <- tabyl(kosky_naive$cdrh3_len)
kosky_cdr3 <- kosky_cdr3[c(1:23),]
colnames(kosky_cdr3) <- c("cdrh3_len", "kosky", "percent")
kosky_cdr3$cdrh3_len <- as.factor(kosky_cdr3$cdrh3_len)

fullnaive$cdrh3_len <- as.factor(fullnaive$cdrh3_len)
full_cdr3 <- tabyl(fullnaive$cdrh3_len)
full_cdr3 <- full_cdr3[,-3]
colnames(full_cdr3) <- c("cdrh3_len", "full")
full_cdr3 <- merge(full_cdr3, kosky_cdr3, by = "cdrh3_len")
full_cdr3 <- full_cdr3[,-4]

full_cdr3$cdrh3_len <- as.factor(full_cdr3$cdrh3_len)
cdr4 <- melt(full_cdr3)
write.table(cdr4, "clipboard", sep =  "\t")


cdr <- read.csv("F:/cdrh3_naive_final2.csv", sep = ",")
cdr$cdrh3_len <- as.factor(cdr$cdrh3_len)
cdr2 <- melt(cdr)

ngs_cdr3_leaf <- ggplot(cdr, aes(x = as.numeric(cdrh3_len)+3)) +
  geom_histogram(aes(y = -..density.., weight = ngs), fill = "#7b3294", binwidth = 1, alpha = .4, color = "black", position = "dodge") +
  geom_histogram(aes(y = ..density.., weight = kosky), fill = "#008837", binwidth = 1, alpha = .4, color = "black", position = "dodge") +
  #geom_density(aes(y = -..density.., weight = value), alpha = 0.4, adjust = .5) +
  theme(panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent", color = NA), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(size = 12, face = "plain", color = "black"),
        axis.title.y = element_blank(),
        axis.ticks = element_line(color = "black"),
        legend.text = element_text(size = 12),
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_blank(),
        legend.title = element_text(size = 12, color = "black"),
        legend.position = c(0.9, 0.9),
        plot.title = element_text(size = 12, hjust = .5),
        aspect.ratio = .5) +
  scale_fill_manual(values = c("#7b3294","#008837")) +
  scale_x_continuous(expand = c(0,0), breaks = seq(0,30,2), limits = c(3.5,27.5)) +
  scale_y_continuous(expand = c(0,0), limits = c(-0.13, 0.13)) +
  ylab("CDRH3 length (aa)") + xlab("Density")

ngs_cdr3_leaf

ngs_cdr3_leaf2 <- ggplot(cdr, aes(x = as.numeric(cdrh3_len)+3)) +
  geom_histogram(aes(y = -..density.., weight = ngs), fill = "#7b3294", binwidth = 1, alpha = .4, color = "black", position = "dodge") +
  geom_histogram(aes(y = ..density.., weight = kosky), fill = "#008837", binwidth = 1, alpha = .4, color = "black", position = "dodge") +
  theme(panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent", color = NA), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(size = 12, face = "plain", color = "black"),
        axis.title.y = element_blank(),
        axis.ticks = element_line(color = "black"),
        legend.text = element_text(size = 12),
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_blank(),
        legend.title = element_text(size = 12, color = "black"),
        legend.position = c(0.9, 0.9),
        plot.title = element_text(size = 12, hjust = .5),
        aspect.ratio = .5) +
  scale_fill_manual(values = c("#7b3294","#008837")) +
  scale_x_continuous(expand = c(0,0), breaks = seq(0,30,2), limits = c(3.5,27.5)) +
  scale_y_continuous(expand = c(0,0), limits = c(-0.13, 0.13)) +
  xlab("CDRH3 length (aa)") + ylab("Density")

ngs_cdr3_leaf2

ggplot(cdr3_in2, aes(x=variable, y = value, fill = variable)) +
  geom_violin(scale = "width", width = 0.8, trim = F, adjust = 2) +
  stat_summary(fun=mean, geom = "point", size =2) +
  #scale_fill_manual(values = c("#485c99", "#c5462c", "#fbbe48")) +
  theme(panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent", color = NA), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(size = 12, face = "plain", color = "black"),
        axis.title.y = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.ticks.length=unit(.25, "cm"),
        legend.text = element_text(size = 12),
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_blank(),
        legend.title = element_text(size = 12, color = "black"),
        legend.position = c(0.9, 0.9),
        plot.title = element_text(size = 12, hjust = .5),
        aspect.ratio = 1.5) +
  scale_fill_manual(values = c("#7b3294","#008837")) +
  scale_y_continuous(expand = c(0,0), limits = c(0,30))



####  Network diagram, HC/LC pairing####
hclc <- tabyl(full, v_gene, v_gene2)
hc <- hclc[,1]
hclc2 <- data.frame(hclc[,c(2:ncol(hclc))])
rownames(hclc2) <- hc
hclc2$sums <- rowSums(hclc2)

hclc5 <- hclc2 %>% slice_max(sums, n=15)
hclc5 <- hclc5[,-62]


lcs <- names(hclc5) 

for (i in seq_along(lcs)){
  hclc5[,i] <- hclc5[,i]/sum(hclc5[i,])*100
}

pheatmap(hclc5, color = viridis(10), cellwidth = 5, 
         cellheight = 5, cluster_rows = FALSE, cluster_cols = FALSE, 
         cutree_cols = 1, cutree_rows = 1, angle_col = 45,
        border_color = "black")

rowSums(hclc)


hclc5 <- hclc5[,-62]

hclcf <- t(hclc5)

hclcf <- data.frame(hclcf[,c(1:ncol(hclcf))])  
hclcf$sums <- rowSums(hclcf)
hclcf <- hclcf %>% slice_max(sums, n = 5)
hclcf <- hclcf[,-6]

hclcf <- t(hclcf)


links <- pairings12 %>% 
  as.data.frame() %>% 
  rownames_to_column(var="source") %>% 
  gather(key="target", value="value", -1) %>%
  filter(value != 0)

nodes <- data.frame(
  name=c(as.character(links$source), as.character(links$target)) %>% 
    unique()
)

links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1

my_color <- 'd3.scaleOrdinal() .domain(["IGHV1-69", "IGHV3-23", "IGHV3-30", "IGHV4-39", "IGHV1-46"]) .range(["#66c2a5","#fc8d62","#8da0cb","#e78ac3","red", "red", "red", "red", "red", "red"])'



my_color2 <- 'd3.scaleOrdinal() .domain(["RBD", "NTD", "S2", "UND", "yes", "1", "2", "3", "4", "5"]) 
.range([ "#3B4992CC", "#EE0000CC", "#008B45CC", "#631879CC", "#cccccc", "#cccccc", "#EE0000CC", "#3B4992CC", "#008B45CC", "#631879CC"])'

links <- links[-c(32:42),]

links$group <- as.factor(c("NTD","NTD","NTD","NTD","NTD","NTD",
                           "RBD","RBD","RBD","RBD","RBD","RBD","RBD",
                           "S2","S2","S2","S2","S2","S2","S2","S2","S2","S2",
                           "UND","UND","UND","UND","UND","UND","UND","UND"))

nodes$group <- as.factor(c("yes","1","1","1","1","1","1","1","1","1","1","2","3","4","5"))
nodes <- nodes[-16,]

p <- sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", 
                   sinksRight=TRUE, nodeWidth = 40, nodePadding = 5, width = 500, height = 500,
                   fontSize = 0, fontFamily = "Arial", colourScale = my_color2,
                   LinkGroup = "group", NodeGroup = "group", iterations = 0)
p


webshot("C:/Users/tomca/Pictures/test.html", "output.pdf",  zoom = 0.5)

p2 <- sankeyNetwork(Links = links_abdab, Nodes = nodes_abdab,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", 
                   sinksRight=TRUE, nodeWidth = , nodePadding = 5,  width = 500, height = 500,
                   fontSize = 20, fontFamily = "Arial", colourScale = my_color,  
                   LinkGroup = "source", )

p2

add.alpha(brewer.pal(6,"Purples"), .8)



ggplot(full, aes(x=donor, y = 1-as.numeric(mut.frac.V.heavy_mean), color = as.numeric(reads_heavy_chain))) +
  geom_jitter(shape = 16, size = 2, width = 0.25) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.text = element_text(size = 10, color = "black"),
        axis.line = element_line(size = .5, color = "black"),
        legend.position = "none",
        aspect.ratio = 1/2) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0.85, 1, 0.05), labels = c("85", "90", "95", "100"), limits = c(0.85, 1.01)) +
  scale_color_gradient(high = "purple", low = "grey80")


ggplot(full, aes(x=as.numeric(CD27_counts))) +
    geom_density(aes(y=..count..)) + scale_x_continuous(limits = c(0,4000))



#####


scim$epitope <- case_when(str_count(scim$Antigenic.Site, "RBD") == TRUE ~ "RBD",
                          str_count(scim$Antigenic.Site, "S2") == TRUE ~ "S2",
                          str_count(scim$Antigenic.Site, "NTD") == TRUE ~ "NTD")


test <- dcast(scim, epitope ~ paste(VH.Germline, VL.Germline ,sep="/"), fun.aggregate = length)
test <- test[-4,]
test2 <- test[,-1]
rownames(test2) <- test[,1]
test3 <- t(test2)

write.table(test3, "clipboard", sep = "\t")
iris.pca <- PCA(test3, graph = FALSE)


  fviz_pca_biplot(iris.pca, repel = TRUE,
                col.var = "blue", # Variables color
                ggtheme = theme_minimal())

  res.hcpc <- HCPC(iris.pca, graph = FALSE)
  data <- res.hcpc$data.clust
  molt <- melt(data)

 que <- molt %>%
    group_by(clust, variable) %>%
    summarise(value = sum(value))
  
fviz_dend(res.hcpc, 
            cex = 0.9,                     # Label size
            palette = "nejm",               # Color palette see ?ggpubr::ggpar
            rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
            rect_border = "nejm",           # Rectangle color
            labels_track_height = 0.8      # Augment the room for labels
  )
  
  
p <- fviz_cluster(res.hcpc,
               show.clust.cent = FALSE,
               repel = FALSE,            # Avoid label overlapping
               palette = "nejm",         # Color palette see ?ggpubr::ggpar
               ggtheme = theme_minimal(),
  )
 
p
data <- p$data # this is all you need
hull_data <-  data %>%
  group_by(cluster) %>%
  slice(chull(x, y))

ggplot(data, aes(x, y, color = cluster)) + 
  geom_point(size = 3, position = position_jitter(h = 0.1, w = 0.1)) +
  geom_polygon(data = hull_data, alpha = 0.5, aes(fill=cluster)) + 
  scale_fill_nejm() + scale_color_nejm() +
  geom_text_repel(aes(label=ifelse(cluster==2,as.character(name),'')),hjust=0, vjust=0, nudge_x = 2.5, nudge_y = -.5) +
  geom_text_repel(aes(label=ifelse(cluster==4,as.character(name),'')),hjust=0, vjust=0) +
  theme(panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent", color = NA), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 20, color = "black"),
        axis.text.y = element_text(size = 20, color = "black"),
        axis.title.x = element_text(size = 20, face = "plain", color = "black"),
        axis.title.y = element_text(size = 20, face = "plain", color = "black"),
        axis.ticks.length = unit(.2,"cm"),
        legend.text = element_text(size = 20),
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_blank(),
        legend.title = element_text(size = 20, color = "black"),
        plot.title = element_text(size = 20, hjust = .5),
        aspect.ratio = 1) +
  xlab("Dim1 (37.0%)") + ylab("Dim2 (33.8%)") +


ggplot(que, aes(x=clust, y = value, fill = variable)) +
  geom_col(position = "stack", color = "black") +
  scale_fill_simpsons() +
  theme(panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent", color = NA), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 20, color = "black"),
        axis.text.y = element_text(size = 20, color = "black"),
        axis.title.x = element_text(size = 20, face = "plain", color = "black"),
        axis.title.y = element_text(size = 20, face = "plain", color = "black"),
        axis.ticks.length = unit(.2,"cm"),
        legend.text = element_text(size = 20),
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_blank(),
        legend.title = element_text(size = 20, color = "black"),
        plot.title = element_text(size = 20, hjust = .5),
        aspect.ratio = 1) +
  scale_y_continuous(expand = c(0,0)) + 
  xlab("Cluster") + ylab("Number of B cells")






#### EPITOPE SPECIFICITY IN SCI IMMUNOL ####
scim$epitope <- case_when(str_count(scim$Antigenic.Site, "RBD") == TRUE ~ "RBD",
                          str_count(scim$Antigenic.Site, "S2") == TRUE ~ "S2",
                          str_count(scim$Antigenic.Site, "NTD") == TRUE ~ "NTD",
                          str_count(scim$Antigenic.Site, "Und") == TRUE ~ "Undefined")


pairs_epitope <- dcast(scim, epitope ~ paste(VH.Germline, VL.Germline ,sep="/"), fun.aggregate = length)
pairs_epitope <- pairs_epitope[-5,]
pairs_epitope2 <- pairs_epitope[,-1]
rownames(pairs_epitope2) <- pairs_epitope[,1]
pairs_epitope3 <- as.data.frame(t(pairs_epitope2))


pairs_epitope3$sum <- rowSums(pairs_epitope3)
pairings12 <- pairs_epitope3 %>% slice_max(sum, n = 10)
pairings12 <- pairings12[, -ncol(pairings12)]
pairings12$pair <- rownames(pairings12)
p <- melt(pairings12)

epitope_pairing <- ggplot(p, aes(x=reorder(pair, value, function(x){ sum(x) }), y = as.numeric(value), fill = variable)) +
  geom_col(position = "stack", color = "black") +
  theme(panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent", color = NA), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(size = 12, face = "plain", color = "black"),
        axis.title.y = element_blank(),
        axis.ticks.length = unit(.2,"cm"),
        legend.text = element_text(size = 12),
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_blank(),
        legend.title = element_text(size = 12, color = "black"),
        legend.position = c(0.7, 0.2),
        plot.title = element_text(size = 12, hjust = .5),
        aspect.ratio = 1.5)  + coord_flip() + 
        scale_fill_manual(values = color42) +
        scale_y_continuous(expand = c(0,0)) +
      ylab("\nNumber of pairs")
  
epitope_pairing
scimund <- scim[scim$Antigenic.Site == "Undefined",]
und <- dcast(scimund, Antigenic.Site ~ paste(VH.Germline, VL.Germline ,sep="/"), fun.aggregate = length)
tund <- as.data.frame(t(und))

donor_pheno_plot | iso_pheno_plot / (vh_usage | vl_usage | epitope_pairing)

p10 <- melt(pairings12)

p10_new = p10 %>% group_by(variable, pair) %>% tally()

cp <- coord_polar(theta = "y")
cp$is_free <- function() TRUE

ggplot(data=p10, aes(x=2, y=value, group=variable, fill=variable)) +
  geom_bar(width = 1, stat = "identity", color = "black") +
  xlim(c(1,3)) +
  cp +
  facet_wrap(~pair) + theme_void() + theme(text = element_text(size = 12), aspect.ratio = 1) 


