library("dada2")
library("seqinr")
library("biomformat")
library("data.table")
library("ggplot2")
library("ggthemes")
library("ampvis2")
library("cowplot")
library("phyloseq")

setwd("~/acer_microbiome_2020")

# phyloseq object output from Dada2
ps <- readRDS("ps_object_sp_acer.rds")
seqtab <- readRDS("seqtab_acer.rds")
# sequences object made from Dada2

# exporting sequences to a fasta file for import into qiime
uniqueSeqs <- as.list(colnames(seqtab))
write.fasta(uniqueSeqs, uniqueSeqs, "uniqueSeqs.fasta")

# import metadata and merge into phyloseq object
mapfile = "Acer_microbiome_metadata_NR.txt"
map = import_qiime_sample_data(mapfile)
sample_data(ps) <- map

# export .tsv asv_table to import into qiime2
otu<-t(as(otu_table(ps),"matrix"))
otu_biom<-make_biom(data=otu)
write_biom(otu_biom,"ps.biom")

# export taxonomy to import into qiime2
tax<-as(tax_table(ps),"matrix")
tax_cols <- c("Kingdom", "Phylum", "Class", "Order","Family","Genus", "Species")
tax<-as.data.frame(tax)
tax$taxonomy<-do.call(paste, c(tax[tax_cols], sep=";"))
for(co in tax_cols) tax[co]<-NULL
write.table(tax, "taxonomy.txt", quote=FALSE, col.names=FALSE, sep="\t")

# summary of data
ps
summary(sample_data(ps))

ntaxa(ps)
nsamples(ps)
rank_names(ps)
sample_names(ps)[1:5]
sample_variables(ps)

# remove mitochondria and chloroplasts, is.na important becuase if not included
# this command will also remove all Family = NA or Order = NA
ps_with_mito = subset_taxa(ps, (Order!="Chloroplast") | is.na(Order))
ps_no_mito = subset_taxa(ps_with_mito, (Family!="Mitochondria") | is.na(Family))
ps_no_Eukaryota = subset_taxa(ps_no_mito, (Kingdom!="Eukaryota") | is.na(Kingdom))

ps = ps_no_Eukaryota #lost 71 taxa

#sara's pruning step
pstax_filt<-filter_taxa(ps,function(x) sum(x>10) > (0.01*length(x)),TRUE)

summary_tab <- data.frame(init=(as.matrix(sample_sums(ps_full)))[,1],
                           chloros_removed=(as.matrix(sample_sums(ps_with_mito)))[,1],
                           mitos_removed=(as.matrix(sample_sums(ps_no_mito)))[,1],
                           euks_removed=(as.matrix(sample_sums(ps_no_Eukaryota)))[,1],
                           pruned=(as.matrix(sample_sums(pstax_filt)))[,1])
write.table(summary_tab, "reads_lost_phyloseq.txt", quote=FALSE, col.names=FALSE, sep="\t")

# pst = fast_melt(ps)
# prevdt = pst[, list(Prevalence = sum(count > 0), 
#                     TotalCounts = sum(count)),
#              by = taxaID]
# keepTaxa = prevdt[(Prevalence >=0 & TotalCounts >15), taxaID]
# ps_pruned = prune_taxa(keepTaxa,ps)
# ps_pruned
# sample_sums(ps_pruned)
# min(sample_sums(ps_pruned))
# ps <- ps_pruned
# save(ps, file = "ps_pruned.RData")

#plot otus by sequencing depth
observed <- estimate_richness(ps, measures = c('Observed'))
explore.df <- cbind(observed, sample_sums(ps), sample_data(ps)$Location)
colnames(explore.df) <- c('Observed', 'Sample_Sums', 'Location')
observed_mean <- mean(explore.df$Observed)
sample_sum_mean <- mean(explore.df$Sample_Sums)
ggplot(data = explore.df, aes(x = Sample_Sums, y = Observed, color = Location)) + 
  geom_point() +
  geom_smooth(method="auto", se=TRUE, fullrange=FALSE, level=0.95, 
              inherit.aes = F, mapping = aes(Sample_Sums, Observed),
              data = explore.df) +
  ylab("Observed OTUs") +
  scale_colour_colorblind()

# Let's use ampvis2 again so we can easily make a rarefaction curve

# Need to convert from phyloseq to ampvis
av2_otutable <- data.frame(OTU = rownames(t(phyloseq::otu_table(ps)@.Data)),
                           t(phyloseq::otu_table(ps)@.Data),
                           phyloseq::tax_table(ps)@.Data,
                           check.names = F
)

#Extract metadata from the phyloseq object:
av2_metadata <- data.frame(phyloseq::sample_data(ps), 
                           check.names = F
)

av2_metadata <- cbind(rownames(av2_metadata), av2_metadata)

#Load the data with amp_load:
av2_obj <- amp_load(av2_otutable, av2_metadata)

# RARE CURVE
rare_plot_amp <- amp_rarecurve(data = av2_obj, color_by = "Region")
rare_curve_plot <- rare_plot_amp + ylab('Observed ASVs (count)') + 
  geom_vline(xintercept=min(sample_sums(ps)), linetype='dashed') +
  scale_colour_colorblind() +
  xlim(c(0, 35000))
plot(rare_plot_amp)
rare_curve_plot


#If we rarefy to 26537
rare_plot_amp <- amp_rarecurve(data = av2_obj, color_by = "Location")
rare_curve_plot <- rare_plot_amp + ylab('Observed ASVs (count)') + 
  geom_vline(xintercept=c(26537), linetype='dashed') +
  scale_colour_colorblind() +
  xlim(c(0, 35000))
#plot(rare_plot_amp)
rare_curve_plot

summary(explore.df$Sample_Sums)
length(explore.df$Sample_Sums)

# plot alpha diversity
plot_richness(ps, x = "Location", color = "Geno", measures = c('Shannon', 'Simpson'))

ps_rarefied <- rarefy_even_depth(ps, sample.size = 26537, rngseed = 999) #546 OTUs removed
ps_rarefied
sum(sample_sums(ps_rarefied))
sample_sums(ps_rarefied)

# # plot alpha diversity
richness.rare <- cbind(estimate_richness(ps_rarefied, 
                                          measures = c('Shannon', 'Simpson')),
                        sample_data(ps_rarefied)$Location)
colnames(richness.rare) <- c('Shannon', 'Simpson', 'Location')
richness.rare$Labels <- rownames(richness.rare)
# 
ggplot(data = richness.rare, aes(x = Shannon, y = Simpson, color = Location)) + 
   geom_point()

ggplot(data = richness.rare, aes(x = Shannon, y = (Simpson-1)*-1, color = Location)) + 
   geom_point() +
     geom_smooth(method="auto", se=TRUE, fullrange=FALSE, level=0.95, 
              inherit.aes = F, mapping = aes(Shannon, Simpson),
                 data = richness.rare) +
   scale_color_colorblind() +
   geom_text(aes(label=ifelse(Shannon<1, Labels, ""), hjust=-0.1),
             show.legend = F) +
   xlim(c(0,max(richness.rare$Shannon))) +
   theme_cowplot()

#relative abundance transform
ps_rel = filter_taxa(ps_rarefied, function(x) mean(x) > 0.1, TRUE)
ps_rel = transform_sample_counts(ps_rel, function(x) x / sum(x) )
ps_rel



# save ps object
save(ps, file = "ps_unrarefied.RData")
save(ps_rarefied, file = "ps_rare_acer.RData")
save(ps_rel, file = "ps_rel.RData")
