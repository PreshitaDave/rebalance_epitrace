

library(dplyr)
library(patchwork)
library(ggpubr)
library(ggrepel) 
library(ggtree)
library(ggplot2)
library(ape)

# Figure 6c
da_tf_with_activity_2 <- readRDS("./Packed_Figure_Reproducibility/Source_Data_for_Figure/Figure6/Figure6_SHARESEQ_brain_multiome_NR2F1_cluster_motif_differential_activity_result.rds")

ggplot(da_tf_with_activity_2[da_tf_with_activity_2$motif_origin %in% 'human',],aes(x=-log(p_val_adj.x,10),y=-log(p_val_adj.y,10))) + geom_point() + ggrepel::geom_label_repel(data=da_tf_with_activity_2[-log(da_tf_with_activity_2$p_val_adj.x,10)>3 & -log(da_tf_with_activity_2$p_val_adj.y,10)>100,],mapping=aes(label=gene_symbol)) + theme_classic() + xlab('-logP Gene Exp') + ylab('-logP TF activity')+ theme(text=element_text(size=10)) + ggtitle('Statistics for TF activity and Gene Exp') -> p1

# Figure 6d
readRDS("./Packed_Figure_Reproducibility/Source_Data_for_Figure/Figure6/Figure6_Brain_Multiome_NR2F1_TFBS_on_UMAP_plot_df.rds") -> NR2F1_tfbs_plot

ggplot(NR2F1_tfbs_plot,aes(x=UMAP_1,y=UMAP_2,fill=Motif_Activity))+geom_point(size=2.2,pch=21,color='darkgray')+xlim(c(-10.2,4.5)) + ylim(c(-9,4)) + scale_fill_gradientn(colors=c('cornflowerblue','beige','red')) + theme_void() + ggtitle(paste0(da_tf_with_activity_2$gene_symbol[da_tf_with_activity_2$motif %in% unique(NR2F1_tfbs_plot$Motif)]%>%unique,' TFBS_ChrAcc ',unique(NR2F1_tfbs_plot$Motif))) + theme(text=element_text(size=10)) -> p2

pdf('./Packed_Figure_Reproducibility/Plots//Figure6c_6d_NR2F1_TCF4_TF_activity.pdf',height=4,width=9)
(p1|p2)
dev.off()

# Figure 6e 6f 
combine_expression_df <- readRDS("./Packed_Figure_Reproducibility/Source_Data_for_Figure/Figure6/Figure6e_6f_Brain_NR2F1_LMO_expression_df.rds")

pdf('./Packed_Figure_Reproducibility/Plots/Figure6e_6f_Brain_NR2F1_LMO3_RNA_Expression.pdf',height=4,width=8)
lapply(unique(combine_expression_df$gene), function(x){
  xdata <- combine_expression_df[combine_expression_df$gene %in% x, ]
  ggplot(xdata,aes(x=UMAP_1,y=UMAP_2,fill=Expression))+geom_point(size=2.2,pch=21,color='darkgray')+xlim(c(-10.2,4.5)) + ylim(c(-9,4)) + scale_fill_gradientn(colors=c('cornflowerblue','beige','red')) + theme_void() + ggtitle(x) + theme(text=element_text(size=20))
}) %>% patchwork::wrap_plots(ncol = 2, guides = 'collect')
dev.off()

# Figure 6g 6h 
nIPC_cytotrace_Epitrace_df <- readRDS('./Packed_Figure_Reproducibility/Source_Data_for_Figure/Figure6//Figure6g_6h_two_population_of_nIPC_Epitrace_and_Cytotrace.rds')

my_comparision <- list(c('LMO3+','NR2F1+'))

ggplot(nIPC_cytotrace_Epitrace_df,aes(x=class,y=EpiTraceAge_iterative)) + geom_violin(aes(fill=class)) + ggbeeswarm::geom_beeswarm(size=3,pch=21) + geom_boxplot(outlier.alpha = 0,width=0.3,fill='gray') +  ggpubr::stat_compare_means(comparisons =my_comparision,label = 'p.signif' ) + theme_classic()  + theme(text=element_text(size=20),axis.title.x = element_blank()) + ylab('EpiTrace Age') + stat_compare_means(label.y = 1.2,size=5) + ggsci::scale_fill_jco() + ggtitle('Age') -> p3

my_comparision <- list(c('LMO3+','NR2F1+'))
ggplot(nIPC_cytotrace_Epitrace_df,aes(x=class,y=cytotrace_rna)) + geom_violin(aes(fill=class)) + ggbeeswarm::geom_beeswarm(size=3,pch=21) + geom_boxplot(outlier.alpha = 0,width=0.3,fill='gray') +  ggpubr::stat_compare_means(comparisons =my_comparision,label = 'p.signif' ) + theme_classic()  + theme(text=element_text(size=20),axis.title.x = element_blank()) + ylab('CytoTRACE by RNA') + stat_compare_means(label.y = 1.2,size=5) + ggsci::scale_fill_jco()  + ggtitle('Stemness') -> p4

pdf('./Packed_Figure_Reproducibility/Plots//Figure6g_6h_compare_Epitrace_and_Cytotrace_two_nIPC_population.pdf',height=4.5,width=9) 
(p3|p4) + plot_layout(guides='collect')
dev.off() 

# Figure 6i 
colorlist <- c('GluN5' = "cadetblue4",
               'IN1'="#A2CD5A",
               'nIPC/GluN1' = 'cornflowerblue',
               'IN2' = "#BCEE68",
               'SP' = "darkgreen",
               'GluN2' = 'cadetblue1',
               'IN3'= "#CAFF70",
               'RG' = "red",
               'GluN4' = 'cadetblue3',
               'GluN3'= 'cadetblue2',
               'Cyc. Prog.' = "#FFA500",
               'mGPC/OPC'= "#68228B",
               'EC/Peric.' = "#8B0000")

brain_meta <- readRDS('./Packed_Figure_Reproducibility/Source_Data_for_Figure/Figure6//Figure6_scMultiome_Brain_metadata.rds')
ggplot(brain_meta,aes(y=celltype,x=EpiTraceAge_iterative)) + geom_violin(scale='width', aes(fill = celltype), alpha = 0.8) + geom_boxplot(width=0.15, outlier.alpha = 0, fill = 'black') + theme_classic() + theme(text=element_text(size=20), legend.position = "none") + scale_fill_manual(values = colorlist) + xlab("\n EpiTrace") + ylab("Cell type") -> pp1
ggplot(brain_meta,aes(y=celltype,x=cytotrace_rna)) + geom_violin(scale='width', aes(fill = celltype), alpha = 0.8) + geom_boxplot(width=0.15, outlier.alpha = 0, fill = 'black') + theme_classic() + theme(text=element_text(size=20), axis.title.y = element_blank(), axis.text.y = element_blank(), legend.position = "none") + scale_fill_manual(values = colorlist) + scale_x_reverse() + xlab("\n CytoTrace") -> pp2

pdf('./Packed_Figure_Reproducibility/Plots//Figure6i_Brain_scMultiome_Epitrace_Cytotrace.pdf',height=9,width=8)
pp1|pp2
dev.off() 

# Figure 6j 
colorset <- c('RG' = "red",'nIPC/GluN1' = 'cornflowerblue','GluN2' = 'cadetblue1','GluN3' = 'cadetblue2','GluN4' = 'cadetblue3','GluN5' = 'cadetblue4','SP' = 'darkgreen')

data.tree_clock <- readRDS("./Packed_Figure_Reproducibility/Source_Data_for_Figure/Figure6//Figure6j_Brain_Epitrace_iterative_tree_result.rds")

data.tree_clock <- ape::root(data.tree_clock,outgroup='RG')
tree_plot_clock <- ggtree::ggtree(data.tree_clock,
                                  layout = "rectangular", ladderize = FALSE) +
  geom_tiplab(aes(), color='black', size = 5, offset = 40) + geom_tippoint(aes(fill=label,color=label),size=7) +
  scale_color_manual(values = colorset)
xmax <- (tree_plot_clock$data$branch.length %>% max(na.rm = T)) *
  1.4
tree_plot_clock <- tree_plot_clock + xlim(c(NA, xmax))

pdf('./Packed_Figure_Reproducibility/Plots//Figure6j_Brain_Epitrace_tree.pdf',width=6,height=3) 
print(tree_plot_clock)
dev.off()
