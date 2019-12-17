library(dplyr)
library(circlize)
library(reshape2)
library(tidyr)
library(data.table)
library(ComplexHeatmap)
options(stringsAsFactors = F)
# ht_opt(ROW_ANNO_PADDING = unit(0.5, "cm"))

snv_data <- function(snv_file){
  snv_info <- read.table(snv_file, header = T, sep = "\t")
  BPT_snv_info <- snv_info[substr(snv_info$sample, 1, 3) == "BPT", ] %>%
    separate_rows(Gene.refGene, sep = ";")
  BPT_snv_info <- BPT_snv_info[(BPT_snv_info$sample != "BPT085"),]
  BPT_snv_info$ExonicFunc.refGene[BPT_snv_info$ExonicFunc.refGene == "nonsynonymous SNV"] <- "missense"
  BPT_snv_info$ExonicFunc.refGene[BPT_snv_info$ExonicFunc.refGene == "nonframeshift deletion"] <- "infradel"
  BPT_snv_info$ExonicFunc.refGene[BPT_snv_info$ExonicFunc.refGene == "frameshift insertion"] <- "frameins"
  BPT_snv_info$ExonicFunc.refGene[BPT_snv_info$ExonicFunc.refGene == "frameshift deletion"] <- "framedel"
  BPT_snv_info$ExonicFunc.refGene[BPT_snv_info$ExonicFunc.refGene == "nonframeshift insertion"] <- "infrains"
  chosen_sample <- unique(BPT_snv_info$sample)
  snv_data_list <- list("BPT_snv_info" = BPT_snv_info, 
                        "sample" = chosen_sample)
  return(snv_data_list)
}

### get the matrix to do main oncoprint
snv_plot_data <- function(snv_plot, snv_genes){
  snv_chosen_gene <- read.table(snv_genes, header = F, sep = "\t")[,1] %>% as.character
  sample_info <- aggregate(ExonicFunc.refGene ~ Gene.refGene + sample, 
                           data = snv_plot, FUN = function(x) paste0(unique(x), collapse = ";"))
  sample_info$ExonicFunc.refGene[unlist(lapply(sample_info$ExonicFunc.refGene, 
                                               function(x) 
                                                 length(strsplit(x, split = ";")[[1]]))) > 1] <- "multiple"
  heatmap_info <- dcast(sample_info, Gene.refGene ~ sample,
                        value.var = "ExonicFunc.refGene", drop = T)
  rownames(heatmap_info) <- heatmap_info[,1]
  heatmap_matrix <- as.matrix(heatmap_info[,-1])
  heatmap_matrix <- heatmap_matrix[rownames(heatmap_matrix) %in% snv_chosen_gene,]
  return(heatmap_matrix)
}

##################
input_info <- function(file_name){
  input_context <- read.table(file_name, row.names = 1, header = F, sep = "\t")
  return(input_context)
}

### deal with top common annotation
top_common_annotation <- function(treat_file, subtype_file, metas_file, chosen_sample){
  treat_info <- input_info(treat_file)
  subtype_info <- input_info(subtype_file)
  metas_info <- input_info(metas_file)
  sample_order <- rownames(subtype_info)[rownames(subtype_info) %in% chosen_sample]
  treat_order <- treat_info[sample_order,]
  metas_order <- metas_info[sample_order,]
  subtype_order <- subtype_info[sample_order,]
  top_list <- list("treat_order" = treat_order, 
                   "subtype_order" = subtype_order, 
                   "metas_order" = metas_order, 
                   "sample_order" = sample_order)
  return(top_list)
}

### get all top annotation data format
snv_top_annotation <- function(snv_plot, treat_info, subtype_info, metas_info){
  mutation_info <- plyr::count(snv_plot, vars = c("sample", "ExonicFunc.refGene")) %>%
    acast(sample ~ ExonicFunc.refGene, value.var = "freq", fill = 0)
  ha = HeatmapAnnotation(cbar = anno_barplot(mutation_info, 
                                             gp = gpar(fill = c("#FFC125", "#BF3EFF", "red", "#00FFFF", "blue", "#008000", "#FF4040"),
                                                       col = c("#FFC125", "#BF3EFF", "red", "#00FFFF", "blue", "#008000", "#FF4040")),
                                             border = F),
                         subtype = subtype_info,
                         metastasis = metas_info,
                         treat = treat_info,
                         col = list(subtype = c("Basal" = "slategray1", "Her 2+" = "#008000",
                                                "Luminal A" = "#FF4040", "Luminal B" = "#9F79EE",
                                                "Unknown" = "steelblue"),
                                    metastasis = c("Primary" = "red", 
                                                   "Metastasis" = "darkblue",
                                                   "Relapse" = "orange"),
                                    treat = c("Treated" = "blue", "Untreated" = "#FFC125")
                         ),
                         annotation_name_side = "left",
                         annotation_name_gp = gpar(fontsize = 8, fontface = "bold"),
                         show_annotation_name = c(cbar = FALSE)
  )
  return(ha)
}

gene_frequency_statistics <- function(total_info, samples, chosen_genes){
  sample_info <- total_info[(total_info$sample %in% samples),]
  dealt_info <- sample_info %>%
    separate_rows(Gene.refGene, sep = ";") %>%
    select(sample, Gene.refGene) %>%
    group_by(Gene.refGene) %>%
    dplyr::summarise(val=paste(unique(sample), collapse=";")) %>%
    as.data.frame
  
  dealt_info$count <- sapply(dealt_info$val,
                             function(x) length(strsplit(as.character(x), split = ";")[[1]]))
  annot_gene_frequency <- dealt_info[(dealt_info$Gene.refGene %in% chosen_genes),]
  annot_gene_frequency <- merge(dealt_info, as.data.frame(chosen_genes),
                                by.x = "Gene.refGene",
                                by.y = "chosen_genes",
                                all.y = T)
  annot_gene_frequency[is.na(annot_gene_frequency)] <- 0
  annot_gene_frequency$frequency <- (annot_gene_frequency$count)/(length(unique(sample_info$sample)))
  rownames(annot_gene_frequency) <- annot_gene_frequency[,1]
  annot_gene_frequency <- annot_gene_frequency[,-1]
  return(annot_gene_frequency)
}

left_annotation_data <- function(total_info, chosen_genes_file, treat_file, heatmap_matrix){
  treat_info <- input_info(treat_file)
  chosen_genes <- read.table(chosen_genes_file, header = F, sep = "\t")[,1] %>% as.character
  treated_samples <- rownames(treat_info)[treat_info[,1] == "Treated"]
  untreated_samples <- rownames(treat_info)[treat_info[,1] == "Untreated"]
  treated_frequency <- gene_frequency_statistics(total_info, treated_samples, chosen_genes)
  untreated_frequency <- gene_frequency_statistics(total_info, untreated_samples, chosen_genes)
  treated_annot <- treated_frequency[rownames(heatmap_matrix),3]
  untreated_annot <- untreated_frequency[rownames(heatmap_matrix),3]
  max_value = round(max(treated_annot, untreated_annot), digits = 1)
  max_lable <- paste((100*max_value), "%", sep = "")
  ra = rowAnnotation(treat = anno_barplot(as.matrix(treated_annot), gp = gpar(fill = "blue", col = "blue"),
                                          border = F,
                                          ylim = c(0, max_value),
                                          width = unit(1, "cm"),
                                          axis_param = list(direction = "reverse",
                                                            side = "top",
                                                            at = c(0, max_value),
                                                            labels = c("0", max_lable))),
                     untreat = anno_barplot(as.matrix(untreated_annot), gp = gpar(fill = "#FFC125", col = "#FFC125"),
                                            border = F,
                                            ylim = c(0, max_value),
                                            width = unit(1, "cm"),
                                            axis_param = list(side = "top",
                                                              at = max_value,
                                                              labels = max_lable)),
                     
                     show_annotation_name = F,
                     annotation_name_offset = unit(5, "mm")
  )
  return(ra)
}

snv_oncoprint <- function(heatmap_matrix, sample_order, dirver_genes, ha = NULL, ra = NULL){
  heatmap_matrix <- heatmap_matrix[, sample_order]
  alter_fun_list = list(
    background = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
                gp = gpar(fill = "#CCCCCC", col = NA))
    },
    stopgain = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
                gp = gpar(fill = "#008000", col = NA))
    },
    stoploss = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
                gp = gpar(fill = "#FF4040", col = NA))
    },
    missense = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
                gp = gpar(fill = "blue", col = NA))
    },
    infradel = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
                gp = gpar(fill = "red", col = NA))
    },
    infrains = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
                gp = gpar(fill = "#00FFFF", col = NA))
    },
    framedel = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
                gp = gpar(fill = "#FFC125", col = NA))
    },
    frameins = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
                gp = gpar(fill = "#BF3EFF", col = NA))
    },
    multiple = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
                gp = gpar(fill = "black", col = NA)) 
    }
  )
  col = c("stopgain" = "#008000", "stoploss" = "#FF4040", "missense" = "blue",
          "infradel" = "red", "infrains" = "#00FFFF", "multiple" = "black",
          "framedel" = "#FFC125", "frameins" = "#BF3EFF") 
  ht1 = oncoPrint(heatmap_matrix, get_type = function(x) strsplit(x, ";")[[1]],
                  alter_fun = alter_fun_list, col = col, show_column_names = T,
                  show_pct = F, row_names_side = "left",
                  pct_side = "right",
                  column_order = sample_order, 
                  top_annotation = ha,
                  right_annotation = ra,
                  row_names_gp = gpar(fontsize = 8, fontface = "bold"),
                  column_names_gp = gpar(fontsize = 8, fontface = "bold"),
                  remove_empty_columns = F, 
                  heatmap_legend_param = list(title = "Alternations", 
                                              at = c("stopgain", "missense", "infradel", "infrains", "framedel", "frameins", "multiple"),
                                              labels = c("Stop gain", "Missense", "Inframe deletion", "Inframe insertion", 
                                                         "Frameshift deletion", "Frameshift insertion", "Multiple mutation")),
                  row_split = ifelse(rownames(heatmap_matrix) %in% dirver_genes, "b", "a"),
                  row_gap = unit(2, "mm"),
                  row_title = NULL)
  return(ht1)
}

cnv_data <- function(cnv_file, cnv_chosen_genes){
  cnv_info <- read.table(cnv_file, header = F, sep = "\t")
  genelist <- read.table(cnv_chosen_genes, header = F, sep = "\t")[,1] %>% as.character
  colnames(cnv_info) <- c("sample", "Chr", "Start", "End", "Ref", "Alt", "Func.refGene",
                          "Gene.refGene", "GeneDetail.refGene", "ExonicFunc.refGene",
                          "AAChange.refGene", "cytoBand", "ratio", "cn")
  BPT_cnv <- cnv_info[substr(cnv_info$sample, 1, 3) == "BPT",]
  BPT_cnv <- BPT_cnv[(BPT_cnv$sample != "BPT085"),]
  BPT_cnv$ratio <- BPT_cnv$ratio * 2
  chosen_sample <- unique(BPT_cnv$sample)
  cnv_pre_info <- BPT_cnv %>%
    separate_rows(Gene.refGene, sep = ";") %>%
    select(sample, Gene.refGene, ratio)
  cnv_total_info <- aggregate(cnv_pre_info$ratio, list(cnv_pre_info$sample,
                                                       cnv_pre_info$Gene.refGene), mean)
  colnames(cnv_total_info) <- c("sample", "Gene.refGene", "Ratio")
  heatmap_matrix <- acast(cnv_total_info, Gene.refGene ~ sample, 
                          value.var = "Ratio")
  heatmap_matrix[is.na(heatmap_matrix)] <- 2
  heatmap_matrix[heatmap_matrix > 4] <- 4
  cnv_heatmap <- heatmap_matrix[rownames(heatmap_matrix) %in% genelist,]
  cnv_data_list <- list("BPT_cnv_info" = BPT_cnv,
                        "cnv_matrix" = cnv_heatmap, 
                        "chosen_sample" = chosen_sample)
  return(cnv_data_list)
}

cnv_top_annotations <- function(treat_info, subtype_info, metas_info){
  cnv_ta <- HeatmapAnnotation(treat = treat_info,
                              metastasis = metas_info, 
                              subtype = subtype_info,
                              col = list(treat = c("Treated" = "blue", "Untreated" = "#FFC125"),
                                         metastasis = c("Primary" = "red", 
                                                        "Metastasis" = "darkblue",
                                                        "Relapse" = "orange"),
                                         subtype = c("Basal" = "slategray1", "Her 2+" = "#008000",
                                                     "Luminal A" = "#FF4040", "Luminal B" = "#9F79EE",
                                                     "Unknown" = "steelblue")),
                              annotation_name_side = "left")
  return(cnv_ta)
}

cnv_Heatmap <- function(cnv_heatmap, sample_order, driver_genes, cnv_ra = NULL, cnv_ta = NULL){
  col_fun = colorRamp2(c(0, 2, 4), c("darkblue", "white", "red"))
  cnv_heatmap_matrix <- cnv_heatmap[, sample_order]
  ht2 = Heatmap(as.matrix(cnv_heatmap_matrix), col = col_fun, name = "Copy Number", 
                row_names_side = "left", column_names_side = "bottom",
                column_order = sample_order,
                top_annotation = cnv_ta,
                right_annotation = cnv_ra,
                cluster_rows = F, cluster_columns = F,
                row_names_gp = gpar(fontsize = 8, fontface = "bold"),
                column_names_gp = gpar(fontsize = 8, fontface = "bold"),
                heatmap_legend_param = list(
                  direction = "vertical",
                  title = "Copy Number", at = c(0,1,2,3,4),
                  labels = c("0", "1", "2", "3", "4 or more")          
                ),
                row_split = ifelse(rownames(cnv_heatmap_matrix) %in% driver_genes, "b", "a"),
                row_gap = unit(2, "mm"),
                row_title = NULL)
  return(ht2)
}




snv_file = "~/Desktop/Ping/snv/merge.txt"
snv_chosen_genes = "~/Desktop/Ping/snv/genelist.txt"
snv_driver_genes = "~/Desktop/Ping/snv/drivergenes.txt"
cnv_file = "~/Desktop/Ping/Control-Freec/cnv_exonic.txt"
cnv_chosen_genes = "~/Desktop/Ping/cnv/genelist.txt"
cnv_driver_genes = "~/Desktop/Ping/cnv/drivergenes.txt"
treat_file = "~/Desktop/Ping/snv/treat.txt"
subtype_file <- "~/Desktop/Ping/snv/subtype.txt"
metas_file <- "~/Desktop/Ping/snv/meta.txt"


snv_data_res <- snv_data(snv_file)
snv_heatmap_matrix <- snv_plot_data(snv_data_res$BPT_snv_info, snv_chosen_genes)
common_annot_list <- top_common_annotation(treat_file, subtype_file, metas_file, snv_data_res$sample)
snv_ha <- snv_top_annotation(snv_data_res$BPT_snv_info, 
                             common_annot_list$treat_order, 
                             common_annot_list$subtype_order, 
                             common_annot_list$metas_order)
snv_ra <- left_annotation_data(snv_data_res$BPT_snv_info, snv_chosen_genes, treat_file, snv_heatmap_matrix)  
snv_driver <- read.table(snv_driver_genes, header = F, sep = "\t")[,1] %>% as.character
snv_plot <- snv_oncoprint(snv_heatmap_matrix, 
                          common_annot_list$sample_order,
                          snv_driver,
                          snv_ha,
                          snv_ra)

cnv_data_res <- cnv_data(cnv_file, cnv_chosen_genes)
cnv_ta <- cnv_top_annotations(common_annot_list$treat_order, 
                              common_annot_list$subtype_order, 
                              common_annot_list$metas_order)
cnv_rannot <- left_annotation_data(cnv_data_res$BPT_cnv_info, cnv_chosen_genes, 
                                   treat_file, cnv_data_res$cnv_matrix)
cnv_driver <- read.table(cnv_driver_genes, header = F, sep = "\t")[,1] %>% as.character
cnv_plot <- cnv_Heatmap(cnv_data_res$cnv_matrix, common_annot_list$sample_order, 
                        cnv_driver,
                        cnv_ra = cnv_rannot)
ht_list = snv_plot %v% cnv_plot
draw(ht_list, ht_gap = unit(c(10, 10), "mm"), 
     heatmap_legend_side = "right",
     annotation_legend_side = "right")
