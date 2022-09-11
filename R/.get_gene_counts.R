.get_gene_counts <- function(df.DEG, df.DEG.ref.matched){

  #save counts of up and down regulated genes before accounting for whether their gene symbol can be matched to a TOA diagnosticity score
  up_reg_count <- table(df.DEG$DEG)["1"]
  down_reg_count <- table(df.DEG$DEG)["-1"]
  reg_count <- down_reg_count + up_reg_count
  gene_count <- rep(c(up_reg_count, down_reg_count, reg_count),2)
  names(gene_count) <- rep(c("upregulated", "downregulated", "regulated"),2)

  #save counts of up and down regulated genes after accounting for whether their gene symbol can be matched to a TOA diagnosticity score
  used_up_reg_count <- table(df.DEG.ref.matched$DEG)["1"]
  used_down_reg_count <- table(df.DEG.ref.matched$DEG)["-1"]
  used_reg_count <- used_down_reg_count + used_up_reg_count
  used_gene_count <- rep(c(used_up_reg_count, used_down_reg_count, used_reg_count),2)
  names(used_gene_count) <- rep(c("upregulated", "downregulated", "regulated"),2)

  return(list(gene_count = gene_count,
              used_gene_count = used_gene_count))
}
