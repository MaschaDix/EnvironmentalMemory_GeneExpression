# GO annotations
# see https://github.com/lyijin/topGO_pipeline for more detail and folder structure

library(topGO)
library(ggplot2)

setwd()

folders <- list.files("./genes_of_interest")
folders

for (i in 1:length(folders)) {
  folder_of_interest = file.path("./genes_of_interest", folders[i],"")
  mult_files = list.files(folder_of_interest, pattern="*.txt")
  for (go_category in c('bp','cc','mf')) {
    annot_filename='./aip_prot_no_duplication.tsv'
    gene_id_to_go = readMappings(file=annot_filename)
    gene_id_to_go = gene_id_to_go[gene_id_to_go != 'no_hit']
    gene_names = names(gene_id_to_go)
    
    for (m in mult_files) {
      print(paste("Current file:", m))
      genes_of_interest_filename = paste0(folder_of_interest, m)
      genes_of_interest = scan(genes_of_interest_filename, character(0), sep="\n")
      
      genelist = factor(as.integer(gene_names %in% genes_of_interest))
      names(genelist) = gene_names
      
      GOdata = try(new("topGOdata", ontology=toupper(go_category), allGenes=genelist, gene2GO=gene_id_to_go, annotationFun=annFUN.gene2GO))
      
      # handle error
      if (class(GOdata) == "try-error") {
        print (paste0("Error for file", m, "!"))
        next
      }
      
      # weight01 is the default algorithm used in Alexa et al. (2006)
      weight01.fisher <- runTest(GOdata, statistic = "fisher")
      
      # generate a results table (for only the top 1000 GO terms)
      #   topNodes: highest 1000 GO terms shown
      #   numChar: truncates GO term descriptions at 1000 chars (basically, disables truncation)
      results_table = GenTable(GOdata, P_value=weight01.fisher, orderBy="P_value", topNodes=1000, numChar=1000)
      
      # write it out into a file for python post-processing
      output_filename = paste0("./topGO_output/", folders[i], "/", go_category, "_", m)
      write.table(results_table, file=output_filename, quote=FALSE, sep='\t')
    }
  }
}



# make dot plots for top 25 up- and down-regulated genes in MEturquoise

results_table <- read.table("MEturquoise_UPDOWN_topGOs.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote="", fill=TRUE)

# Convert P_value to numeric if it isn't
str(results_table)
results_table$P_value <- as.numeric(results_table$P_value)

# Convert terms into factor with levels in the current order so ggplot doesn't put them in alphabetically
results_table$Term <- factor(results_table$Term, levels = rev(results_table$Term))

term_colors <- ifelse(results_table$Direction == "UP", "#C5081A", "#00008B")
names(term_colors) <- results_table$Term

# Generate dot plot using ggplot
pdf("MEturquoise_UPDOWN_topGOs.pdf", width = 13, height = 10)
ggplot(results_table, aes(x = Term, 
                          y = as.numeric(Significant/Annotated), 
                          size = Annotated, 
                          color = -log10(as.numeric(P_value)))) +
  geom_point() +
  scale_color_viridis_c(option = "plasma") +
  scale_x_discrete(labels = function(x) {
    sapply(x, function(term) {
      paste0("<span style='color:", term_colors[term], "'>", term, "</span>")
    })
  }) +
  theme_bw() +
  coord_flip() +
  labs(title = "MEturquoise All Up", 
       y = "Significant/Annotated ratio", x = "GO Term",
       size = "Annotated", color = "-log10 of p-value") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.y = ggtext::element_markdown(size = 14),
    axis.title = element_text(size = 16, color = "#202020")
  )
dev.off()
