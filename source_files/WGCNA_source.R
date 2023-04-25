##WGCNA functions

## Check data for missing values
missing <- function(data){
  if(goodSamplesGenes(data, verbose=0)$allOK == F){
    good <- which(goodSamplesGenes(data, verbose=0)$goodGenes == T)
    data[,c(good)]
  } else{data}
}

## Plot sample clustering
sampleclustering <- function(data){
  par(cex=0.7)
  h <- hclust(dist(data), method="average")
  plot(h, main="Sample Clustering for Outliers in log FPKM", sub="", xlab="")
}

## Determining soft threshold
softthreshold <- function(soft){
  
  # Plot correlation with scale free topology
  plot.data <- soft$fitIndices %>%
    dplyr::mutate(R.2=-sign(slope) * SFT.R.sq)
  
  p <- ggplot(plot.data) +
    geom_text(aes(x=Power, y=R.2, label=powers)) +
    geom_hline(yintercept=0.9, lty=2, color="red") +
    labs(x="Soft Threshold Power", y="Signed R-Squared Value") + basic_theme
  
  # Plot mean connectivity
  q <- ggplot(plot.data) +
    geom_text(aes(x=Power, y=mean.k., label=powers)) +
    labs(x="Soft Threshold Power", y="Mean Connectivity") + basic_theme
  
  cat("\n###","Scale Free Topology"," \n")
  print(p)
  cat("\n \n")
  
  cat("\n###", "Mean Connectivity"," \n")
  print(q)
  cat("\n \n")
  
}


## Creation of bar plots for modules
module_barplot <- function(module){
  module.names <- as.character(unique(module$Module))
  name <- str_split(module.names,"_")
  samples <-c()
  for(i in 1:length(name)){
    samples[[i]] <- name[[i]][2]
  }  
  samples <- unlist(samples)
  count <- modules %>%
    count(Module) %>%
    group_by(Module) %>%          # now required with changes to dplyr::count()
    mutate(prop = prop.table(n))
  
  n <- sum(count$prop)
  p <- ggplot(data=count,aes(x=Module,y=n))
  p <- p + geom_bar(color="black", fill=sort(samples), stat="identity",position="identity") + theme_classic()
  p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + scale_x_discrete(labels=sort(samples)) + xlab("Modules")
  p <- p + theme(axis.text = element_text(size = 14)) 
  print(p)
}

##function to create plot that shows relationship between modules and traits

modules_traits <- function(eigens, phenotype){
  name <- str_split(names(eigens),"_")
  samples <-c()
  for(i in 1:length(name)){
    samples[[i]] <- name[[i]][2]
  } 
  nSamples <- nrow(eigens)
  colnames(eigens) <-samples
  
  MEs = orderMEs(eigens)
  moduleTraitCor = cor(MEs, phenotype, use = "p")
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
  
  sizeGrWindow(10,6)
  # Will display correlations and their p-values
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor)
  #par(mar = c(6, 8.5, 3, 3));
  
  # Display the correlation values within a heatmap plot
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = names(phenotype),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = T,
                 colors = colorRampPalette(brewer.pal(12, "Paired"))(10),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.2,
                 zlim = c(-1,1),
                 cex.lab.x = .8,
                 cex.lab.y = .6,
                 main = paste("Module-trait relationships"))
  
}

Data_setup <- function(data, by.x=NULL, by.y=NULL, file="file.RData", folder=NULL){
  ensembl_location<-readRDS(here("Data","Ensembl_gene_id_and_location.RData"))
  WGCNA <- data
  #match(WGCNA$Metabolite,ensembl_location$ensembl_gene_id)
  matched_annotation <- merge(WGCNA,ensembl_location,by.x="Gene",by.y="ensembl_gene_id")
  matched_annotation <- matched_annotation[order(matched_annotation$X..Module.),]
  saveRDS(matched_annotation,here("Data",folder,file))
}

## Function that finds pathways
pathways<-function(annotation, organism=NULL, src_filter=NULL, pathwayfile = "file.RData", genefile = "file.RData", folder = NULL){
  gprofiler2::set_base_url("https://biit.cs.ut.ee/gprofiler_archive3/e104_eg51_p15/")
  pathway=list()
  pathway.genes=list()
  Matched.module <- split(annotation$Gene, annotation$X..Module.)
  for(i in 1:length(Matched.module)){
    name <- str_split(names(Matched.module)[i],"_")
    name <- str_split(name[[1]][2],"\"")
    gene <- as.vector(Matched.module[[i]])
    profiler <- gost(gene,organism = "mmusculus", correction_method = "fdr", sources = c("KEGG","REAC"), evcodes = T)
    #profiler <- gprofiler(gene,organism = "mmusculus", src_filter = c("KEGG","REAC"))
    
    if(length(profiler) == 0) {pathway[[i]] <- 0
    pathway.genes[[i]] <- 0}
    else {pathway[[i]] <- profiler$result
    pathway.genes[[i]] <- as.list(profiler$result$intersection)
    pathway.genes[[i]] <- stringr::str_split(pathway.genes[[i]],",")
    names(pathway.genes[[i]]) <- profiler$result$term_name}
  }
  
  saveRDS(pathway,here("Data",folder,pathwayfile))
  saveRDS(pathway.genes,here("Data",folder,genefile))
}

## GOST PLOTS

gost.plot<-function(annotation, organism=NULL, src_filter=NULL, folder = NULL){
  gprofiler2::set_base_url("https://biit.cs.ut.ee/gprofiler_archive3/e104_eg51_p15/")
  Matched.module <- split(annotation$Gene, annotation$X..Module.)
  for(i in 1:length(Matched.module)){
    name <- str_split(names(Matched.module)[i],"_")
    name <- str_split(name[[1]][2],"\"")
    gene <- as.vector(Matched.module[[i]])
    
    cat("\n###",name[[1]],"\n")
    profiler <- gost(gene,organism = "mmusculus", correction_method = "fdr", sources = c("KEGG","REAC"), evcodes = T)
    
    if(length(profiler) == 0) {next}
    else {print(htmltools::tagList(gostplot(profiler)))}
    cat("\n \n")
  }
}

## Function that lists pathways
pathways.list <-function(pathway, matched){
  module.names <- unique(matched$X..Module.)
  name <- str_split(module.names,"_")
  samples <-c()
  for(i in 1:length(name)){
    samples[[i]] <- name[[i]][2]
  }  
  name <- str_split(samples,"\"")
  
  pathway.print <- list()
  for(i in 1:length(name)){
    cat("\n###",name[[i]], "Module","\n")
    if(class(pathway[[i]]) == "numeric"){next}
    else {print(knitr::kable(pathway[[i]][c(11,3:6)], format = "markdown", col.names = c("Pathway","FDR p-value","Term Size", "Query Size", "Overlap Size")))
    data <- as.data.frame(pathway[[i]][c(11,3:6)])
    names(data) <- c("Pathway","FDR p-value","Term Size", "Query Size", "Overlap Size") 
    output_name <- paste(name[[i]][1],"Module Pathways")
    print(data %>%
            download_this(
              output_name = output_name,
              output_extension = ".csv",
              button_label = "Download data as csv",
              button_type = "info",
              has_icon = TRUE,
              icon = "fa fa-save"
            ))}
    cat("\n \n")
    
  }
}

## Function that lists genes in each pathway
genes.pathways.list <-function(pathways, matched){
  module.names <- unique(matched$X..Module.)
  name <- str_split(module.names,"_")
  
  samples <-c()
  for(i in 1:length(name)){
    samples[[i]] <- name[[i]][2]
  }  
  name <- str_split(samples,"\"")
  
  for(i in 1:length(name)){
    cat("\n###",name[[i]], "Module","\n")
    
    print(knitr::kable(pathways[[i]]$intersection, col.names = "Genes in Pathway"))
    cat("\n \n")
  }
}

## Plot Module membership vs. Gene significance
MMvGS <- function(modulename, geneModuleMembership, geneTraitSignificance, MEs){
  module <- modulename
  modNames <- substring(names(MEs), 1)
  column <- match(module, modNames)
  moduleGenes = name==module
  external_gene <- Matched %>% 
    filter((Gene %in% rownames(geneModuleMembership[moduleGenes,]))) %>%
    dplyr::select(Gene, external_gene_name) %>% 
    column_to_rownames(var = "Gene") 
  external_gene <- full_join(rownames_to_column(external_gene), rownames_to_column(geneModuleMembership[moduleGenes,]), by = "rowname") 
  external_gene <- external_gene %>%
    column_to_rownames(var = "rowname")
  
  scattermembershiptrait <- as.data.frame(cbind(abs(geneModuleMembership[moduleGenes, column]),abs(geneTraitSignificance[moduleGenes, 1])))
  rownames(scattermembershiptrait) <- rownames(geneModuleMembership[moduleGenes,])
  scattermembershiptrait <- full_join(rownames_to_column(external_gene), rownames_to_column(scattermembershiptrait), by = "rowname")
  scattermembershiptrait <- scattermembershiptrait %>% 
    dplyr::select("rowname","external_gene_name", "V1", "V2")
  scattermembershiptrait <- scattermembershiptrait %>% filter(!is.na(V1)) %>% filter(!is.na(V2))
  
  module2 <- c()
  for(i in 1:length(scattermembershiptrait[,4])){
    if(scattermembershiptrait[,4][i] <= 0.05){module2[i] <- "lawngreen"}
    else {module2[i] <- module}
  }
  sizeGrWindow(7, 7)
  par(mfrow = c(1,1))
  p <- plotly::plot_ly(scattermembershiptrait, x = scattermembershiptrait[,3], y = scattermembershiptrait[,4]) %>% 
    plotly::add_trace(marker = list(size = 5, color = module2, line = list(color = "grey50", width = 1)), type = "scatter", mode = "markers", text = scattermembershiptrait[,2], hovertext = scattermembershiptrait[,1], hovertemplate = paste('Gene name: %{text}',
                                                                                                                                                                                                                                       '<br>Gene: %{hovertext}',
                                                                                                                                                                                                                                       '<br>MM: %{x}<br>',
                                                                                                                                                                                                                                       'P-value: %{y}')) %>%
    plotly::layout(title = paste("Module membership vs. gene significance\n"), xaxis = list(title = paste("Module Membership in", module, "module")), yaxis = list( title = paste("Gene significance for", names(geneTraitSignificance))), showlegend =FALSE)
  return(p)
}

## Plot barplots of eigengenes and RNA expression

eigen.expression <-function(eigens, expression){
  module.names <- colnames(eigens)[-1]
  name <- stringr::str_split(module.names,"_")
  
  samples <-c()
  for(i in 1:length(name)){
    samples[[i]] <- name[[i]][2]
  }  
  name <- stringr::str_split(samples,"\"")
  
  for(i in 1:length(name)){
    module.genes <- modules %>% 
      filter(Module == module.names[i])
    
    expression2 <- expression %>%
      dplyr::select(match(module.genes$Gene,colnames(expression)))
    
    expression$Treatment <- gsub("None", "Control", expression$Treatment)
    expression$Tissue<- gsub("Pre-frontal Cortex", "Prefrontal Cortex", expression$Tissue)
    expression$Tissue<- gsub("Hypothanamus","Hypothalamus", expression$Tissue)
    
    cat("\n###",name[[i]], "Module","\n")
    breaks<-seq(0,11,length.out = 50)
    expression2 <- t(expression2)
    anno <- as.data.frame(expression[,c("Time", "Treatment", "Tissue")])  
    rownames(anno) <- colnames(expression2)
    colors <- list(Tissue = c("Spleen" = "coral3", "Kidney" = "deeppink4", "Liver" = "palegreen3", "Prefrontal Cortex" = "royalblue4", "Heart" = "darkorange", "Hippocampus" = "darkgoldenrod3", "Hypothalamus" = "thistle3", "Skeletal Muscle" = "firebrick", "Small Intestine" = "sienna3"), Treatment = c("Control" = "aquamarine3", "2DG" = "darkolivegreen"), Time = c("4 wks" = "orchid4", "96 hrs" = "deepskyblue2"))
    g <- pheatmap(expression2, cluster_rows = T, cluster_cols = T, show_rownames = F, scale = "none", breaks = breaks, color = colorRampPalette(brewer.pal(12, "Paired"))(50), annotation_col = anno, annotation_colors = colors, fontsize_col = 8)
    p <- ggplot(eigens) + geom_bar(aes(eigens$X,unlist(eigens[i+1])), stat = "identity", color=name[[i]], fill = name[[i]]) + basic_theme + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) + labs(x = "Samples", y = "Module Summary Eigengene") + geom_hline(yintercept = 0, color = "black", size = .25) 
      
    print(g)
    print(p)
    cat("\n \n")
  }
}

##Dotplot
dot.plot <-function(eigens){
  module.names <- colnames(eigens)[-1]
  name <- str_split(module.names,"_")
  
  samples <-c()
  for(i in 1:length(name)){
    samples[[i]] <- name[[i]][2]
  }  
  name <- stringr::str_split(samples,"\"")
  
  for(i in 1:length(name)){
    module.genes <- modules %>% 
      filter(Module == module.names[i])
    
    eigenlength <- stringr::str_split(eigens$X,":")
    n <- length(eigenlength[[1]])
    
    eigens.name <- stringr::str_split_fixed(eigens$X,":", length(eigenlength[[1]]))
    if(dim(eigens.name)[2] == 2) {colnames(eigens.name) <- c("Sample","Treatment")
    } else if(dim(eigens.name)[2] == 3) {colnames(eigens.name) <- c("Sample","Time","Treatment")
    } else if(dim(eigens.name)[2] == 4) {colnames(eigens.name) <- c("Sample","Time","Treatment","Tissue")}
    
    if(n == 3) {assign(paste0(colnames(eigens.name)[2], "_by_", colnames(eigens.name)[3]), paste(eigens.name[,2],eigens.name[,3]))
      eigens2 <- cbind(eigens.name[,-1],Time_by_Treatment,eigens[-1])
    } else if(n == 4) {assign(paste0(colnames(eigens.name)[2], "_by_", colnames(eigens.name)[3]), paste(eigens.name[,2],eigens.name[,3]))
      assign(paste0(colnames(eigens.name)[2], "_by_", colnames(eigens.name)[4]), paste(eigens.name[,2],eigens.name[,4]))
      assign(paste0(colnames(eigens.name)[3], "_by_", colnames(eigens.name)[4]), paste(eigens.name[,3],eigens.name[,4]))
      assign(paste0(colnames(eigens.name)[2], "_by_", colnames(eigens.name)[3], "_by_", colnames(eigens.name)[4]), paste(eigens.name[,2],eigens.name[,3],eigens.name[,4]))
      eigens2 <- cbind(eigens.name[,-1],Time_by_Treatment,Time_by_Tissue,Treatment_by_Tissue,Time_by_Treatment_by_Tissue,eigens[-1])
    }     
    
    if(n == 2) {n = 1
    } else if(n == 3) {n = 3
    } else if(n == 4) {n = 7}
    cat("\n###",name[[i]], "Module","{.tabset .tabset-fade .tabset-pills}","\n")
    for(j in 1:n){
      cat("\n####",colnames(eigens2)[j],"\n")
      p <- ggplot(eigens2,aes(x = as.factor(eigens2[,j]),y = unlist(eigens2[i+n]))) + 
        geom_dotplot(binaxis='y', stackdir='center', stackratio=1.5, dotsize=.5, 
                     color=name[[i]], fill = name[[i]]) + basic_theme + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14), 
              axis.text.y = element_text(size = 14), axis.title.x = element_text(size = 16),
              axis.title.y = element_text(size = 16)) + 
        labs(x = "Samples", y = "Module Summary Eigengene") + 
        geom_hline(yintercept = 0, color = "black", size = .25) +
        theme(plot.background = element_rect(fill = "linen"),
              panel.background = element_rect(fill = "linen",
                                              colour = "linen",
                                              size = 0.5, linetype = "solid"),
              panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                              colour = "grey87")) +
                scale_x_discrete(labels=c("2DG" = "2DG", "None" = "Control"))
      
      print(p)
      cat("\n \n")
    }
    cat("\n \n")
  }
}

##Plot the number of times a pathway appears across all modules

count.pathways <- function(pathways){
  #filter.pathways <- list.filter(pathways, sum(query.number) >= 1)
  
  list.pathways <- list()
  for(i in 1:length(pathways)){
    if(class(pathways[[i]]) == "numeric") {next}
    else {list.pathways[[i]] <- pathways[[i]]$term_name}
  }
  
  list.pathways <- unlist(list.pathways)
  
  print(length(list.pathways))
  print(length(unique(list.pathways)))
  
  count <- as.data.frame(table(list.pathways)) %>% arrange(desc(Freq))
}
  


# Plot pathways shared genes between tissues in jaccard overlaps

genes.jaccard.pathways <- function(Matched, module){
  gprofiler2::set_base_url("https://biit.cs.ut.ee/gprofiler_archive3/e104_eg51_p15")
  genes.module <- Matched[which(str_detect(names(Matched),module))]
  genes.module <- map_df(genes.module, ~as.data.frame(.))
  gene <- as.vector(genes.module$Genes)
  profiler <- gost(gene,organism = "mmusculus", correction_method = "fdr", sources = c("KEGG","REAC"))
}
  



##Plot the number of times a pathway appears across all modules related to trait

count.pathways.traits <- function(pathways){
  #filter.pathways <- list.filter(pathways, sum(query.number) >= 1)
  
  list.pathways <- list()
  for(i in 1:length(pathways)){
    if(class(pathways[[i]]) == "numeric") {next}
    else {list.pathways[[i]] <- pathways[[i]]$term_name}
  }
  
  list.pathways <- unlist(list.pathways)
  
  print(length(list.pathways))
  print(length(unique(list.pathways)))
  
  count <- as.data.frame(table(list.pathways)) %>% arrange(desc(Freq))
  return(count)
}


## Hub genes
hub.genes <- function(net.deg, module.labels,modules,gene.name){
  hub <- list()
  for(m in module.labels) {
    name <- str_split(m,"_")[[1]][2]
    hub.met <- net.deg[modules$Gene,] %>%
      dplyr::mutate(Gene=modules$Gene) %>%
      dplyr::filter(modules$Module == m) %>%
      dplyr::arrange(desc(kWithin)) %>%
      dplyr::select(Gene, kWithin)
    gene <- gene.name$external_gene_name[match(hub.met$Gene,gene.name$ensembl_gene_id)]
    
    hub[[m]] <- cbind(gene,hub.met)
  }
  return(hub)
}

## Hub genes ANOVA and barplots
hub.anova <- function(net.deg, module.labels, modules,gene.name, data, factors = NULL){
  for(m in module.labels) {
    name <- str_split(m,"_")[[1]][2]
    hub.met <- net.deg[modules$Gene,] %>%
      dplyr::mutate(Gene=modules$Gene) %>%
      dplyr::filter(modules$Module == m) %>%
      dplyr::arrange(desc(kWithin)) %>%
      dplyr::select(Gene, kWithin)
    gene <- gene.name$external_gene_name[match(hub.met$Gene,gene.name$ensembl_gene_id)]
    cat("\n###",name, "\n")
    for(i in 1:10){
      cat("\n####",hub.met$Gene[i],"\n")
      n <- length(factors)
      m <- art(data = data, data[,hub.met$Gene[i]] ~ as.factor(Time)*as.factor(Treatment))
      model <- anova(m)
      panderOptions('knitr.auto.asis', FALSE)
      print(pander(model))
      cat("\n \n")
    }
    cat("\n \n")
  }
}

## plot interaction
plot.interaction <- function(data, var1, var2, var3=NULL, resp) {
  plot.data <- data %>%
    dplyr::group_by(!!sym(var1), !!sym(var2)) %>%
    dplyr::summarise(Mean=mean(!!sym(resp)), .groups="drop")
  
  p <- ggplot(plot.data) +
    geom_line(aes(x=!!sym(var1), y=Mean, color=!!sym(var2), group=!!sym(var2))) +
    geom_point(aes(x=!!sym(var1), y=Mean, color=!!sym(var2))) +
    ylab(resp) + ylim(min(data[,resp]), max(data[,resp])) +
    basic_theme
  
  return(p)
}

## eigenmetabolite
eigenmetabolite <- function(factors, data){
  for (trait in factors) {
    eigen.strat <- as.data.frame(do.call(cbind, lapply(
      split(module.eigens[rownames(data),], data[,trait]), 
      function(x) colSums(x) / nrow(x)
    )))
    eigen.strat$Module = sapply(str_split(rownames(eigen.strat),"_"),"[",2)
    
    plot.data <- tidyr::gather(eigen.strat, key=!!trait, "Value", -Module)
    
    abs.max <- max(abs(plot.data$Value))
    
    p <- ggplot(plot.data) +
      geom_tile(aes_string(x=trait, y="Module", fill="Value")) +
      scale_fill_distiller(palette="PRGn", limits=c(-abs.max, abs.max)) +
      basic_theme + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      theme(axis.text.y = element_text(size = 14))
    
    cat("\n###",trait,"\n")
    print(p)
    cat("\n \n")
  }
}

## Connectivity
label_facet <- function(original_var, custom_name){
  lev <- levels(as.factor(original_var))
  lab <- paste0(custom_name, ": ", lev)
  names(lab) <- lev
  return(lab)  
}

connectivity <- function(factors, data, modules){
  for (trait in factors) {
    gs <- as.data.frame(matrix(nrow=nrow(modules), ncol=2))
    gs[,1] = modules$Gene
    for (i in 1:nrow(modules)) {
      metab <- modules$Gene[i]
      gs[i,2] = -log10(summary(aov(data[,metab] ~ data[,trait]))[[1]][["Pr(>F)"]][1])
    }
    gs <- gs %>%
      dplyr::select(Gene=1, Log.10.P.Value=2) %>%
      dplyr::mutate(Module=modules$Module) %>%
      dplyr::mutate(Within.Degree=net.deg[modules$Gene,"kWithin"])
    
    plot.data <- gs %>%
      dplyr::group_by(Module) %>%
      dplyr::summarise(Mean=mean(Log.10.P.Value), .groups="drop") %>%
      dplyr::mutate(Color=unlist(lapply(str_split(Module, "_"), function(x) x[2])))
    
    p <- ggplot(plot.data) +
      geom_bar(aes(x=Module, y=Mean), fill=plot.data$Color, color="black", stat="identity") +
      ylab(paste0("Mean -Log10(P-Value) with ", trait)) +
      basic_theme + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
    
    plot.data <- gs %>%
      dplyr::mutate(Color=unlist(lapply(str_split(Module, "_"), function(x) x[2])))
    
    q <- ggplot(plot.data, aes(x=Within.Degree, y=Log.10.P.Value)) +
      geom_point(color="black", size=1.5) +
      geom_point(color=plot.data$Color, size=1) +
      geom_smooth(formula="y~x", method="lm") +
      #stat_cor(aes(label = ..p.label..), label.x = 0, label.y=-0.75) +
      xlab("Within-Module Connectivity") +
      ylab(paste0("-Log10(P-Value) with ", trait)) +
      facet_wrap(~Color, ncol=2) +
      basic_theme + theme(strip.background=element_rect(color="black", fill="white"), strip.text=element_text(face="bold"))
    
    cat("\n###", trait, "\n")
    print(p)
    
    #cat("\n####",Color,"\n")
    #print(q)
    #cat("\n \n")
    #print(suppressWarnings(plot_grid(p, q, rel_widths=c(1, 4))))
    cat("\n \n")
    #ggsave2(paste0("results/4a_module_analysis_combined/metabolite_connectivity/metabolite_connectivity_", tolower(trait), "_association.png"), width=12, height=8, dpi=600)
    #ggsave2(paste0("results/4a_module_analysis_combined/metabolite_connectivity/metabolite_connectivity_", tolower(trait), "_association.svg"), width=12, height=8, dpi=600)
  }
}

## gene significance
gene.significance <- function(factors, data, modules){
  for (trait in factors) {
    gs <- as.data.frame(matrix(nrow=nrow(modules), ncol=3))
    genes <- modules$Gene
    genes <- as.data.frame(genes, col.names = "Gene")
    external_gene <- Matched %>% 
      filter((genes %in% genes)) %>%
      dplyr::select(Gene, external_gene_name) 
    external_gene <- left_join(genes, external_gene, by = c("genes" = "Gene"))
    gs[,1] = external_gene$external_gene_name
    
    for (i in 1:nrow(modules)) {
      metab <- modules$Gene[i]
      gs[i,2] = -log10(summary(aov(data[,metab] ~ data[,trait]))[[1]][["Pr(>F)"]][1])
    }
    gs <- gs %>%
      dplyr::select(Gene=1, Log.10.P.Value=2) %>%
      dplyr::mutate(Module=modules$Module) %>%
      dplyr::mutate(Within.Degree=net.deg[modules$Gene,"kWithin"])
    
    plot.data <- gs %>%
      dplyr::mutate(Color=unlist(lapply(str_split(Module, "_"), function(x) x[2]))) %>%
      dplyr::mutate(Gene=Gene)
    plot.data$Color <- as.factor(plot.data$Color)
    
    bar.data <- gs %>%
      dplyr::group_by(Module) %>%
      dplyr::summarise(Mean=mean(Log.10.P.Value), .groups="drop") %>%
      dplyr::mutate(Color=unlist(lapply(str_split(Module, "_"), function(x) x[2])))
    bar.data$Color <- as.factor(bar.data$Color)
    
    cat("\n###",trait,"{.tabset .tabset-fade .tabset-pills}","\n")
    
    for(i in 1:nlevels(plot.data$Color)){
      data2 <- plot.data %>%
        filter(Color == levels(Color)[i])
      data3 <- bar.data %>%
        filter(Color == levels(Color)[i])
      
      p <- ggplot(data3) +
        geom_bar(aes(x=Module, y=Mean), fill=data3$Color, color="black", stat="identity") +
        ylab(paste0("Mean -Log10(P-Value) with ", trait)) +
        basic_theme + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
      
      q <- ggplot(data2, aes(x=Within.Degree, y=Log.10.P.Value)) + #, text = paste("Gene:", Gene))) +
        geom_point(color="black", size=1.5) +
        geom_point(color=data2$Color, size=1, aes(text = paste("Gene:", Gene))) +
        geom_smooth(formula="y~x", method="lm") +
        #stat_cor(aes(label = ..p.label..), label.x = 0, label.y=-0.75) +
        labs(x = "Within-Module Connectivity", y = paste0("-Log10(P-Value) with ", trait)) +
        #facet_wrap(~Color, ncol=4) +
        basic_theme #+ theme(strip.background=element_rect(color="black", fill="white"), strip.text=element_text(face="bold"))
      w <- ggplotly(q)
      
      #return(w)
      cat("\n####",levels(data2$Color)[i],"\n")
      print(htmltools::tagList(subplot(p, w, shareX = F, shareY = F, margin = .05, widths = c(0.3, 0.7)) %>%
                                 layout(yaxis = list(tickfont = list(size = 15), title = paste0("Mean -Log10(P-Value) with ", trait)), 
                                        xaxis = list(tickfont = list(size = 15), title = paste(levels(data2$Color)[i],"Module")), 
                                        xaxis2 = list(tickfont = list(size = 15), title = "Within-Module Connectivity"), 
                                        yaxis2 = list(tickfont = list(size = 15), title = paste0("-Log10(P-Value) with ", trait)))))
      # print(knitr::kable(data2[,c(1,2,4)]))
      # download <- data2[,c(1,2,4)]
      # names(download) <- c("Gene","Log 10 p-value","Connectivity")
      # output_name <- paste(levels(data2$Color)[i],"Module Hub Genes")
      # print(download %>%
      #   download_this(
      #     output_name = output_name,
      #     output_extension = ".csv",
      #     button_label = "Download data as csv",
      #     button_type = "info",
      #     has_icon = TRUE,
      #     icon = "fa fa-save"
      #   ))
      cat("\n \n")
      
    }
    # q <- ggplot(plot.data, aes(x=Within.Degree, y=Log.10.P.Value, text = paste("Gene:", Gene))) +
    #   geom_point(color="black", size=1.5) +
    #   geom_point(color=plot.data$Color, size=1) +
    #   geom_smooth(formula="y~x", method="lm") +
    #   #stat_cor(aes(label = ..p.label..), label.x = 0, label.y=-0.75) +
    #   xlab("Within-Module Connectivity") +
    #   ylab(paste0("-Log10(P-Value) with ", trait)) +
    #   facet_wrap(~Color, ncol=4) +
    #   basic_theme + theme(strip.background=element_rect(color="black", fill="white"), strip.text=element_text(face="bold"))
    # w <- ggplotly(q) 
    # 
    # return(w)
    cat("\n \n")
  }
}

## ANOVA for genes in pathways

Pathway.gene.ANOVA <- function(modules, WGCNA.gene) {
  
  h <- c()
  add <- c()
  for(m in 1:length(modules)) {
    cat("\n###",modules[m], "{.tabset .tabset-fade .tabset-pills}","\n")
    
    for(j in 1:length(WGCNA.gene[[m]])) {
      
      cat("\n####",names(WGCNA.gene[[m]][j]),"{.tabset .tabset-fade .tabset-pills}","\n")
      
      h[j] <- length(WGCNA.gene[[m]][[j]])
      add[j] <- sum(h) - h[j]
      
      
      for(k in 1:length(WGCNA.gene[[m]][[j]])) {
        if(WGCNA.gene[[m]][[j]] == 0) {"NA"}
        
        else{
          cat("\n#####",names(stat.table[k + add[j]]),"{.tabset .tabset-fade .tabset-pills}","\n")
          
          table.gene <- as.data.frame(stat.table[[k + add[j]]])
          
          print(knitr::kable(table.gene) %>%
                  column_spec(column = 8, color = "white", background = spec_color(table.gene[,8], end = 0.05, scale_from = c(0,0.05))))
          
          output_name <- paste("Genes in", names(stat.table[k + add[j]]), "ANOVA")
          print(table.df %>%
                  download_this(
                    output_name = output_name,
                    output_extension = ".csv",
                    button_label = "Download data as csv",
                    button_type = "info",
                    has_icon = TRUE,
                    icon = "fa fa-save"
                    ))
          cat("\n \n")
        }
      }
      cat("\n \n")
    }
    cat("\n \n")
  }
}


## Create gene set for GSVA

#ensembl_location<-readRDS(here("Data","Ensembl_gene_id_and_location.RData"))

#profiler <- gprofiler2::gost(ensembl_location$ensembl_gene_id,organism = "mmusculus", correction_method = "fdr", sources = c("KEGG","REAC"), evcodes = T, significant = F)

#saveRDS(profiler, here("Data","gene_set_from_ensembl_gprofiler_KEGG_REACTOME.RData"))

#ensembl.pathway <- readRDS(here("Data","gene_set_from_ensembl_gprofiler_KEGG_REACTOME.RData"))

#pathway.genes <- as.list(ensembl.pathway$result$intersection)
#pathway.genes <- stringr::str_split(pathway.genes, ",")
#names(pathway.genes) <- gene.pathway$result$term_name

#saveRDS(pathway.genes, here("Data","gene_set_from_ensembl_gprofiler_KEGG_REACTOME_gene_annotation.RData"))

## Run GSVA

GSVA.modules <- function(modules, logdata, data){ ## data contains sample info
  gene.pathway <- readRDS(here("Data","gene_set_from_ensembl_gprofiler_KEGG_REACTOME_gene_annotation.RData"))

  Matched.module <- split(modules$Gene, modules$Module)
  module.gene <- list()
  for(i in 1:length(Matched.module)){
    module <- as.vector(Matched.module[[i]])
    module.gene[[i]] <- logdata[match(module,rownames(logdata)),] #%>% 
      #select(one_of(module))
    }
  names(module.gene) <- names(Matched.module)

  tES <- list()
  for(s in 1:length(module.gene)){
    ES <- gsva(as.matrix(module.gene[[s]]), gset.idx.list = gene.pathway, method = "plage", verbose = F)
    tES1 <- t(ES)
    tES[[s]] <- tES1[,colMeans(tES1) < 1]
    }
  names(tES) <- names(module.gene)
  
  
  table.module <-list()
  for(j in 1:length(tES)){
    name.tES <- colnames(tES[[j]])
    tES[[j]] <- as.data.frame(tES[[j]])
    colnames(tES[[j]]) <- name.tES
    
  name <- str_split(names(tES)[j],"_")[[1]][2]
  
  table.df <- data.frame(matrix(NA, 0, 5))
  
  for(f in 1:ncol(tES[[j]])){
    m <- art(data = data, tES[[j]][,f] ~ Time*Treatment)
    model <- anova(m)
    adjust <- p.adjust(model$`Pr(>F)`, method = "BH")
    
    pathway.name <- str_split(colnames(as.data.frame(tES[j]))[f],"_")[[1]][2]
    
    avg <- as.data.frame(cbind("Treatment" = as.factor(data$Treatment), tES[[j]][,f])) %>%
      group_by(Treatment) %>%
      summarise(median(V2))
    
    DG <- as.numeric(avg[1,2])
    control <- as.numeric(avg[2,2])
    
    expression <- if(DG < control) { "down"
    } else {"up"}
    
    table.df[f,] <- cbind(pathway.name, adjust[1], adjust[2], adjust[3], expression)
    
  }
  
  
  colnames(table.df) <- c("Pathway","Time", "Treatment","Time by Treatment","2DG expression compared to Control")
  
  table.df$Time <- as.numeric(table.df$Time)
  table.df$Treatment <- as.numeric(table.df$Treatment)
  table.df$`Time by Treatment` <- as.numeric(table.df$`Time by Treatment`)
  table.module[[j]] <- table.df
  }
  
  names(table.module) <- names(module.gene)
  
  return(list(table.module = table.module, tES = tES))
}

print.gsva.table <- function(w = 1, table.module, tES){
  name <- str_split(names(tES)[w],"_")[[1]][2]
  
  table.module2 <- table.module[[w]]
  
  if(length(unique(table.module2$`2DG expression compared to Control`)) == 1 && unique(table.module2$`2DG expression compared to Control`) == "up"){ 
    DT::datatable(table.module2, extensions = 'Buttons',
                                           rownames = FALSE, 
                                           filter="top",
                                           options = list(dom = 'Blfrtip',
                                                          buttons = c('copy', 'csv', 'excel'),
                                                          lengthMenu = list(c(10,25,50,-1),
                                                                            c(10,25,50,"All")), 
                                                          scrollX= TRUE), class = "display") %>%
                               formatStyle(
                                 'Time',
                                 color = styleInterval(c(0,0.05), c('white',"white","black")),
                                 backgroundColor = styleInterval(0.05, c('#440154FF',"white"))) %>%

                               formatStyle(
                                 'Treatment',
                                 color = styleInterval(c(0,0.05), c('white',"white","black")),
                                 backgroundColor = styleInterval(0.05, c("#440154FF","white"))) %>%

                               formatStyle(
                                 'Time by Treatment',
                                 color = styleInterval(c(0,0.05), c('white',"white","black")),
                                 backgroundColor = styleInterval(0.05, c('#440154FF',"white"))) %>%

                               formatStyle(
                                 '2DG expression compared to Control',
                                 color = "white",
                                 backgroundColor = styleEqual(unique(table.module2$`2DG expression compared to Control`), c("#CDB1AD")))

  } else if(length(unique(table.module2$`2DG expression compared to Control`)) == 1 && unique(table.module2$`2DG expression compared to Control`) == 
            "down"){ 
    DT::datatable(table.module2, extensions = 'Buttons',
                                           rownames = FALSE, 
                                           filter="top",
                                           options = list(dom = 'Blfrtip',
                                                          buttons = c('copy', 'csv', 'excel'),
                                                          lengthMenu = list(c(10,25,50,-1),
                                                                            c(10,25,50,"All")), 
                                                          scrollX= TRUE), class = "display") %>%
                               formatStyle(
                                 'Time',
                                 color = styleInterval(c(0,0.05), c('white',"white","black")),
                                 backgroundColor = styleInterval(0.05, c('#440154FF',"white"))) %>%

                               formatStyle(
                                 'Treatment',
                                 color = styleInterval(c(0,0.05), c('white',"white","black")),
                                 backgroundColor = styleInterval(0.05, c("#440154FF","white"))) %>%

                               formatStyle(
                                 'Time by Treatment',
                                 color = styleInterval(c(0,0.05), c('white',"white","black")),
                                 backgroundColor = styleInterval(0.05, c('#440154FF',"white"))) %>%

                               formatStyle(
                                 '2DG expression compared to Control',
                                 color = "white",
                                 backgroundColor = styleEqual(unique(table.module2$`2DG expression compared to Control`), c("#2A385B")))

      } else if(length(unique(table.module2$`2DG expression compared to Control`)) == 2){ 
    DT::datatable(table.module2, extensions = 'Buttons',
                  rownames = FALSE, 
                  filter="top",
                  options = list(dom = 'Blfrtip',
                                 buttons = c('copy', 'csv', 'excel'),
                                 lengthMenu = list(c(10,25,50,-1),
                                                   c(10,25,50,"All")), 
                                 scrollX= TRUE), class = "display") %>%
      formatStyle(
        'Time',
        color = styleInterval(c(0,0.05), c('white',"white","black")),
        backgroundColor = styleInterval(0.05, c('#440154FF',"white"))) %>%

      formatStyle(
        'Treatment',
        color = styleInterval(c(0,0.05), c('white',"white","black")),
        backgroundColor = styleInterval(0.05, c("#440154FF","white"))) %>%

      formatStyle(
        'Time by Treatment',
        color = styleInterval(c(0,0.05), c('white',"white","black")),
        backgroundColor = styleInterval(0.05, c('#440154FF',"white"))) %>%

      formatStyle(
        '2DG expression compared to Control',
        color = "white",
        backgroundColor = styleEqual(unique(table.module2$`2DG expression compared to Control`), c("#2A385B","#CDB1AD")))

  } else if(unique(table.module2$`2DG expression compared to Control`) == "up") {
    DT::datatable(table.module2, extensions = 'Buttons',
                  rownames = FALSE, 
                  filter="top",
                  options = list(dom = 'Blfrtip',
                                 #buttons = c('copy', 'csv', 'excel'),
                                 lengthMenu = list(c(10,25,50,-1),
                                                   c(10,25,50,"All")), 
                                 scrollX= TRUE),class = "display") %>%
      formatStyle(
        'Time',
        color = styleInterval(c(0,0.05), c('white',"white","black")),
        backgroundColor = styleInterval(0.05, c('#440154FF',"white"))) %>%

      formatStyle(
        'Treatment',
        color = styleInterval(c(0,0.05), c('white',"white","black")),
        backgroundColor = styleInterval(0.05, c("#440154FF","white"))) %>%

      formatStyle(
        'Time by Treatment',
        color = styleInterval(c(0,0.05), c('white',"white","black")),
        backgroundColor = styleInterval(0.05, c('#440154FF',"white"))) %>%

      formatStyle(
        '2DG expression compared to Control',
        color = "white",
        backgroundColor = styleEqual(unique(table.module2$`2DG expression compared to Control`), "#CDB1AD"))

  } else if(unique(table.module2$`2DG expression compared to Control`) == "down") {
    DT::datatable(table.module2, extensions = 'Buttons',
                  rownames = FALSE, 
                  filter="top",
                  options = list(dom = 'Blfrtip',
                                 #buttons = c('copy', 'csv', 'excel'),
                                 lengthMenu = list(c(10,25,50,-1),
                                                   c(10,25,50,"All")), 
                                 scrollX= TRUE),class = "display") %>%
      formatStyle(
        'Time',
        color = styleInterval(c(0,0.05), c('white',"white","black")),
        backgroundColor = styleInterval(0.05, c('#440154FF',"white"))) %>%

      formatStyle(
        'Treatment',
        color = styleInterval(c(0,0.05), c('white',"white","black")),
        backgroundColor = styleInterval(0.05, c("#440154FF","white"))) %>%

      formatStyle(
        'Time by Treatment',
        color = styleInterval(c(0,0.05), c('white',"white","black")),
        backgroundColor = styleInterval(0.05, c('#440154FF',"white"))) %>%

      formatStyle(
        '2DG expression compared to Control',
        color = "white",
        backgroundColor = styleEqual(unique(table.module2$`2DG expression compared to Control`), "#2A385B"))

  }
  # output_name <- paste(name, "Module GSVA ANOVA")
  # print(table.module2 %>%
  #         download_this(
  #           output_name = output_name,
  #           output_extension = ".csv",
  #           button_label = "Download data as csv",
  #           button_type = "info",
  #           has_icon = TRUE,
  #           icon = "fa fa-save"
  #         ))
#}
}

## ANOVA for genes within WGCNA pathways time and treatment

ANOVA.gene.pathways <- function(WGCNA.gene, logdata){
  ensembl.location <- readRDS(here("Data","Ensembl_gene_id_and_location.RData"))
  
  logdata$Time <- as.factor(logdata$Time)
  logdata$Treatment <- as.factor(logdata$Treatment)
  logdata$Tissue <- as.factor(logdata$Tissue)
  
  for (i in 1:length(WGCNA.gene)){
  cat("\n###",modules[i],"{.tabset .tabset-fade .tabset-pills}"," \n")
    if(is.list(WGCNA.gene[[i]])==FALSE){print("No significantly overrepresented pathways reported")}
    else{
  for(j in 1:length(WGCNA.gene[[i]])){

    cat("\n####",names(WGCNA.gene[[i]])[j]," \n")
    frame <- matrix(data = NA, nrow = length(WGCNA.gene[[i]][[j]]), ncol = 8)
    
    for(k in 1:length(WGCNA.gene[[i]][[j]])){
      gene <- WGCNA.gene[[i]][[j]][k]
      column <- which(names(logdata)==gene)
      if(is.integer(which(names(logdata)==gene))==FALSE){
        frame[k,] <- c(gene, "NA","NA","NA","NA","NA","NA","NA")}
      else if(length(logdata[gene][,1])==0){next}
        else {
        m <- art(data = logdata, unlist(logdata[gene][,1]) ~ Time*Treatment)
        
        model <- anova(m)
        adjust <- p.adjust(model$`Pr(>F)`, method = "BH")
        num <- which(names(logdata)==gene)
        gene2 <- ensembl.location$external_gene_name[match(names(logdata)[num],ensembl.location$ensembl_gene_id)]
        
        frame[k,] <- c(gene2, names(logdata)[num], model[,7][2], adjust[2], model[,7][1], adjust[1], model[,7][3],adjust[3])
      }
    }
    
    frame <- as.data.frame(frame)
    
    print(knitr::kable(frame, col.names = c("External Gene Name", "Gene ID", "Treatment", "BH Treatment",
                                            "Time", "BH Time", "Treatment by time", "BH Treatment by time")) %>% 
            kableExtra::kable_styling(bootstrap_options = c("striped", "hover", "responsive")))
    
    output_name <- paste("Genes in", names(WGCNA.gene[[i]])[j], " pathway for Heart", modules[i])
    
    print(frame %>%
            download_this(
              output_name = output_name,
              output_extension = ".csv",
              button_label = "Download data as csv",
              button_type = "info",
              has_icon = TRUE,
              icon = "fa fa-save"
            ))
    
    cat("\n \n")
  }             
  cat("\n \n")
}
  }
}

## ANOVA for genes within WGCNA pathways tissue, time, and treatment

ANOVA.gene.tissues.pathways <- function(WGCNA.gene, logdata){
  ensembl.location <- readRDS(here("Data","Ensembl_gene_id_and_location.RData"))
  
  logdata$Time <- as.factor(logdata$Time)
  logdata$Treatment <- as.factor(logdata$Treatment)
  logdata$Tissue <- as.factor(logdata$Tissue)
  
  for (i in 1:length(WGCNA.gene)){
    cat("\n###",modules[i],"{.tabset .tabset-fade .tabset-pills}"," \n")
    if(is.list(WGCNA.gene[[i]])==FALSE){next}
    else{
      for(j in 1:length(WGCNA.gene[[i]])){
        
        cat("\n####",names(WGCNA.gene[[i]])[j]," \n")
        frame <- matrix(data = NA, nrow = length(WGCNA.gene[[i]][[j]]), ncol = 16)
        
        for(k in 1:length(WGCNA.gene[[i]][[j]])){
          gene <- WGCNA.gene[[i]][[j]][k]
          column <- which(names(logdata)==gene)
          if(is.integer(which(names(logdata)==gene))==FALSE){
            frame[k,] <- c(gene, "NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA")}
          else if(length(logdata[gene][,1])==0){next}
          else {
            m <- art(data = logdata, unlist(logdata[gene][,1]) ~ Time*Treatment*Tissue)
            
            model <- anova(m)
            adjust <- p.adjust(model$`Pr(>F)`, method = "BH")
            num <- which(names(logdata)==gene)
            gene2 <- ensembl.location$external_gene_name[match(names(logdata)[num],ensembl.location$ensembl_gene_id)]
            
            frame[k,] <- c(gene2, names(logdata)[num], model[,7][2], adjust[2], model[,7][1], adjust[1], model[,7][3],adjust[3], model[,7][4], adjust[4],model[,7][5], adjust[5],model[,7][6], adjust[6],model[,7][7], adjust[7])
          }
        }
        
        frame <- as.data.frame(frame)
        
        print(knitr::kable(frame, col.names = c("External Gene Name", "Gene ID", "Treatment", "BH Treatment","Time", "BH Time", "Tissue", "BH Tissue", "Treatment by time", "BH Treatment by time", "Time by Tissue","BH Time by Tissue", "Treatment by Tissue", "BH Treatment by Tissue", "Time by Treatment by Tissue", "BH Time by Treatment by Tissue")) %>%
          kable_styling(bootstrap_options = c("striped", "hover", "responsive")))
        
        output_name <- paste("Genes in", names(WGCNA.gene[[i]])[j], " pathway for Heart", modules[i])
        
        print(frame %>%
                download_this(
                  output_name = output_name,
                  output_extension = ".csv",
                  button_label = "Download data as csv",
                  button_type = "info",
                  has_icon = TRUE,
                  icon = "fa fa-save"
                ))
        cat("\n \n")
      }             
      cat("\n \n")
    }
  }
}

## Create wrapper for DT table with basic theme
DT.table <- function(data){DT::datatable(data, extensions = 'Buttons',
              rownames = FALSE, 
              filter="top",
              options = list(dom = 'Blfrtip',
                             buttons = c('copy', 'csv', 'excel'),
                             lengthMenu = list(c(10,25,50,-1),
                                               c(10,25,50,"All")), 
                             scrollX= TRUE), class = "display") 
}

## Create wrapper for DT table with column names for hub genes
DT.table.hub <- function(data){DT::datatable(data, extensions = 'Buttons',
                                         rownames = FALSE, 
                                         filter="top",
                                         options = list(dom = 'Blfrtip',
                                                        buttons = c('copy', 'csv', 'excel'),
                                                        lengthMenu = list(c(10,25,50,-1),
                                                                          c(10,25,50,"All")), 
                                                        scrollX= TRUE), class = "display", colnames = c("Gene name", "Ensembl Gene", "kWithin")) 
}

## Create wrapper for DT table with column names for Jaccard genes
DT.table.jaccard <- function(data){DT::datatable(data, extensions = 'Buttons',
                                             rownames = FALSE, 
                                             filter="top",
                                             options = list(dom = 'Blfrtip',
                                                            buttons = c('copy', 'csv', 'excel'),
                                                            lengthMenu = list(c(10,25,50,-1),
                                                                              c(10,25,50,"All")), 
                                                            scrollX= TRUE), class = "display", colnames = c("Pathway", "Term Size", "Query Size", "Intersection Size", "FDR"))
}

## Create wrapper for DT table with column names for Jaccard pathways
DT.table.jaccard.pathways <- function(data){DT::datatable(data, extensions = 'Buttons',
                                                 rownames = FALSE, 
                                                 filter="top",
                                                 options = list(dom = 'Blfrtip',
                                                                buttons = c('copy', 'csv', 'excel'),
                                                                lengthMenu = list(c(10,25,50,-1),
                                                                                  c(10,25,50,"All")), 
                                                                scrollX= TRUE), class = "display", colnames = c("Pathway", "Heart", "Hippo","Hypo","Kidney",
                                                                                                                "Liver","Cortex","Muscle","Intestine","Spleen", "Freq"))
}

## Create wrapper for DT table for frequency plot
DT.table.freq <- function(data){DT::datatable(data, extensions = 'Buttons',
                                         rownames = FALSE, 
                                         filter="top",
                                         options = list(dom = 'Blfrtip',
                                                        buttons = c('copy', 'csv', 'excel'),
                                                        lengthMenu = list(c(10,25,50,-1),
                                                                          c(10,25,50,"All")), 
                                                        scrollX= TRUE), class = "display", colnames = c("Pathway", "Frequency"))
}

## Create wrapper for DT table for pathways
DT.table.path <- function(data){DT::datatable(data, extensions = 'Buttons',
                                              rownames = FALSE, 
                                              filter="top",
                                              options = list(dom = 'Blfrtip',
                                                             buttons = c('copy', 'csv', 'excel'),
                                                             lengthMenu = list(c(10,25,50,-1),
                                                                               c(10,25,50,"All")), 
                                                             scrollX= TRUE), class = "display", colnames = c("Pathway","FDR p-value","Term Size", "Query Size", "Overlap Size"))
}

