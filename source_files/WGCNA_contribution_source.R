## time contribution
time.contribution <- function(phenotype, log.tdata.FPKM.subset, MEs){
# Define variable weight containing the weight column of datTrait
Time <- as.data.frame(phenotype$Time);
names(Time) <- "Time"

# names (colors) of the modules
modNames <- substring(names(MEs), 1)
geneModuleMembership <- as.data.frame(cor(log.tdata.FPKM.subset, MEs, use = "p"));
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) <- paste("MM", modNames, sep="");
names(MMPvalue) <- paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(log.tdata.FPKM.subset, Time, use = "p"));
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) <- paste("GS.", names(Time), sep="");
names(GSPvalue) <- paste("p.GS.", names(Time), sep="");                 

modulename <- modNames

for(i in 1:length(modulename)){
  module <- modulename[i]
  cat("###", modulename[i],"\n")
  print(htmltools::tagList(MMvGS(module, geneModuleMembership, geneTraitSignificance, MEs)))
  cat("\n \n")
}
}

## tissue contribution
tissue.contribution <- function(phenotype, log.tdata.FPKM.subset, MEs){
# Define variable weight containing the weight column of datTrait
Tissue <- as.data.frame(phenotype$Tissue);
names(Tissue) <- "Tissue"

# names (colors) of the modules
modNames <- substring(names(MEs), 1)
geneModuleMembership <- as.data.frame(cor(log.tdata.FPKM.subset, MEs, use = "p"));
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) <- paste("MM", modNames, sep="");
names(MMPvalue) <- paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(log.tdata.FPKM.subset, Tissue, use = "p"));
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) <- paste("GS.", names(Tissue), sep="");
names(GSPvalue) <- paste("p.GS.", names(Tissue), sep="");                 

modulename <- modNames

for(i in 1:length(modulename)){
  module <- modulename[i]
  cat("###", modulename[i],"\n")
  print(htmltools::tagList(MMvGS(module, geneModuleMembership, geneTraitSignificance, MEs)))
  cat("\n \n")
}
}

## treatment contribution
treat.contribution <- function(phenotype, log.tdata.FPKM.subset, MEs){
  # Define variable weight containing the weight column of datTrait
  Treatment <- as.data.frame(phenotype$Treatment);
  names(Treatment) <- "Treatment"
  
  # names (colors) of the modules
  modNames <- substring(names(MEs), 1)
  geneModuleMembership <- as.data.frame(cor(log.tdata.FPKM.subset, MEs, use = "p"));
  MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
  
  names(geneModuleMembership) <- paste("MM", modNames, sep="");
  names(MMPvalue) <- paste("p.MM", modNames, sep="");
  geneTraitSignificance = as.data.frame(cor(log.tdata.FPKM.subset, Treatment, use = "p"));
  GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
  names(geneTraitSignificance) <- paste("GS.", names(Treatment), sep="");
  names(GSPvalue) <- paste("p.GS.", names(Treatment), sep="");                 

  modulename <- modNames
  
  for(i in 1:length(modulename)){
    module <- modulename[i]
    cat("###", modulename[i],"\n")
    print(htmltools::tagList(MMvGS(module, geneModuleMembership, geneTraitSignificance, MEs)))
    cat("\n \n")
  }

}

## tissue time contribution
tissue.time.contribution <- function(phenotype, log.tdata.FPKM.subset, MEs){
  # Define variable weight containing the weight column of datTrait
  Tissue.Time <- as.data.frame(phenotype$Tissue.Time);
  names(Tissue.Time) <- "Tissue.Time"
  
  # names (colors) of the modules
  modNames <- substring(names(MEs), 1)
  geneModuleMembership <- as.data.frame(cor(log.tdata.FPKM.subset, MEs, use = "p"));
  MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
  
  names(geneModuleMembership) <- paste("MM", modNames, sep="");
  names(MMPvalue) <- paste("p.MM", modNames, sep="");
  geneTraitSignificance = as.data.frame(cor(log.tdata.FPKM.subset, Tissue.Time, use = "p"));
  GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
  names(geneTraitSignificance) <- paste("GS.", names(Tissue.Time), sep="");
  names(GSPvalue) <- paste("p.GS.", names(Tissue.Time), sep="");                 

  modulename <- modNames
  
  for(i in 1:length(modulename)){
    module <- modulename[i]
    cat("###", modulename[i],"\n")
    print(htmltools::tagList(MMvGS(module, geneModuleMembership, geneTraitSignificance, MEs)))
    cat("\n \n")
  }

}

## Treatment Time contribution
treat.time.contribution <- function(phenotype, log.tdata.FPKM.subset, MEs){
  # Define variable weight containing the weight column of datTrait
  Treat.Time <- as.data.frame(phenotype$Treat.Time);
  names(Treat.Time) <- "Treat.Time"
  
  # names (colors) of the modules
  modNames <- substring(names(MEs), 1)
  geneModuleMembership <- as.data.frame(cor(log.tdata.FPKM.subset, MEs, use = "p"));
  MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
  
  names(geneModuleMembership) <- paste("MM", modNames, sep="");
  names(MMPvalue) <- paste("p.MM", modNames, sep="");
  geneTraitSignificance = as.data.frame(cor(log.tdata.FPKM.subset, Treat.Time, use = "p"));
  GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
  names(geneTraitSignificance) <- paste("GS.", names(Treat.Time), sep="");
  names(GSPvalue) <- paste("p.GS.", names(Treat.Time), sep="");                 

  modulename <- modNames
  
  for(i in 1:length(modulename)){
    module <- modulename[i]
    cat("###", modulename[i],"\n")
    print(htmltools::tagList(MMvGS(module, geneModuleMembership, geneTraitSignificance, MEs)))
    cat("\n \n")
  }

}

## Treatment tissue contribution

treat.tissue.contribution <- function(phenotype, log.tdata.FPKM.subset, MEs){
# Define variable weight containing the weight column of datTrait
Treat.Tissue <- as.data.frame(phenotype$Treat.Tissue);
names(Treat.Tissue) <- "Treat.Tissue"

# names (colors) of the modules
modNames <- substring(names(MEs), 1)
geneModuleMembership <- as.data.frame(cor(log.tdata.FPKM.subset, MEs, use = "p"));
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) <- paste("MM", modNames, sep="");
names(MMPvalue) <- paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(log.tdata.FPKM.subset, Treat.Tissue, use = "p"));
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) <- paste("GS.", names(Treat.Tissue), sep="");
names(GSPvalue) <- paste("p.GS.", names(Treat.Tissue), sep="");                 

modulename <- modNames

for(i in 1:length(modulename)){
  module <- modulename[i]
  cat("###", modulename[i],"\n")
  print(htmltools::tagList(MMvGS(module, geneModuleMembership, geneTraitSignificance, MEs)))
  cat("\n \n")
}
}

## Time Treatment Tissue contribution

time.treat.tissue.contribution <- function(phenotype, log.tdata.FPKM.subset, MEs){
  
  # Define variable weight containing the weight column of datTrait
  Treat.Tissue.Time <- as.data.frame(phenotype$Treat.Tissue.Time);
  names(Treat.Tissue.Time) <- "Treat.Tissue.Time"
  
  # names (colors) of the modules
  modNames <- substring(names(MEs), 1)
  geneModuleMembership <- as.data.frame(cor(log.tdata.FPKM.subset, MEs, use = "p"));
  MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
  
  names(geneModuleMembership) <- paste("MM", modNames, sep="");
  names(MMPvalue) <- paste("p.MM", modNames, sep="");
  geneTraitSignificance = as.data.frame(cor(log.tdata.FPKM.subset, Treat.Tissue.Time, use = "p"));
  GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
  names(geneTraitSignificance) <- paste("GS.", names(Treat.Tissue.Time), sep="");
  names(GSPvalue) <- paste("p.GS.", names(Treat.Tissue.Time), sep="");                 

  modulename <- modNames
  
  for(i in 1:length(modulename)){
    module <- modulename[i]
    cat("###", modulename[i],"\n")
    print(htmltools::tagList(MMvGS(module, geneModuleMembership, geneTraitSignificance, MEs)))
    cat("\n \n")
  }

}