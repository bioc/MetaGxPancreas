## validate with `AnnotationHub::readMetadataFromCsv("TENxBrainData")`
## (above pkg directory)

infoString = c("BADEA", "BALAGURANATHAN" ,"CHEN", "COLLISSON", "GRUTZMANN", "ICGCSEQ", "ICGCMICRO",
               "KIRBY", "OUH", "PCSI", "PEI", "TCGA",
               "UNC", "WINTER", "ZHANG", "HAMIDI","YANG", "LUNARDI", "JANKY", "BAUER", "HAIDER")


main.data <- data.frame(
  Title = c("BADEA", "BALAGURANATHAN" ,"CHEN", "COLLISSON", "GRUTZMANN", "ICGCSEQ", "ICGCMICRO",
            "KIRBY", "OUH", "PCSI", "PEI", "TCGA",
            "UNC", "WINTER", "ZHANG", "HAMIDI","YANG", "LUNARDI", "JANKY", "BAUER", "HAIDER"
  ),
  Description = paste(
    "Pancreas cancer gene expression data ",
    c(
      "used in, Badea et al, Hepatogastroenterology 2008",
      "used in Balagurunathan et al, Mol Cancer Ther 2008",
      "used in, Chen et al, PLoS One 2015",
      "used in, Collisson et al.,Nat Med 2011",
      "used in, Grutzmann et al, Neoplasia, 2004",
      "used in, Nones et al, Int. J. Cancer, 2014",
      "used in, Bailey et al, Nature, 2016",
      "used in, Kirby et al., Mol Oncol 2016",
      "used in, Sandhu et al, Mol Onc, 2015",
      "used in, Notta et al, Nature 2016 ",
      "used in, Pei et al, Cancer Cell 2009",
      "used in, TCGA Research Network, Cancer Cell 2017",
      "used in, Moffitt et al, Nat Genet 2015",
      "used in, Winter et al, PLoS Comput Biol, 2012",
      "used in, Zhang et al, PLoS One 2012",
      "from Hamidi et al, not used in a paper",
      "used in, Yang et al, 2016, Cancer Research",
      "used in, Lunardi S et al, 2014, Oncotarget",
      "used in, Janky et al, BMC Cancer 2016",
      "used in, Bauer et al, 2016, Gastroenterology",
      "used in, haider et al, Genome medicine, 2014"
    )
  ),
  RDataPath = c(
    "MetaGxPancreas/BADEA.rds", "MetaGxPancreas/BALAGURANATHAN.rds", "MetaGxPancreas/CHEN.rds",
    "MetaGxPancreas/COLLISSON.rds", "MetaGxPancreas/GRUTZMANN.rds", "MetaGxPancreas/ICGCSEQ.rds",
    "MetaGxPancreas/ICGCMICRO.rds", "MetaGxPancreas/KIRBY.rds", "MetaGxPancreas/OUH.rds",
    "MetaGxPancreas/PCSI.rds", "MetaGxPancreas/PEI.rds", "MetaGxPancreas/TCGA.rds",
    "MetaGxPancreas/UNC.rds", "MetaGxPancreas/WINTER.rds", "MetaGxPancreas/ZHANG.rds",
    "MetaGxPancreas/HAMIDI.rds","MetaGxPancreas/YANG.rds", "MetaGxPancreas/LUNARDI.rds", 
    "MetaGxPancreas/JANKY.rds", "MetaGxPancreas/BAUER.rds", "MetaGxPancreas/HAIDER.rds"
  ),
  BiocVersion="3.7",
  Genome=c("Affymetrix Human Genome U133 Plus 2.0 Array", "Human 1A Microarray G4110A-G4110B", "Affymetrix,Rosetta-Merck RSTA Custom 2.0", "Affymetrix, array U133 Plus 2.0","Affymetrix GeneChip Human Genome HG-U133B",
           "Illumina HumanHT-12 V4.0 expression beadchip","Illumina HiSeq 2000-2500","Illumina, RNA seqencing HiSeq","Agilent-028004 SurePrint G3 Human GE 8x60K Microarray","Illumina HiSeq 2000-2500","Affymetrix Human Genome U133 Plus 2.0 Array","IlluminaHiSeq_RNASeqV2","Agilent-014850 Whole Human Genome Microarray 4x44K G4112F Agilent-014850 Whole Human Genome Microarray 4x44K G4112F ",
           "A-AFFY-44-Affymetrix GeneChip Human Genome U133 Plus 2.0 [HG-U133_Plus_2]","Affymetrix GeneChip Human Gene 1.0 ST arrays", "Agilent-012097 Human 1A Microarray (V2) G4110B", "Affymetrix GeneChip Human Gene 1.0 ST", " Agilent-014850 Whole Human Genome Microarray 4x44K G4112F",
           "Affymetrix Human Genome U219 Array", "Illumina human WG6 Expression BeadChip", "Affymetrix Human Exon 1.0 ST Array"
  ),
  SourceType="RData",
  SourceUrl=c("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse15471", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE11838",
              "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE57495", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE17891",
              "https://www.ebi.ac.uk/arrayexpress/experiments/E-MEXP-950/?query=pilarsky&s_page=1&s_pagesize=50", "http://icgc.org/icgc/cgp/68/304/798",
              "http://icgc.org/icgc/cgp/68/304/798", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE79670",
              "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60980", "Private",
              "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse16515", "https://portal.gdc.cancer.gov/projects/TCGA-PAAD",
              "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71729", "http://www.ebi.ac.uk/arrayexpress/experiments/E-MEXP-2780/",
              "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE28735", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE77858",
              "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62452", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE55643",
              "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62165", "https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-1791/",
              "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE56560"),
  SourceVersion="July 09 2020",
  Species="Homo sapien",
  TaxonomyId=9606,
  Coordinate_1_based=FALSE,
  DataProvider=c(rep("GEO", 4), "Array Express", rep("ICGC", 2), rep("GEO", 2), "Private", "GEO", "GDC", "GEO", "Array Express", "GEO", rep("GEO", 4), "Array Express", "GEO"),
  Maintainer="Michael Zon <michaelzon7@gmail.com>",
  RDataClass="SummarizedExperiment",
  DispatchClass="RDS",
  Tags = "Pancreatic cancer expression",
  ResourceName = c(paste0(infoString, ".rda"))
)

write.csv(file="metadata.csv", main.data, row.names=FALSE)

