#' @export
#' @importFrom igraph as_graphnel
#' @importFrom igraph graph.data.frame
#' @import MASS
#' @import pvclust
#' @import RedeR
#' @import gRain
#' @import ggplot2
#' @import gbm
#' @import limma
#' @import DMwR
#' @import EBcoexpress
#' @import matrixStats
#'
#'@name{SRG}
#'@title{Identification of differentially expressed and synergistic response genes (SRGs) in transcriptomics profile}
#'@description{We use the "SRG" function for identification of synergy or interaction effects between genes in transcriptomics profile in response to combination of mutations, drugs or environmental exposure. SRG returns lists of synergistic response genes and differentially expressed genes, which can be found in home directory of package as text file under title of "List_SRGs" and "List_DEGs", respectively.}
#'@author{Matthew McCall, Isar Nassiri}
#'@examples{
#'data(Transcriptomics)
#'SRG(0.01)}

SRG <- NULL
SRG <- function(pvalue) {

  e2 <- as.matrix(Transcriptomics)

  gnames <- e2[,1]
  samplename <- as.character(e2[1,])
  samplename <- samplename[-1]
  e2 <- e2[-1,-1]

  e2 <- apply(e2, c(1, 2), as.numeric)

  n1 <- length(which(samplename == "c"))
  n2 <- length(which(samplename == "q1"))
  n3 <- length(which(samplename == "q2"))
  n4 <- length(which(samplename == "q3"))

  iP53 <- as.numeric(rep(c(0,1,0,1),c(n1,n2,n3,n4)))
  iRas <- as.numeric(rep(c(0,1,1,0),c(n1,n2,n3,n4)))

  #model.matrix creat the design matix as input of lmFit
  design <- model.matrix( ~ iP53 * iRas)
  colnames(design) <- c("none","dRas","dP53","dBoth")

  #lmFit use to identify which genes are differentially expressed between the double mutant samples and YAMC
  fit <- lmFit(e2, design)
  fit <- eBayes(fit)
  res <- topTable(fit, coef = 4, number = nrow(e2))

  ## which genes have significant interaction term
  dP53 <- NULL
  dRas <- NULL
  dBoth <- NULL
  contrast.matrix <-
    makeContrasts(dP53 + dRas + dBoth, levels = design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  res2 <- topTable(fit2, coef = 1, number = nrow(e2))

  ## which genes are both differentially expressed and show an interaction
  ## could alter the p-value threshold to be more or less restrictive
  iInteraction <-
    as.numeric(rownames(res)[res$adj.P.Val < pvalue])
  iDiffExprs <-
    as.numeric(rownames(res2)[res2$adj.P.Val < pvalue])
  ind <- intersect(iInteraction, iDiffExprs)
  newSRGs <- unique(sort(gnames[ind]))

  selected <- NULL
  selected <- data.frame()
  n <- NULL

  for (n in 1:length(newSRGs))
  {
    selected[n,1] <- n;
    selected[n,2] <- newSRGs[n];
  }

  colnames(selected) <- c("Row","geneSymbol")
  s <- n #number of SRGs

  DEG <- NULL
  DEG <- data.frame()

  newDEGs <- unique(sort(gnames[iInteraction]))

  for (n in 1:length(iDiffExprs))
  {
    DEG[n,1] <- newDEGs[n];

    DEG[n,2] <-
      log10(mean(as.numeric(e2[iDiffExprs[n],which(samplename == "q3")] / e2[iDiffExprs[n],which(samplename == "c")])))

    DEG[n,3] <-
      log10(mean(as.numeric(e2[iDiffExprs[n],which(samplename == "q2")] / e2[iDiffExprs[n],which(samplename == "c")])))

    DEG[n,4] <-
      log10(mean(as.numeric(e2[iDiffExprs[n],which(samplename == "q1")] / e2[iDiffExprs[n],which(samplename == "c")])))
  }

  colnames(DEG) <- c("Gene_names","a","b","a+b")   #in example files: "Gene_names","mp53","Ras","mp53/Ras"
  DEG <- DEG[!duplicated(DEG$Gene_names),]
  DEG <- DEG[!is.na(DEG$Gene_names),]
  rownames(DEG) <- DEG$Gene_names

  #-- save the results --
  Destiny_Folder = system.file(package = "SRGnet")
  Destiny_Folder = paste(Destiny_Folder, "/PLCRG.txt", sep = "")       #List of SRGs

  write.table(
    selected, Destiny_Folder, row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE
  )

  Destiny_Folder = system.file(package = "SRGnet")
  Destiny_Folder = paste(Destiny_Folder, "/Differentially_expressed_genes.txt", sep = "")      #List of DEGs

  write.table(
    DEG,Destiny_Folder, row.names = FALSE, quote = FALSE, sep = "\t"
  )
}

#'@name{SRGnet}
#'@title{Gene regulatory network inference based on predifined list of differentially expressed genes, list of synergistic response genes and transcriptomics profile}
#'@description{The "SRGnet" can be applied if user has transcriptomic profile, list of differentially expressed genes and synergistic response genes as inputs. The function can be ran in two mode of Slow or Fast. In fast mode, step of expectation maximization for estimation of hyperparameters is omitted. User can run the function in fast or slow mode by using the "F" or "S" as input of "PL()" function, respectively [e.g. SRGs_identification(“F”)]. SRGs_identification returns the topology of SRMs network and ranked list of genes in network based on differential connectivity score, which can be found in home directory of package under title of "DC_score" and "Topology_of_integrated_network" as text files.}
#'@author{Isar Nassiri, Matthew McCall}
#'@examples{
#'data(Differentially_expressed_genes)
#'data(Transcriptomics)
#'data(PLCRG)
#'SRGnet("F")
#'}
#'@export

SRGnet <- NULL
SRGnet <- function(type_of_run) #SRGnet("F"): fast run; SRGnet("S"): slow run
{
  #Load data
  
  DEGt <- Differentially_expressed_genes
  
  colnames <- as.matrix(DEGt[1,])
  DEGt <- DEGt[-1,]
    
  DEG <- data.frame()
  
  for (num_row in 1:dim(DEGt)[1])
  {
    DEG[num_row,1] <-
      as.character(levels(DEGt[num_row,1])[DEGt[num_row,1]])
    DEG[num_row,2] <-
      as.numeric(levels(DEGt[num_row,2])[DEGt[num_row,2]])
    DEG[num_row,3] <-
      as.numeric(levels(DEGt[num_row,3])[DEGt[num_row,3]])
    DEG[num_row,4] <-
      as.numeric(levels(DEGt[num_row,4])[DEGt[num_row,4]])
  }
  
  colnames(DEG) <- colnames
  rownames(DEG) <- DEG[,1]
  #-----
  #Load data
 
  reading_expresson_profile <- as.matrix(Transcriptomics)
  
  samplename <- as.character(reading_expresson_profile[1,])
  
  colnames(reading_expresson_profile) <- samplename
  
  reading_expresson_profile <- reading_expresson_profile[-1,]
  
  samplename <- samplename[-1]
  
  e <- reading_expresson_profile
  
  DEG_Gene_names <- c(match(DEG[,1],e[,1]))
  DEG_Gene_names <- DEG_Gene_names[!is.na(DEG_Gene_names)]
  
  e <- e[DEG_Gene_names,]
  
  names <- as.matrix(e[,1])
  
  e <- e[,-1]
  e <- apply(e, c(1, 2), as.numeric)
  e1 <- matrix(,nrow = dim(e)[1],ncol = dim(e)[2])
  e[is.na(e)] <- 0
  
  for (n in 1:dim(e)[1])
  {
    for (n1 in 1:dim(e)[2])
    {
      e1[n,n1] <- log10(e[n,n1])
    }
  }
  
  rownames(e1) <- names
  colnames(e1) <- samplename
  
  #-------------
 
  selected <- data.frame()
  CRM_CRG <- data.frame()
    
  list_CRGs <- (PLCRG)
  
  colnames(list_CRGs) <- as.matrix(list_CRGs[1,])
  list_CRGs <- list_CRGs[-1,]
  
  inf <- 1;
  CC <- 1;
   
  for (inf in 1:length(list_CRGs[,1]))
  {
    a1 <- length(intersect(DEG[,1],list_CRGs$Gene_symbol[inf]))
    a2 <- length(intersect(names,list_CRGs$Gene_symbol[inf]))
    
    if (1 == (a1) && 1 == (a2))
    {
      selected[CC,1] <- list_CRGs[inf,2]
      selected[CC,2] <- list_CRGs[inf,1]
      CC = CC + 1;
      a1 <- 0
    }
  }
  
  colnames(selected) <- c("Synergy_score","geneSymbol")
  s <- length(selected$geneSymbol)
  
  #-- DC analysis --
  n1 <- length(grep("c", samplename))
  tinyCond <- as.numeric(rep(c(1,2),c(n1,length(samplename) - n1)))
  tinyPat <- ebPatterns(c("1,1","1,2"))
  
  D <- makeMyD(e1, tinyCond, useBWMC = F, gpsep = "~")
  
  set.seed(3)
  initHP <- initializeHP(D[,1], tinyCond)
  
  if (type_of_run == "S")
  {
    fout <-
      ebCoexpressOneStep(D, tinyCond, tinyPat, initHP, controlOptions = list())
  } else {
    fout <-
      ebCoexpressZeroStep(D, tinyCond, tinyPat, initHP, controlOptions = list())
  }
  
  softThresh <- crit.fun(fout$POSTPROB[,1], 0.05)
  ppbDC <- 1 - fout$POSTPROB[,1]
  kept_s2 <- ppbDC[ppbDC >= softThresh]
  
  graph_e <- data.frame()
  
  for (i in 1:length(names(kept_s2)))
  {
    e_ <- unlist(strsplit(names(kept_s2)[i], "~"))
    graph_e[i,1] <- e_[1]
    graph_e[i,2] <- e_[2]
  }
  
  interactions_base_NRR2 <- as.matrix(graph_e)
  colnames(interactions_base_NRR2) <- c("geneSymbol", "geneSymbol2")
  
  #-------------
  Fcol <- c(interactions_base_NRR2[,1])
  Scol <- c(interactions_base_NRR2[,2])
  
  Fcol <- type.convert(Fcol, numerals = c("allow.loss"))
  Scol <- type.convert(Scol, numerals = c("allow.loss"))
  
  matrix_of_interactions <- NULL
  matrix_of_interactions <- data.frame(names[Fcol,1],names[Scol,1])
  
  #-- LDA Analysis --
  
  SRG_Set <- data.frame(); i <- 1; n <- 1;
  
  for (n in 1:s)
  #for each SRG
  {
    #create three set
    
    num_row <- c(match(selected[n,2], rownames(DEG)));
    
    if (!identical(NULL, num_row))
    {
      S_p53 <- DEG[num_row,2]
      S_ras <- DEG[num_row,3]
      S_p53_ras <- DEG[num_row,4]
      
      SRG_Set[i,1] <-
        selected[n,2]; SRG_Set[i,2] <-
        DEG[num_row,2];SRG_Set[i,3] <- c("S_p53"); i <- i + 1;
      SRG_Set[i,1] <-
        selected[n,2]; SRG_Set[i,2] <-
        DEG[num_row,3];SRG_Set[i,3] <- c("S_ras"); i <- i + 1;
      SRG_Set[i,1] <-
        selected[n,2]; SRG_Set[i,2] <-
        DEG[num_row,4];SRG_Set[i,3] <- c("S_p53_ras"); i <- i + 1;
    }
  }
  colnames(SRG_Set) <- c("gene_symbol", "FC_value", "Set")
  Set <- NULL
  
  S_p53 <- subset(SRG_Set, Set == "S_p53");
  S_ras <- subset(SRG_Set, Set == "S_ras");
  S_p53_ras <- subset(SRG_Set, Set == "S_p53_ras");
  
  par(mfrow = c(1,1))
  boxplot(
    DEG[,2],DEG[,3],DEG[,4], col = (c("gold","darkgreen","yellow")),main = "FC of all DEGs", names =
      c("p53","Ras","p53_Ras"),ylab = "FC"
  )
  boxplot(
    S_p53$FC_value,S_ras$FC_value,S_p53_ras$FC_value, col = (c("gold","darkgreen","yellow")),main =
      "FC of selected SRGs", names = c("p53","Ras","p53_Ras"),ylab = "FC"
  )
  
  summary(S_p53_ras$FC_value)
  summary(DEG[,4])
  
  #Training of LDA model [MASS package]
  lda.fit = lda(Set ~ FC_value,data = SRG_Set)
  plot(lda.fit, dimen = 1, type = "both") # fit from lda
  
  #-- extraction of interactions related to each SRGs --
  
  n <- 1; i <- 1; i2 <- 1; p <- 1;
  
  actors <- data.frame()
  
  relations <- data.frame()
  prediction <- data.frame();
  network_all <- data.frame()
  
  recordimprovement <- data.frame();
  rm <- 1
  
  for (n in 1:s)
    #for each SRG
  {
    p <- 1;
    
    num_row <- c(match(selected[n,2], rownames(DEG)));
    
    S_p53 <- DEG[num_row,2]
    S_ras <- DEG[num_row,3]
    S_p53_ras <- DEG[num_row,4]
    
    actors[i,1] <- (name = selected[n,2]); i <- i + 1;
    
    #List of genes that interact directly with indicated SGRs
    
    l <- grep(selected[n,2],matrix_of_interactions[,1])
    L1 <- matrix_of_interactions[l,2]
    l <- grep(selected[n,2],matrix_of_interactions[,2])
    L2 <- matrix_of_interactions[l,1]
    LIN <- length(L1) + length(L2)#number of interaction of s
    #--
    
    if (0 < length(L1))
    {
      for (m in 1:length(L1))
      {
        num_row <- c(match(L1[m], DEG[,1]))
        
        S_p53 <- DEG[num_row,2]
        S_ras <- DEG[num_row,3]
        S_p53_ras <- DEG[num_row,4]
        
        prediction[p,1] <-
          DEG[num_row,1]; prediction[p,2] <-
          S_p53;prediction[p,3] <- c("S_p53"); p <- p + 1;
        prediction[p,1] <-
          DEG[num_row,1]; prediction[p,2] <-
          S_ras;prediction[p,3] <- c("S_ras"); p <- p + 1;
        prediction[p,1] <-
          DEG[num_row,1]; prediction[p,2] <-
          S_p53_ras;prediction[p,3] <- c("S_p53_ras"); p <- p + 1;
      }
    }
    
    if (0 < length(L2))
    {
      for (m in 1:length(L2))
      {
        num_row <- c(match(L2[m], DEG[,1]));
        
        S_p53 <- DEG[num_row,2]
        S_ras <- DEG[num_row,3]
        S_p53_ras <- DEG[num_row,4]
        
        prediction[p,1] <-
          DEG[num_row,1]; prediction[p,2] <-
          S_p53;prediction[p,3] <- c("S_p53"); p <- p + 1;
        prediction[p,1] <-
          DEG[num_row,1]; prediction[p,2] <-
          S_ras;prediction[p,3] <- c("S_ras"); p <- p + 1;
        prediction[p,1] <-
          DEG[num_row,1]; prediction[p,2] <-
          S_p53_ras;prediction[p,3] <- c("S_p53_ras"); p <- p + 1;
      }
    }
    
    #-- Refinement step --
    
    if (length(prediction))
    {
      #use the FC values of related genes to the given SRG as test set
      colnames(prediction) <- c("gene_symbol", "FC_value", "Set")
      lda.pred = predict(lda.fit,prediction)
      ct <-
        table(lda.pred$class,prediction$Set)#table() returns a contingency table
      mean(lda.pred$class == prediction$Set)
      diag(prop.table(ct, 1)) # percent correct for each category
      
      record_first <-
        sum(diag(prop.table(ct))) # total percent correct
      recordimprovement[rm,1] <- record_first;
      #--
      
      #Some data structures before refinemene
      record <- data.frame(); n1 <- -2; t <- 1;
      prediction_final <- data.frame();
      prediction_test <- data.frame();
      prediction_final <- prediction
      
      boundary <- (dim(prediction)[1] - 3) / 3
      #--
      
      for (u in 1:boundary)
      {
        n1 <- n1 + 3;
        
        prediction_test <- prediction_final[-(n1:(n1 + 2)),]
        
        colnames(prediction_test) <-
          c("gene_symbol", "FC_value", "Set")
        
        lda.pred = predict(lda.fit,prediction_test)
        ct <-
          table(lda.pred$class,prediction_test$Set)#table() returns a contingency table
        mean(lda.pred$class == prediction_test$Set)
        diag(prop.table(ct, 1)) # percent correct for each category of G
        
        record[t,1] <- sum(diag(prop.table(ct)));
        
        a <- record_first
        b <- record[t,1]
        c <- (a <= b)
        
        if (c)
        {
          prediction_final <- prediction_test;
          n1 <- n1 - 3;
          
        } else {
          actors[i,1] <-
            as.character(prediction_final$gene_symbol[n1]);i <-
              i + 1;
            
            relations[i2,1] <- selected[n,2];
            relations[i2,2] <-
              prediction_final$gene_symbol[n1]; i2 <- i2 + 1;
              
        }
        t <- t + 1; #total percent correct
      }
      
      #Save the results of prediction
      recordimprovement[rm,2] <- max(record);
      rm <- rm + 1;
      
      #Boxplot before refinement
      S_p53 <- subset(prediction, Set == "S_p53")
      S_ras <- subset(prediction, Set == "S_ras")
      S_p53_ras <- subset(prediction, Set == "S_p53_ras")
      
      #Boxplot after refinement
      S_p53_2 <- subset(prediction_final, Set == "S_p53")
      S_ras_2 <- subset(prediction_final, Set == "S_ras")
      S_p53_ras_2 <- subset(prediction_final, Set == "S_p53_ras")
      
      #Boxplot all SRGs
      S_p53_3 <- subset(SRG_Set, Set == "S_p53")
      S_ras_3 <- subset(SRG_Set, Set == "S_ras")
      S_p53_ras_3 <- subset(SRG_Set, Set == "S_p53_ras")
      
      op <- par(mfrow = c(3,1))
      par(mar = c(5.1,4.1,4.1,2.1))
      
      boxplot(
        S_p53_ras_3$FC_value,S_ras_3$FC_value,S_p53_3$FC_value, col = terrain.colors(4),
        main = "Training sets", names = c("p53_Ras","Ras","p53"),ylab =
          "FC", cex = 0.9
      )
      boxplot(
        S_p53_ras$FC_value,S_ras$FC_value,S_p53$FC_value, col = terrain.colors(4),
        main = "Test sets before Refinement",
        xlab = selected[n,2], names = c("p53_Ras","Ras","p53"), ylab = "FC", cex = 0.9
      )
      boxplot(
        S_p53_ras_2$FC_value,S_ras_2$FC_value,S_p53_2$FC_value, col = terrain.colors(4),
        main = "Test sets after Refinement",xlab =
          selected[n,2], names = c("p53_Ras","Ras","p53"), ylab = "FC", cex = 0.9
      )
      on.exit(par(op))
    }
  }
  
  recordimprovement$V1
  op <- par(mfrow = c(1,1))
  boxplot(
    recordimprovement$V1 * 100,recordimprovement$V2 * 100, col = (c("light yellow","light green")),
    main = paste("Distribution of total percent correct fit"),names =
      c("Before Refinement","After Refinement"),ylab = "Percentage"
  )
  
  dim(recordimprovement)
  ttest <- t.test(recordimprovement$V1,recordimprovement$V2)
  qt(c(.025, .975), ttest$parameter)
  
  #igraph format
  actors <- unique(actors)
  actors <- actors[!is.na(actors)]
  relations <- unique(relations)
  relations2 <- na.omit(relations)
  colnames(relations2) <- c("from","to")
  g <-
    graph.data.frame(relations2, directed = TRUE, vertices = actors) #igraph
  
  graph.density(g, loops = TRUE)   #self loops are considered to be possible.
  #---
  
  Destiny_Folder <- system.file(package = "SRGnet")
  Destiny_Folder = paste(Destiny_Folder, "/Topology_of_integrated_network.txt", sep = "")
  
  write.table(
    relations2, Destiny_Folder, sep = "\t", row.names = F, quote = FALSE
  )
  
  #--- Visulization of network ---
  #rdp <- RedPort()
  #calld(rdp)
  #addGraph(rdp, g, layout.kamada.kawai(g))
  #relax(rdp,p2 = 400,p5 = 30,ps = TRUE)
  
  #--- ranking SRGs ---
  
  #Load data
  
  reading_expresson_profile <- as.matrix(Transcriptomics)
  
  samplename <- as.character(reading_expresson_profile[1,])
  
  colnames(reading_expresson_profile) <- samplename
  
  reading_expresson_profile <- reading_expresson_profile[-1,]
  
  samplename <- samplename[-1]
  
  e <- reading_expresson_profile
  
  DEG_Gene_names <- c(match(actors,e[,1]))
  DEG_Gene_names <- DEG_Gene_names[!is.na(DEG_Gene_names)]
  
  e <- e[DEG_Gene_names,]

  names <- as.matrix(e[,1])
  
  e <- e[,-1]
  e <- apply(e, c(1, 2), as.numeric)
  e1 <- matrix(,nrow = dim(e)[1],ncol = dim(e)[2])
  e[is.na(e)] <- 0
  
  for (n in 1:dim(e)[1])
  {
    for (n1 in 1:dim(e)[2])
    {
      e1[n,n1] <- log10(e[n,n1])
    }
  }
  
  rownames(e1) <- names
  colnames(e1) <- samplename
  
  #-------------
  n1 <- length(grep("c", samplename))
  
  tinyCond <- as.numeric(rep(c(1,2),c(n1,length(samplename) - n1)))
  tinyPat <- ebPatterns(c("1,1","1,2"))
  
  D <- makeMyD(e1, tinyCond, useBWMC = F, gpsep = "~")
  
  set.seed(3)
  initHP <- initializeHP(D[,1], tinyCond)
  
  if (type_of_run == "S")
  {
    fout <-
      ebCoexpressOneStep(D, tinyCond, tinyPat, initHP, controlOptions = list())
  } else {
    fout <-
      ebCoexpressZeroStep(D, tinyCond, tinyPat, initHP, controlOptions = list())
  }
  
  hubs <-
    rankMyGenes(fout, thresh = 0.95, sep = "~")       #A threshold for determining whether a pair is DC.
  print(
    "Prioritization based on Differential Connectivity (DC) score: Name of gene, degree"
  )
  print(hubs)
  
  hubs <- data.frame(hubs)
  colnames(hubs) <- c("Gene_name DC_score")
  
  Destiny_Folder <- system.file(package = "SRGnet")
  Destiny_Folder = paste(Destiny_Folder, "/DC_score.txt", sep = "")
  
  write.table(
    hubs, Destiny_Folder, sep = "\t", row.names = TRUE, quote = FALSE
  )
  
}