#' Extraction of cooperation response genes and differentially expressed genes based on transcriptomic profile
#'
#' This function allows you to extract cooperation response genes and differentially expressed genes.
#' @export
#' @importFrom qpgraph qpNrr
#' @importFrom qpgraph qpHist
#' @importFrom qpgraph qpGraphDensity
#' @importFrom qpgraph qpGraph
#' @importFrom igraph as_graphnel
#' @importFrom igraph graph.data.frame
#' @import MASS
#' @import pvclust
#' @import RedeR
#' @import bnlearn
#' @import gRain
#' @import ggplot2
#' @import gbm
#' @import limma
#' @import DMwR
#' @examples 
#' CRGs_extraction()

CRGs_extraction <- NULL
CRGs_extraction <- function() {
  
  data(Transcriptomics)
  
  #Load the text file in Destiny_Folder
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
  
  iP53 <- as.numeric(rep(c(0,0,1,1),c(n1,n2,n3,n4)))
  iRas <- as.numeric(rep(c(0,1,0,1),c(n1,n2,n3,n4)))
  
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
    as.numeric(rownames(res)[res$adj.P.Val < 0.01])
  iDiffExprs <-
    as.numeric(rownames(res2)[res2$adj.P.Val < 0.01])
  ind <- intersect(iInteraction, iDiffExprs)
  newCRGs <- unique(sort(gnames[ind]))
  
  selected <- NULL
  selected <- data.frame()
  n <- NULL
  
  for (n in 1:length(newCRGs))
  {
    selected[n,1] <- n;
    selected[n,2] <- newCRGs[n];
  }
  
  s <- n #number of CRGs
  
  DEG <- NULL
  DEG <- data.frame()
  
  newDEGs <- unique(sort(gnames[iInteraction]))
  
  for (n in 1:length(iInteraction))
  {
    DEG[n,1] <- newDEGs[n];
    DEG[n,2] <-
      log10(mean(as.numeric(e2[iInteraction[n],which(samplename == "q1")] / e2[iInteraction[n],which(samplename == "c")])))
    DEG[n,3] <-
      log10(mean(as.numeric(e2[iInteraction[n],which(samplename == "q2")] / e2[iInteraction[n],which(samplename == "c")])))
    DEG[n,4] <-
      log10(mean(as.numeric(e2[iInteraction[n],which(samplename == "q3")] / e2[iInteraction[n],which(samplename == "c")])))
  }
  
  colnames(DEG) <- c("Gene_names","S1","S2","S1 & S2")
  DEG <- DEG[!duplicated(DEG$Gene_names),]
  DEG <- DEG[!is.na(DEG$Gene_names),]
  rownames(DEG) <- DEG$Gene_names
  
  Destiny_Folder = system.file(package = "NetCRG")
  Destiny_Folder = paste(Destiny_Folder, "/List_CRGs.txt", sep = "")
  
  write.table(
    selected, Destiny_Folder, row.names = FALSE, quote = FALSE, sep = "\t", col.names = FALSE
  )
  
  Destiny_Folder = system.file(package = "NetCRG")
  Destiny_Folder = paste(Destiny_Folder, "/List_DEGs.txt", sep = "")
  
  write.table(
    DEG,Destiny_Folder, row.names = FALSE, quote = FALSE, sep = "\t"
  )
}

#' Use the predefinded list of CRGs for Network inference
#'
#' This function allows you to use the list of known cooperation response genes as input.
#' @export
#' @examples 
#' PL()

PL <- NULL
PL <- function()
  {
  #Load data
  data(Transcriptomics)
  
  reading_expresson_profile <- as.matrix(Transcriptomics)

  samplename <- as.character(reading_expresson_profile[1,])
  
  colnames(reading_expresson_profile) <- samplename
  
  reading_expresson_profile <- reading_expresson_profile[-1,]
  
  n1 <- (which(samplename == "c"))
  
  e <- reading_expresson_profile[,-c(n1)]
 
  #Load data
  data(Differentially_expressed_genes)
  DEGt <- (Differentially_expressed_genes)
  
  colnames <- as.matrix(DEGt[1,])
  DEGt <- DEGt[-1,]
  
  DEG <- data.frame()
  
  for(num_row in 1:dim(DEGt)[1])
  {
    DEG[num_row,1] <- as.character(levels(DEGt[num_row,1])[DEGt[num_row,1]])
    DEG[num_row,2] <- as.numeric(levels(DEGt[num_row,2])[DEGt[num_row,2]])
    DEG[num_row,3] <- as.numeric(levels(DEGt[num_row,3])[DEGt[num_row,3]])
    DEG[num_row,4] <- as.numeric(levels(DEGt[num_row,4])[DEGt[num_row,4]])
  }

  colnames(DEG) <- colnames
  rownames(DEG) <- DEG[,1]
  
  DEG_Gene_names<-c(match(e[,1],DEG$Gene_names))
  DEG_Gene_names<-DEG_Gene_names[!is.na(DEG_Gene_names)]
  DEG <- DEG[DEG_Gene_names,]

  selected <- data.frame()
  CRM_CRG <- data.frame()
  
  #Load data
  data(PLCRG)
  list_CRGs <- (PLCRG)
  
  colnames(list_CRGs) <- as.matrix(list_CRGs[1,])
  list_CRGs <- list_CRGs[-1,]
  
  inf <- 1;
  CC <- 1;
  
  for (inf in 1:length(list_CRGs[,1]))
  {
    a1 <- length(intersect(DEG[,1],list_CRGs$Gene_symbol[inf]))
    a2 <- length(intersect(e[,1],list_CRGs$Gene_symbol[inf]))
    
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
  
  #-- Inference of Network topology -- 
  
  names <- as.matrix(e[,1])
   
  e <- e[,-1]
  e <- apply(e, c(1, 2), as.numeric)
  colnames(e) <- samplename[-c(1,n1)]
  e1 <- matrix(,nrow = dim(e)[1],ncol = dim(e)[2])
  e[is.na(e)] <- 0
  
  for (n in 1:dim(e)[1])
    {
      for (n1 in 1:dim(e)[2])
      {
        if (0 < e[n,n1]) {
          e1[n,n1] <-
            log10(e[n,n1])
        } else {
          e[n,n1] <- -(e[n,n1]);e1[n,n1] <- -log10(e[n,n1])
        }
      }
    }
  
  #calculation of non-rejection rates [qpgraph]
  nrr.estimates <- qpNrr(e1, q = 3)
  #--
  
  #-- Selection of most significant relations --
  interactions_base_NRR <-
    qpGraph(nrr.estimates, threshold = 0.1, return.type = "edge.list") #output is list of
  #---

  interactions_base_NRR2 <-
    matrix(
      interactions_base_NRR,nrow = dim(interactions_base_NRR)[1],ncol = 2
    )
  colnames(interactions_base_NRR2) <- c("geneSymbol", "geneSymbol2")
  
  Fcol <- c(interactions_base_NRR2[,1])
  Scol <- c(interactions_base_NRR2[,2])
  
  Fcol <- type.convert(Fcol, numerals = c("allow.loss"))
  Scol <- type.convert(Scol, numerals = c("allow.loss"))
  
  matrix_of_interactions <- NULL
  matrix_of_interactions <- data.frame(names[Fcol,1],names[Scol,1])
  
  #-- Linear Discriminant Analysis
  
  CRG_Set <- data.frame(); i <- 1; n <- 1;
  
  for (n in 1:s)
    #for each CRG
  {
    #create three set for CRG
    
    num_row <- c(match(selected[n,2], rownames(DEG)));
    
    if (!identical(NULL, num_row))
    {
      S_p53 <- DEG[num_row,2]
      S_ras <- DEG[num_row,3]
      S_p53_ras <- DEG[num_row,4]
      
      CRG_Set[i,1] <-
        selected[n,2]; CRG_Set[i,2] <-
        DEG[num_row,2];CRG_Set[i,3] <- c("S_p53"); i <- i + 1;
      CRG_Set[i,1] <-
        selected[n,2]; CRG_Set[i,2] <-
        DEG[num_row,3];CRG_Set[i,3] <- c("S_ras"); i <- i + 1;
      CRG_Set[i,1] <-
        selected[n,2]; CRG_Set[i,2] <-
        DEG[num_row,4];CRG_Set[i,3] <- c("S_p53_ras"); i <- i + 1;
    }
  }
  colnames(CRG_Set) <- c("gene_symbol", "FC_value", "Set")
  Set <- NULL
  ## visualizations
  S_p53 <- subset(CRG_Set, Set == "S_p53");
  S_ras <- subset(CRG_Set, Set == "S_ras");
  S_p53_ras <- subset(CRG_Set, Set == "S_p53_ras");
  
  par(mfrow = c(1,1))
  boxplot(
    DEG[,2],DEG[,3],DEG[,4], col = (c("gold","darkgreen","yellow")),main = "FC of all DEGs", names =
      c("p53","Ras","p53_Ras"),ylab = "FC"
  )
  boxplot(
    S_p53$FC_value,S_ras$FC_value,S_p53_ras$FC_value, col = (c("gold","darkgreen","yellow")),main =
      "FC of selected CRGs", names = c("p53","Ras","p53_Ras"),ylab = "FC"
  )
  
  summary(S_p53_ras$FC_value)
  summary(DEG[,4])
  
  #Training of Discriminant Analysis model [MASS]
  lda.fit = lda(Set ~ FC_value,data = CRG_Set)
  plot(lda.fit, dimen = 1, type = "both") # fit from lda
  
  #-- extraction of interactions related to each CRGs --
  
  n <- 1; i <- 1; i2 <- 1; p <- 1;
  
  actors <- data.frame()
  
  relations <- data.frame()
  prediction <- data.frame();
  network_all <- data.frame()
  
  recordimprovement <- data.frame();
  rm <- 1
  
  for (n in 1:s)
    #for each CRG
  {
    p <- 1;
    
    num_row <- c(match(selected[n,2], rownames(DEG)));
    
    S_p53 <- DEG[num_row,2]
    S_ras <- DEG[num_row,3]
    S_p53_ras <- DEG[num_row,4]
    
    actors[i,1] <- (name = selected[n,2]); i <- i + 1;
    
    #The FC of given CGR (s) in three groups of experiments
    
    CT1 <- c(S_p53_ras)
    CT2 <- c(S_p53)
    CT3 <- c(S_ras)
    
    #List of genes that interact directly with given CGRs
    
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
        
        ST1 <- c(CT1, S_p53_ras)
        ST2 <- c(CT2, S_p53)
        ST3 <- c(CT3, S_ras)
        
        CT1 <- c(S_p53_ras, CT1)
        CT2 <- c(S_p53, CT2)
        CT3 <- c(S_ras, CT3)
        
        prediction[p,1] <-
          DEG[num_row,1]; prediction[p,2] <-
          S_p53;prediction[p,3] <- c("S_p53"); p <- p + 1;
        prediction[p,1] <-
          DEG[num_row,1]; prediction[p,2] <-
          S_ras;prediction[p,3] <- c("S_ras"); p <- p + 1;
        prediction[p,1] <-
          DEG[num_row,1]; prediction[p,2] <-
          S_p53_ras;prediction[p,3] <- c("S_p53_ras"); p <- p + 1;
        
        rm(ST1, ST2, ST3)
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
        
        ST1 <- c(CT1, S_p53_ras)
        ST2 <- c(CT2, S_p53)
        ST3 <- c(CT3, S_ras)
        
        CT1 <- c(S_p53_ras, CT1)
        CT3 <- c(S_p53, CT2)
        CT2 <- c(S_ras, CT3)
        
        prediction[p,1] <-
          DEG[num_row,1]; prediction[p,2] <-
          S_p53;prediction[p,3] <- c("S_p53"); p <- p + 1;
        prediction[p,1] <-
          DEG[num_row,1]; prediction[p,2] <-
          S_ras;prediction[p,3] <- c("S_ras"); p <- p + 1;
        prediction[p,1] <-
          DEG[num_row,1]; prediction[p,2] <-
          S_p53_ras;prediction[p,3] <- c("S_p53_ras"); p <- p + 1;
        
        rm(ST1, ST2, ST3)
      }
    }
    
    #-- The score of extracted set of given CRG + its related DEGs
    
    if (length(prediction))
    {
      #use the FC values of related genes to the given CRG as test set
      colnames(prediction) <- c("gene_symbol", "FC_value", "Set")
      lda.pred = predict(lda.fit,prediction)
      ct <-
        table(lda.pred$class,prediction$Set)#table() returns a contingency table
      mean(lda.pred$class == prediction$Set)
      diag(prop.table(ct, 1)) # percent correct for each category of G
      
      record_first <-
        sum(diag(prop.table(ct))) # total percent correct
      #save the results of prediction
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
        #print(n1)
        prediction_test <- prediction_final[-(n1:(n1 + 2)),]
        #head(prediction_test)
        
        colnames(prediction_test) <-
          c("gene_symbol", "FC_value", "Set")
        #prediction_test=na.omit(prediction_test)
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
          #save the given CRG and its related gene as graph
          
          actors[i,1] <-
            as.character(prediction_final$gene_symbol[n1]);i <- i + 1;
            
            relations[i2,1] <- selected[n,2];
            relations[i2,2] <-
              prediction_final$gene_symbol[n1]; i2 <- i2 + 1;
              
        }
        t <- t + 1; #for saving total percent correct
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
      
      #Boxplot all CRGs
      S_p53_3 <- subset(CRG_Set, Set == "S_p53")
      S_ras_3 <- subset(CRG_Set, Set == "S_ras")
      S_p53_ras_3 <- subset(CRG_Set, Set == "S_p53_ras")
      
      op <- par(mfrow = c(3,1))
      # boxplot(DEG[,4],DEG[,3],DEG[,2], col=(c("gold","darkgreen","yellow")),
      # main="Distribution of FCs in training sets - DEGs", names=c("p53_Ras","Ras","p53"),ylab ="FC")
      boxplot(
        S_p53_ras_3$FC_value,S_ras_3$FC_value,S_p53_3$FC_value, col = (c("gold","darkgreen","yellow")),
        main = "Distribution of FCs in training sets", names = c("p53_Ras","Ras","p53"),ylab =
          "FC"
      )
      boxplot(
        S_p53_ras$FC_value,S_ras$FC_value,S_p53$FC_value, col = (c("gold","darkgreen","yellow")),
        main = "Distribution of FCs in test sets before Refinement",sub =
          selected[n,2], names = c("p53_Ras","Ras","p53"),ylab = "FC",cex = 0.1
      )
      boxplot(
        S_p53_ras_2$FC_value,S_ras_2$FC_value,S_p53_2$FC_value, col = (c("gold","darkgreen","yellow")),
        main = "Distribution of FCs in test sets after Refinement",sub =
          selected[n,2], names = c("p53_Ras","Ras","p53"),ylab = "FC"
      )
      
      on.exit(par(op))
      
    }
    rm(CT1, CT2, CT3)
  }
  
  #visulization of results
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
  
  # covert relation of CRGs and DEGs to igraph format (the direction of archs are arbitrary)
  actors <- unique(actors)
  actors <- actors[!is.na(actors)]
  relations <- unique(relations)
  relations2 <- na.omit(relations)
  colnames(relations2) <- c("from","to")
  g <-
    graph.data.frame(relations2, directed = TRUE, vertices = actors) #igraph
  #---
  
  Destiny_Folder <- system.file(package = "NetCRG")
  Destiny_Folder = paste(Destiny_Folder, "/Topology_of_integrated_network.txt", sep = "")
  
  write.table(
    relations2, Destiny_Folder, sep = "\t", row.names = F, quote = FALSE
  )
  
  #---
  rdp <- RedPort ()
  calld(rdp)
  addGraph(rdp, g, layout.kamada.kawai(g))
  relax(rdp,p2 = 400,p5 = 30,ps = TRUE)
  
}

#' Extraction of cooperation response genes and Network inference
#'
#' This function allows you to extract cooperation response genes and inference related the gene regulatory network.
#' @export
#' @examples 
#' CRGs_extraction_inf()

CRGs_extraction_inf <- NULL
CRGs_extraction_inf <- function() {
  
  data(Transcriptomics)
  
  #Load the text file in Destiny_Folder
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
  
  iP53 <- as.numeric(rep(c(0,0,1,1),c(n1,n2,n3,n4)))
  iRas <- as.numeric(rep(c(0,1,0,1),c(n1,n2,n3,n4)))
  
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
    as.numeric(rownames(res)[res$adj.P.Val < 0.01])
  iDiffExprs <-
    as.numeric(rownames(res2)[res2$adj.P.Val < 0.01])
  ind <- intersect(iInteraction, iDiffExprs)
  newCRGs <- unique(sort(gnames[ind]))
  
  selected <- NULL
  selected <- data.frame()
  n <- NULL
  
  for (n in 1:length(newCRGs))
  {
    selected[n,1] <- n;
    selected[n,2] <- newCRGs[n];
  }
  
  s <- n #number of CRGs
  
  DEG <- NULL
  DEG <- data.frame()
  
  newDEGs <- unique(sort(gnames[iInteraction]))
  
  for (n in 1:length(iInteraction))
  {
    DEG[n,1] <- newDEGs[n];
    DEG[n,2] <-
      log10(mean(as.numeric(e2[iInteraction[n],which(samplename == "q1")] / e2[iInteraction[n],which(samplename == "c")])))
    DEG[n,3] <-
      log10(mean(as.numeric(e2[iInteraction[n],which(samplename == "q2")] / e2[iInteraction[n],which(samplename == "c")])))
    DEG[n,4] <-
      log10(mean(as.numeric(e2[iInteraction[n],which(samplename == "q3")] / e2[iInteraction[n],which(samplename == "c")])))
  }
  
  colnames(DEG) <- c("Gene_names","dRas","dP53","dBoth")
  DEG <- DEG[!duplicated(DEG$Gene_names),]
  DEG <- DEG[!is.na(DEG$Gene_names),]
  rownames(DEG) <- DEG$Gene_names
  
  #-- Inference of Network topology -- 
  
  data(Transcriptomics)
  
  #Load the text file in Destiny_Folder
  reading_expresson_profile <- as.matrix(Transcriptomics)
  gnames <- reading_expresson_profile[,1]
  samplename <- as.character(reading_expresson_profile[1,])
  colnames(reading_expresson_profile) <- samplename
  reading_expresson_profile <- reading_expresson_profile[-1,]
  n1 <- (which(samplename == "c"))
  e <- as.data.frame(reading_expresson_profile[,-c(n1)])
  gnames <- as.matrix(e[,1])
  e <- e[,-1]
  e <- apply(e, c(1, 2), as.numeric)
  colnames(e) <- samplename[-c(1,n1)]
  
  DEG_Gene_names <- c(match(DEG$Gene_names,gnames))
  DEG_Gene_names <- DEG_Gene_names[!is.na(DEG_Gene_names)]
  e <- e[DEG_Gene_names,]
  
  names <- (gnames)
   
  e1 <- matrix(,nrow = dim(e)[1],ncol = dim(e)[2])
  
  e[is.na(e)] <- 0
  
  for (n in 1:dim(e)[1])
  {
    for (n1 in 1:dim(e)[2])
    {
      if (0 < e[n,n1]) {
        e1[n,n1] <-
          log10(e[n,n1])
      } else {
        e[n,n1] <- -(e[n,n1]);e1[n,n1] <- -log10(e[n,n1])
      }
    }
  }
  
  #calculation of non-rejection rates [qpgraph]
  nrr.estimates <- qpNrr(e1, q = 3)
  qpHist(nrr.estimates, A = NULL,titlehist = "all estimated\n\tnon-rejection rates", freq =
           TRUE)
  abline(v = 0.1, col = "red", lty = 3)
  #--
  
  #-- Selection of most significant relations --

  interactions_base_NRR <-
    qpGraph(nrr.estimates, threshold = 0.1, return.type = "edge.list") #output is list of
  
  #---
  
  interactions_base_NRR2 <-
    matrix(
      interactions_base_NRR,nrow = dim(interactions_base_NRR)[1],ncol = 2
    )
  colnames(interactions_base_NRR2) <- c("geneSymbol", "geneSymbol2")
  
  Fcol <- c(interactions_base_NRR2[,1])
  Scol <- c(interactions_base_NRR2[,2])
  
  Fcol <- type.convert(Fcol, numerals = c("allow.loss"))
  Scol <- type.convert(Scol, numerals = c("allow.loss"))
  
  matrix_of_interactions <- NULL
  matrix_of_interactions <- data.frame(names[Fcol,1],names[Scol,1])
  
  #-- Linear Discriminant Analysis
  
  CRG_Set <- data.frame(); i <- 1; n <- 1;
  
  for (n in 1:s)
    #for each CRG
  {
    #create three set for CRG
    
    num_row <- c(match(selected[n,2], rownames(DEG)));
    
    if (!identical(NULL, num_row))
    {
      S_p53 <- DEG[num_row,2]
      S_ras <- DEG[num_row,3]
      S_p53_ras <- DEG[num_row,4]
      
      CRG_Set[i,1] <-
        selected[n,2]; CRG_Set[i,2] <-
        DEG[num_row,2];CRG_Set[i,3] <- c("S_p53"); i <- i + 1;
      CRG_Set[i,1] <-
        selected[n,2]; CRG_Set[i,2] <-
        DEG[num_row,3];CRG_Set[i,3] <- c("S_ras"); i <- i + 1;
      CRG_Set[i,1] <-
        selected[n,2]; CRG_Set[i,2] <-
        DEG[num_row,4];CRG_Set[i,3] <- c("S_p53_ras"); i <- i + 1;
    }
  }
  colnames(CRG_Set) <- c("gene_symbol", "FC_value", "Set")
  Set <- NULL
  ## visualizations
  S_p53 <- subset(CRG_Set, Set == "S_p53");
  S_ras <- subset(CRG_Set, Set == "S_ras");
  S_p53_ras <- subset(CRG_Set, Set == "S_p53_ras");
  
  par(mfrow = c(1,1))
  boxplot(
    DEG[,2],DEG[,3],DEG[,4], col = (c("gold","darkgreen","yellow")),main = "FC of all DEGs", names =
      c("p53","Ras","p53_Ras"),ylab = "FC"
  )
  boxplot(
    S_p53$FC_value,S_ras$FC_value,S_p53_ras$FC_value, col = (c("gold","darkgreen","yellow")),main =
      "FC of selected CRGs", names = c("p53","Ras","p53_Ras"),ylab = "FC"
  )
  
  summary(S_p53_ras$FC_value)
  summary(DEG[,4])
  
  #Training of Discriminant Analysis model [MASS]
  lda.fit = lda(Set ~ FC_value,data = CRG_Set)
  plot(lda.fit, dimen = 1, type = "both") # fit from lda
  
  #-- extraction of interactions related to each CRGs --
  
  n <- 1; i <- 1; i2 <- 1; p <- 1;
  
  actors <- data.frame()
  
  relations <- data.frame()
  prediction <- data.frame();
  network_all <- data.frame()
  
  recordimprovement <- data.frame();
  rm <- 1
  
  for (n in 1:s)
    #for each CRG
  {
    p <- 1;
    
    num_row <- c(match(selected[n,2], rownames(DEG)));
    
    S_p53 <- DEG[num_row,2]
    S_ras <- DEG[num_row,3]
    S_p53_ras <- DEG[num_row,4]
    
    actors[i,1] <- (name = selected[n,2]); i <- i + 1;
    
    #The FC of given CGR (s) in three groups of experiments
    
    CT1 <- c(S_p53_ras)
    CT2 <- c(S_p53)
    CT3 <- c(S_ras)
    
    #List of genes that interact directly with given CGRs
    
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
        
        ST1 <- c(CT1, S_p53_ras)
        ST2 <- c(CT2, S_p53)
        ST3 <- c(CT3, S_ras)
        
        CT1 <- c(S_p53_ras, CT1)
        CT2 <- c(S_p53, CT2)
        CT3 <- c(S_ras, CT3)
        
        prediction[p,1] <-
          DEG[num_row,1]; prediction[p,2] <-
          S_p53;prediction[p,3] <- c("S_p53"); p <- p + 1;
        prediction[p,1] <-
          DEG[num_row,1]; prediction[p,2] <-
          S_ras;prediction[p,3] <- c("S_ras"); p <- p + 1;
        prediction[p,1] <-
          DEG[num_row,1]; prediction[p,2] <-
          S_p53_ras;prediction[p,3] <- c("S_p53_ras"); p <- p + 1;
        
        rm(ST1, ST2, ST3)
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
        
        ST1 <- c(CT1, S_p53_ras)
        ST2 <- c(CT2, S_p53)
        ST3 <- c(CT3, S_ras)
        
        CT1 <- c(S_p53_ras, CT1)
        CT3 <- c(S_p53, CT2)
        CT2 <- c(S_ras, CT3)
        
        prediction[p,1] <-
          DEG[num_row,1]; prediction[p,2] <-
          S_p53;prediction[p,3] <- c("S_p53"); p <- p + 1;
        prediction[p,1] <-
          DEG[num_row,1]; prediction[p,2] <-
          S_ras;prediction[p,3] <- c("S_ras"); p <- p + 1;
        prediction[p,1] <-
          DEG[num_row,1]; prediction[p,2] <-
          S_p53_ras;prediction[p,3] <- c("S_p53_ras"); p <- p + 1;
        
        rm(ST1, ST2, ST3)
      }
    }
    
    #-- The score of extracted set of given CRG + its related DEGs
    
    if (length(prediction))
    {
      #use the FC values of related genes to the given CRG as test set
      colnames(prediction) <- c("gene_symbol", "FC_value", "Set")
      lda.pred = predict(lda.fit,prediction)
      ct <-
        table(lda.pred$class,prediction$Set)#table() returns a contingency table
      mean(lda.pred$class == prediction$Set)
      diag(prop.table(ct, 1)) # percent correct for each category of G
      
      record_first <-
        sum(diag(prop.table(ct))) # total percent correct
      #save the results of prediction
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
        #print(n1)
        prediction_test <- prediction_final[-(n1:(n1 + 2)),]
        #head(prediction_test)
        
        colnames(prediction_test) <-
          c("gene_symbol", "FC_value", "Set")
        #prediction_test=na.omit(prediction_test)
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
          #save the given CRG and its related gene as graph
          
          actors[i,1] <-
            as.character(prediction_final$gene_symbol[n1]);i <- i + 1;
            
            relations[i2,1] <- selected[n,2];
            relations[i2,2] <-
              prediction_final$gene_symbol[n1]; i2 <- i2 + 1;
              
        }
        t <- t + 1; #for saving total percent correct
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
      
      #Boxplot all CRGs
      S_p53_3 <- subset(CRG_Set, Set == "S_p53")
      S_ras_3 <- subset(CRG_Set, Set == "S_ras")
      S_p53_ras_3 <- subset(CRG_Set, Set == "S_p53_ras")
      
      op <- par(mfrow = c(3,1))
      # boxplot(DEG[,4],DEG[,3],DEG[,2], col=(c("gold","darkgreen","yellow")),
      # main="Distribution of FCs in training sets - DEGs", names=c("p53_Ras","Ras","p53"),ylab ="FC")
      boxplot(
        S_p53_ras_3$FC_value,S_ras_3$FC_value,S_p53_3$FC_value, col = (c("gold","darkgreen","yellow")),
        main = "Distribution of FCs in training sets", names = c("p53_Ras","Ras","p53"),ylab =
          "FC"
      )
      boxplot(
        S_p53_ras$FC_value,S_ras$FC_value,S_p53$FC_value, col = (c("gold","darkgreen","yellow")),
        main = "Distribution of FCs in test sets before Refinement",sub =
          selected[n,2], names = c("p53_Ras","Ras","p53"),ylab = "FC",cex = 0.1
      )
      boxplot(
        S_p53_ras_2$FC_value,S_ras_2$FC_value,S_p53_2$FC_value, col = (c("gold","darkgreen","yellow")),
        main = "Distribution of FCs in test sets after Refinement",sub =
          selected[n,2], names = c("p53_Ras","Ras","p53"),ylab = "FC"
      )
      
      on.exit(par(op))
      
    }
    rm(CT1, CT2, CT3)
  }
  
  #visulization of results
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
  
  # covert relation of CRGs and DEGs to igraph format (the direction of archs are arbitrary)
  actors <- unique(actors)
  actors <- actors[!is.na(actors)]
  relations <- unique(relations)
  relations2 <- na.omit(relations)
  colnames(relations2) <- c("from","to")
  g <-
    graph.data.frame(relations2, directed = TRUE, vertices = actors) #igraph
  #---
  
  Destiny_Folder <- system.file(package = "NetCRG")
  Destiny_Folder = paste(Destiny_Folder, "/Topology_of_integrated_network.txt", sep = "")
  
  write.table(
    relations2, Destiny_Folder, sep = "\t", row.names = F, quote = FALSE
  )
  
  #---
  rdp <- RedPort ()
  calld(rdp)
  addGraph(rdp, g, layout.kamada.kawai(g))
  relax(rdp,p2 = 400,p5 = 30,ps = TRUE)
  
}