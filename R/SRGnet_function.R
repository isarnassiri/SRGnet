#' @export
#' @import igraph
#' @import MASS
#' @import pvclust
#' @import RedeR
#' @import gRain
#' @import gbm
#' @import limma
#' @import DMwR
#' @import EBcoexpress
#' @import matrixStats
#' @import Hmisc
#'
#'@title{Synergistic response to gene mutations specific network}
#'@description{The "SRGnet" can be applied if user has transcriptomic profile, list of differentially expressed genes and synergistic response genes as inputs. The function can be ran in two mode of Slow or Fast. In fast mode, step of expectation maximization for estimation of hyperparameters is omitted. User can run the function in fast or slow mode by using the "F" or "S" as input of "PL()" function, respectively [e.g. SRGs_identification(“F”)]. SRGs_identification returns the topology of SRMs network and ranked list of genes in network based on differential connectivity score, which can be found in home directory of package under title of "DC_score" and "Topology_of_integrated_network" as text files.}
#'@author{Isar Nassiri, Matthew McCall}
#'@param type_of_run A character, "F": Fast or "S": Slow.
#'@examples{
#'data(Differentially_expressed_genes)
#'data(Transcriptomics)
#'data(PLCRG)
#'SRGnet("F") #Fast run  
#'}
#'@export

SRGnet <- NULL
SRGnet <- function(type_of_run) 
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
   
  e <- as.matrix(Transcriptomics)
  
  samplename <- as.character(colnames(e))

  DEG_Gene_names <- c(match(DEG[,1],rownames(e)))
  DEG_Gene_names <- DEG_Gene_names[!is.na(DEG_Gene_names)]
  
  e <- e[DEG_Gene_names,]
  
  names <- as.matrix(rownames(e))
  
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
  
  samplename <- colnames(e1)
  
  n1 <- length(which(samplename == "c"))
  n2 <- length(which(samplename == "q1"))
  n3 <- length(which(samplename == "q2"))
  n4 <- length(which(samplename == "q3"))
  
  tinyCond <- as.numeric(rep(c(1, 2, 3, 4), c(n1, n2, n3, n4)))
  length(tinyCond)
  
  #-- pattern --
  
  p1 <- paste(c(1, 1, 1, 1), collapse=", ")  #m_c=m_both=m_ras=m_p53
  p2 <- paste(c(1, 2, 0, 0), collapse=", ")
  tinyPat <- ebPatterns(c(p1, p2))
  
  D2 <- makeMyD(e1, tinyCond, useBWMC = F, gpsep = "~")
     
  set.seed(123)
  initHP <- initializeHP(D2[,1], tinyCond, seed = 123)
  
  if (type_of_run == "S")
  {
    fout <-
      ebCoexpressOneStep(D2, tinyCond, tinyPat, initHP, controlOptions = list(applyTransform = TRUE))
  } else {
    fout <-
      ebCoexpressZeroStep(D2, tinyCond, tinyPat, initHP, controlOptions = list(applyTransform = TRUE))
  }
  
  #-- Selection of differential coexpressed pairs of genes based on the FDR
  
  ppbDC <- 1 - fout$POSTPROB[,1]
  kept_h <- ppbDC[ppbDC >= 0.999]
  length(names(kept_h))
  
  graph_e <- data.frame()
  
  for (i in 1:length(names(kept_h)))
  {
    e_ <- unlist(strsplit(names(kept_h)[i], "~"))
    graph_e[i,1] <- e_[1]
    graph_e[i,2] <- e_[2]
  }
  
  matrix_of_interactions <- as.matrix(graph_e)
  colnames(matrix_of_interactions) <- c("geneSymbol", "geneSymbol2")
  print("Number of differential coexpressed pairs of genes: ")
  print(dim(matrix_of_interactions))
  
  #-- LDA Analysis --
  
  SRG_Set <- data.frame(); i <- 1; n <- 1;
  
  for (n in 1:s)
  {
    
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
  {
    p <- 1;
    
    num_row <- c(match(selected[n,2], rownames(DEG)));
    
    S_p53 <- DEG[num_row,2]
    S_ras <- DEG[num_row,3]
    S_p53_ras <- DEG[num_row,4]
    
    actors[i,1] <- (name = selected[n,2]); i <- i + 1;
    
    l <- grep(selected[n,2],matrix_of_interactions[,1])
    L1 <- matrix_of_interactions[l,2]
    l <- grep(selected[n,2],matrix_of_interactions[,2])
    L2 <- matrix_of_interactions[l,1]
    LIN <- length(L1) + length(L2)	
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
      
      record <- data.frame(); n0 <- -2; t <- 1;
      prediction_final <- data.frame();
      prediction_test <- data.frame();
      prediction_final <- prediction
      
      boundary <- (dim(prediction)[1] - 3) / 3
      #--
      
      for (u in 1:boundary)
      {
        n0 <- n0 + 3;
        
        prediction_test <- prediction_final[-(n0:(n0 + 2)),]
        
        colnames(prediction_test) <-
          c("gene_symbol", "FC_value", "Set")
        
        lda.pred = predict(lda.fit,prediction_test)
        ct <-
          table(lda.pred$class,prediction_test$Set)
        mean(lda.pred$class == prediction_test$Set)
        diag(prop.table(ct, 1)) 
        
        record[t,1] <- sum(diag(prop.table(ct)));
        
        a <- record_first
        b <- record[t,1]
        c <- (a <= b)
        
        if (c)
        {
          prediction_final <- prediction_test;
          n0 <- n0 - 3;
          
        } else {
          actors[i,1] <-
            as.character(prediction_final$gene_symbol[n0]);i <-
              i + 1;
            
            relations[i2,1] <- selected[n,2];
            relations[i2,2] <-
              prediction_final$gene_symbol[n0]; i2 <- i2 + 1;
              
        }
        t <- t + 1; #total % correct
      }
      
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
  
  op <- par(mfrow = c(1,1))
  boxplot(
    recordimprovement$V1 * 100,recordimprovement$V2 * 100, col = (c("light yellow","light green")),
    main = paste("Distribution of total percent correct fit"),names =
      c("Before Refinement","After Refinement"),ylab = "Percentage"
  )
  

  actors <- unique(c(relations[,1],relations[,2]))
  relations <- unique(relations)
  relations2 <- na.omit(relations)
  colnames(relations2) <- c("from","to")
  
  #igraph format [use for visulization]
  #g1 <- graph.data.frame(All_CRGs_in_graph_edges, directed=TRUE, vertices=unique(c(All_CRGs_in_graph_edges[,1],All_CRGs_in_graph_edges[,2])))   #igraph, if you conver the g2 to the matrix it will be asymmetrical
  
  #-- Add the interactions between DEGs/DEGs (non-CRGs) 
  
  All_DEGs_in_graph <- setdiff(unique(relations2[,2]), selected[,2])  #column 2 includes the interacting DEGs with CRGs
  length(All_DEGs_in_graph)
  
  for(i in 1:dim(matrix_of_interactions)[1])
  {
    a <- which(matrix_of_interactions[i,1] %in% All_DEGs_in_graph)
    b <- which(matrix_of_interactions[i,2] %in% All_DEGs_in_graph)
    
    if(0<length(a) && 0<length(b))
    {
      relations2 <- rbind(relations2, matrix_of_interactions[i,])
    }
  }
  
  #-- Direction of interactions --
  
  cor <- rcorr(t(e1[,which(samplename == "q1")]), type="pearson")
  R <- as.matrix(cor$r)  #it is same with D
  
  for(i in 1:dim(relations2)[1])
  {
    sign <- R[which(rownames(R) %in% relations2[i,1]), which(rownames(R) %in% relations2[i,2])]
    relations2[i,3] <- ifelse(0<sign, "+", "-")
  }
  
  relations2 <- unique(relations2)
  relations2 <- na.omit(relations2)
  dim(relations2)
  colnames(relations2) <- c("Source_node", "Target_node", "Type_of_interactions")
  
  print("Number of molecules")
  print(length(actors))
  print("Number of interactions")
  print(dim(relations2))
  print("Type of interactions (+: activators, -: inhibitors): ")
  print(table(relations2[,3]))
  
#--- Visulization of network ---

# g2 <- graph.data.frame(relations2[,1:2], directed=TRUE, vertices=actors)   #igraph, if you conver the g2 to the matrix it will be asymmetrical 
# 
# if(visualization)
# {
#   rdp <- RedPort()
#   calld(rdp)
#   
#   addGraph(rdp, g1, gcoord=c(10,25), gscale=20, isNest=TRUE, theme='tm1', zoom=30)
#   addGraph(rdp, g2, gcoord=c(50,70), gscale=50, isNest=TRUE, theme='tm1', zoom=30)
#   
#   nestNodes(rdp, nodes=V(g1)$name, parent="N1", theme='tm2')
#   mergeOutEdges(rdp)
# }
  
  #--- Save the topology of network ---
  
  Destiny_Folder <- system.file(package = "SRGnet")
  Destiny_Folder = paste(Destiny_Folder, "/Topology_of_integrated_network.txt", sep = "")
  
  write.table(
    relations2, Destiny_Folder, sep = "\t", row.names = F, quote = FALSE
  )
  
  print("You can find the results at: ")
  system.file(package="SRGnet")
  
}