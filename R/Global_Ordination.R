#' @importFrom factoextra get_eig
#' @importFrom dplyr %>% select inner_join
#' @importFrom tibble column_to_rownames column_to_rownames
#' @importFrom vegan adonis vegdist
#' @importFrom stats setNames
#' @import convert
#'
#'
#' @title Dimensionality Reduction Analysis: Principal Component Analysis(PCA)
#'
#' @description
#' The common Dimensionality reduction method is PCA. It's suitable for non-zero sparse matrix such as metabolites' profile.
#'
#' @details 12/2/2021 Guangzhou China
#' @author  Hua Zou
#'
#' @param Expression, ExpressionSet; (Required) ExpressionSet object.
#' @param Normalization, Logical; Normalizing feature(default: TRUE).
#' @param Group_info, Character; the group for plot(default: "Group").
#'
#' @return
#' a list object:
#'   PCA score
#'   Result of PERMANOVA
#'
#' @usage Ord_PCA(dataset=ExpressionSet, Normalization=TRUE, Group_info="Group")
#' @examples
#'
#' data(ExprSet_species)
#'
#' PCA_res <- Ord_PCA(dataset=ExprSet_species, Normalization=TRUE, Group_info="Group")
#' PCA_res$PCA
#'
Ord_PCA <- function(dataset=ExprSet_species,
                    Normalization=TRUE,
                    Group_info="Group"){

  metadata <- pData(dataset)
  colnames(metadata)[which(colnames(metadata) == Group_info)] <- "Group"
  profile <- exprs(dataset)

  if(Normalization){
    profile_norm <- apply(t(profile), 1, NormalizeFun, "Zscore")
    rownames(profile_norm) <- rownames(profile)
  }else{
    profile_norm <- profile
  }

  # pca
  pca <- prcomp(t(profile_norm))
  eig <- factoextra::get_eig(pca)
  explains <- paste0(paste0("PC", seq(2)), "(", paste0(round(eig[1:2, 2], 2), "%"), ")")
  score <- dplyr::inner_join(pca$x %>% data.frame() %>%
                        dplyr::select(c(1:2)) %>%
                        stats::setNames(paste0("Axis", seq(2))) %>%
                        tibble::rownames_to_column("SampleID"),
                      metadata %>% tibble::rownames_to_column("SampleID") %>%
                        dplyr::select(all_of(c("SampleID", "Group"))),
                      by = "SampleID")
  # PERMANOVA
  set.seed(123)
  if(any(profile < 0)){
    res_adonis <- vegan::adonis(vegan::vegdist(t(profile), method = "manhattan") ~ metadata$Group, permutations = 999)
  }else{
    res_adonis <- vegan::adonis(vegan::vegdist(t(profile), method = "bray") ~ metadata$Group, permutations = 999)
  }
  adn_pvalue <- res_adonis[[1]][["Pr(>F)"]][1]
  adn_rsquared <- round(res_adonis[[1]][["R2"]][1],3)
  #use the bquote function to format adonis results to be annotated on the ordination plot.
  signi_label <- paste(cut(adn_pvalue,
                           breaks=c(-Inf, 0.001, 0.01, 0.05, Inf),
                           label=c("***", "**", "*", ".")))
  PERMANOVA <- bquote(atop(atop("PERMANOVA", R^2==~.(adn_rsquared)),
                                atop("p-value="~.(adn_pvalue)~.(signi_label), phantom())))

  res <- list(PCA=score,
              epn=explains,
              PER=PERMANOVA)

  return(res)
}


#' @importFrom ape pcoa
#' @importFrom dplyr %>% select inner_join
#' @importFrom tibble column_to_rownames column_to_rownames
#' @importFrom vegan adonis vegdist
#' @importFrom stats setNames
#' @import convert
#'
#'
#' @title Dimensionality Reduction Analysis: Principal Coordinate Analysis(PCoA)
#'
#' @description
#' PCoA uses the distance among samples, which calculates through the multiple variables.
#'
#' @details 12/2/2021 Guangzhou China
#' @author  Hua Zou
#'
#' @param Expression, ExpressionSet; (Required) ExpressionSet object.
#' @param Method, Character; method to calculate the distance(default: bray).
#' @param Group_info, Character; the group for plot(default: "Group").
#'
#' @return
#' a list object:
#'   PCoA score
#'   Result of PERMANOVA
#'
#' @usage Ord_PCoA(dataset=ExpressionSet, Method="bray", Group_info="Group")
#' @examples
#'
#' data(ExprSet_species)
#'
#' PCoA_res <- Ord_PCoA(dataset=ExprSet_species, Method="bray", Group_info="Group")
#' PCoA_res$PCoA
#'
Ord_PCoA <- function(dataset=ExprSet_species,
                     Method="bray",
                     Group_info="Group"){

  metadata <- pData(dataset)
  colnames(metadata)[which(colnames(metadata) == Group_info)] <- "Group"
  profile <- exprs(dataset)

  # pcoa
  pcoa <- ape::pcoa(vegan::vegdist(t(profile), method = Method))
  eig <- pcoa$values[, "Eigenvalues"]
  eig_var <- eig[1:2]
  eig_var_explain <- round(eig_var/sum(eig), 4) * 100
  # explains variable
  explains <- paste0(paste0("PCoA", seq(2)), " (", paste0(eig_var_explain, "%"), ")")
  score <- dplyr::inner_join(pcoa$vectors[, c(1:2)] %>% data.frame() %>%
                               dplyr::select(c(1:2)) %>%
                               stats::setNames(paste0("Axis", seq(2))) %>%
                               tibble::rownames_to_column("SampleID"),
                             metadata %>% tibble::rownames_to_column("SampleID") %>%
                               dplyr::select(all_of(c("SampleID", "Group"))),
                             by = "SampleID")
  # PERMANOVA
  set.seed(123)
  if(any(profile < 0)){
    res_adonis <- vegan::adonis(vegan::vegdist(t(profile), method = "manhattan") ~ metadata$Group, permutations = 999)
  }else{
    res_adonis <- vegan::adonis(vegan::vegdist(t(profile), method = "bray") ~ metadata$Group, permutations = 999)
  }
  adn_pvalue <- res_adonis[[1]][["Pr(>F)"]][1]
  adn_rsquared <- round(res_adonis[[1]][["R2"]][1],3)
  #use the bquote function to format adonis results to be annotated on the ordination plot.
  signi_label <- paste(cut(adn_pvalue,
                           breaks=c(-Inf, 0.001, 0.01, 0.05, Inf),
                           label=c("***", "**", "*", ".")))
  PERMANOVA <- bquote(atop(atop("PERMANOVA", R^2==~.(adn_rsquared)),
                           atop("p-value="~.(adn_pvalue)~.(signi_label), phantom())))

  res <- list(PCoA=score,
              epn=explains,
              PER=PERMANOVA)

  return(res)
}


#' @importFrom MASS isoMDS sammon
#' @importFrom dplyr %>% select inner_join
#' @importFrom tibble column_to_rownames column_to_rownames
#' @importFrom vegan adonis vegdist
#' @importFrom stats setNames dist cmdscale
#' @import convert
#'
#'
#' @title Dimensionality Reduction Analysis: Multidimensional Scaling(MDS)
#'
#' @description
#' MDS use the distance among samples, which calculates through the multiple variables.
#'
#' @details 12/2/2021 Guangzhou China
#' @author  Hua Zou
#'
#' @param Expression, ExpressionSet; (Required) ExpressionSet object.
#' @param Method, Character; method for MDS (default: Classic).
#' @param Group_info, Character; the group for plot(default: "Group").
#'
#' @return
#' a list object:
#'   MDS score
#'   Result of PERMANOVA
#'
#' @usage Ord_MDS(dataset=ExpressionSet, Method="Classic", Group_info="Group")
#' @examples
#'
#' data(ExprSet_species)
#'
#' MDS_res <- Ord_MDS(dataset=ExprSet_species, Method="Classic", Group_info="Group")
#' MDS_res$MDS
#'
Ord_MDS <- function(dataset=ExprSet_species,
                    Method="Classic",
                    Group_info="Group"){

  metadata <- pData(dataset)
  colnames(metadata)[which(colnames(metadata) == Group_info)] <- "Group"
  profile <- exprs(dataset)

  # MDS
  if(is.element(Method, "Classic")){
    # classical (metric) multidimensional scaling
    mds <- t(profile) %>% stats::dist() %>% stats::cmdscale() %>% tibble::as_tibble()
  }else if(is.element(Method, "Non-metric")){
    # Kruskal’s non-metric multidimensional scaling
    mds <- t(profile) %>% stats::dist() %>% MASS::isoMDS() %>% .$points %>% tibble::as_tibble()
  }else if(is.element(Method, "Non-Sammon")){
    # sammon’s non-linear mapping
    mds <- t(profile) %>% stats::dist() %>% MASS::sammon() %>% .$points %>% tibble::as_tibble()
  }
  rownames(mds) <- colnames(profile)
  score <- dplyr::inner_join(mds %>% stats::setNames(paste0("Axis", seq(2))) %>%
                               tibble::rownames_to_column("SampleID"),
                             metadata %>% tibble::rownames_to_column("SampleID") %>%
                               dplyr::select(all_of(c("SampleID", "Group"))),
                             by = "SampleID")
  # PERMANOVA
  set.seed(123)
  if(any(profile < 0)){
    res_adonis <- vegan::adonis(vegan::vegdist(t(profile), method = "manhattan") ~ metadata$Group, permutations = 999)
  }else{
    res_adonis <- vegan::adonis(vegan::vegdist(t(profile), method = "bray") ~ metadata$Group, permutations = 999)
  }
  adn_pvalue <- res_adonis[[1]][["Pr(>F)"]][1]
  adn_rsquared <- round(res_adonis[[1]][["R2"]][1],3)
  #use the bquote function to format adonis results to be annotated on the ordination plot.
  signi_label <- paste(cut(adn_pvalue,
                           breaks=c(-Inf, 0.001, 0.01, 0.05, Inf),
                           label=c("***", "**", "*", ".")))
  PERMANOVA <- bquote(atop(atop("PERMANOVA", R^2==~.(adn_rsquared)),
                           atop("p-value="~.(adn_pvalue)~.(signi_label), phantom())))

  res <- list(MDS=score,
              epn=c("MDS1", "MDS2"),
              PER=PERMANOVA)

  return(res)
}


#' @importFrom Rtsne Rtsne
#' @importFrom dplyr %>% select inner_join
#' @importFrom tibble column_to_rownames column_to_rownames
#' @importFrom vegan adonis vegdist
#' @importFrom stats setNames
#' @import convert
#'
#'
#' @title Dimensionality Reduction Analysis: t-Distributed Stochastic Neighbor Embedding(t-SNE)
#'
#' @description
#' Visualization of High Dimensional Data using t-SNE with R.
#'
#' @details 12/2/2021 Guangzhou China
#' @author  Hua Zou
#'
#' @param Expression, ExpressionSet; (Required) ExpressionSet object.
#' @param Perplexity, Numeric; (Required) numeric; Perplexity parameter (should not be bigger than 3 * perplexity < nrow(X) - 1, see details for interpretation).
#' @param Group_info, Character; the group for plot(default: "Group").
#'
#' @return
#' a list object:
#'   t-SNE score
#'   Result of PERMANOVA
#'
#' @usage Ord_tsne(dataset=ExpressionSet, Perplexity=20, Group_info="Group")
#' @examples
#'
#' data(ExprSet_species)
#'
#' tsne_res <- Ord_tsne(dataset=ExprSet_species, Perplexity=20, Group_info="Group")
#' tsne_res$tsne
#'
Ord_tsne <- function(dataset=ExprSet_species,
                     Perplexity=20,
                     Group_info="Group"){

  metadata <- pData(dataset)
  colnames(metadata)[which(colnames(metadata) == Group_info)] <- "Group"
  profile <- exprs(dataset)

  # tsne
  Rtsne <- Rtsne::Rtsne(t(profile),
                 dims=2,
                 perplexity=Perplexity,
                 verbose=FALSE,
                 max_iter=500,
                 eta=200)

  point <- Rtsne$Y %>% data.frame() %>%
    dplyr::select(c(1:2)) %>%
    stats::setNames(c("Axis1", "Axis2"))
  rownames(point) <- colnames(profile)
  score <- dplyr::inner_join(point %>% tibble::rownames_to_column("SampleID"),
                             metadata %>% tibble::rownames_to_column("SampleID") %>%
                               dplyr::select(all_of(c("SampleID", "Group"))),
                             by = "SampleID")
  # PERMANOVA
  set.seed(123)
  if(any(profile < 0)){
    res_adonis <- vegan::adonis(vegan::vegdist(t(profile), method = "manhattan") ~ metadata$Group, permutations = 999)
  }else{
    res_adonis <- vegan::adonis(vegan::vegdist(t(profile), method = "bray") ~ metadata$Group, permutations = 999)
  }
  adn_pvalue <- res_adonis[[1]][["Pr(>F)"]][1]
  adn_rsquared <- round(res_adonis[[1]][["R2"]][1],3)
  #use the bquote function to format adonis results to be annotated on the ordination plot.
  signi_label <- paste(cut(adn_pvalue,
                           breaks=c(-Inf, 0.001, 0.01, 0.05, Inf),
                           label=c("***", "**", "*", ".")))
  PERMANOVA <- bquote(atop(atop("PERMANOVA", R^2==~.(adn_rsquared)),
                           atop("p-value="~.(adn_pvalue)~.(signi_label), phantom())))

  res <- list(tsne=score,
              epn=c("tSNE1", "tSNE2"),
              PER=PERMANOVA)

  return(res)
}

#' @importFrom umap umap
#' @importFrom dplyr %>% select inner_join
#' @importFrom tibble column_to_rownames column_to_rownames
#' @importFrom vegan adonis vegdist
#' @importFrom stats setNames
#' @import convert
#'
#'
#' @title Dimensionality Reduction Analysis: Uniform Manifold Approximation and Projection (UMAP)
#'
#' @description
#' UMAP: a non-linear dimensionality reduction algorithm.
#'
#' @details 12/2/2021 Guangzhou China
#' @author  Hua Zou
#'
#' @param Expression, ExpressionSet; (Required) ExpressionSet object.
#' @param Group_info, Character; the group for plot(default: "Group").
#'
#' @return
#' a list object:
#'   umap score
#'   Result of PERMANOVA
#'
#' @usage Ord_umap(dataset=ExpressionSet, Group_info="Group")
#' @examples
#'
#' data(ExprSet_species)
#'
#' umap_res <- Ord_umap(dataset=ExprSet_species, Group_info="Group")
#' umap_res$umap
#'
Ord_umap <- function(dataset=ExprSet_species,
                     Group_info="Group"){

  metadata <- pData(dataset)
  colnames(metadata)[which(colnames(metadata) == Group_info)] <- "Group"
  profile <- exprs(dataset)

  # tsne
  Umap <- umap::umap(t(profile))
  point <- Umap$layout %>% data.frame() %>%
    dplyr::select(c(1:2)) %>%
    stats::setNames(c("Axis1", "Axis2"))
  rownames(point) <- colnames(profile)
  score <- dplyr::inner_join(point %>% tibble::rownames_to_column("SampleID"),
                             metadata %>% tibble::rownames_to_column("SampleID") %>%
                               dplyr::select(all_of(c("SampleID", "Group"))),
                             by = "SampleID")
  # PERMANOVA
  set.seed(123)
  if(any(profile < 0)){
    res_adonis <- vegan::adonis(vegan::vegdist(t(profile), method = "manhattan") ~ metadata$Group, permutations = 999)
  }else{
    res_adonis <- vegan::adonis(vegan::vegdist(t(profile), method = "bray") ~ metadata$Group, permutations = 999)
  }
  adn_pvalue <- res_adonis[[1]][["Pr(>F)"]][1]
  adn_rsquared <- round(res_adonis[[1]][["R2"]][1],3)
  #use the bquote function to format adonis results to be annotated on the ordination plot.
  signi_label <- paste(cut(adn_pvalue,
                           breaks=c(-Inf, 0.001, 0.01, 0.05, Inf),
                           label=c("***", "**", "*", ".")))
  PERMANOVA <- bquote(atop(atop("PERMANOVA", R^2==~.(adn_rsquared)),
                           atop("p-value="~.(adn_pvalue)~.(signi_label), phantom())))

  res <- list(umap=score,
              epn=c("UMAP1", "UMAP2"),
              PER=PERMANOVA)

  return(res)
}


#' @importFrom tibble column_to_rownames column_to_rownames
#' @import ggplot2
#' @importFrom cowplot insert_xaxis_grob insert_yaxis_grob ggdraw
#'
#'
#' @title Plot Ordination result with scatterplot
#'
#' @description
#' Show the result of Ordination using scatterplot.
#'
#' @details 12/2/2021 Guangzhou China
#' @author  Hua Zou
#'
#' @param Score, Data.frame; (Required) Result of Ordination.
#' @param Axis_name, Character; the labels of two axis.
#' @param Pvalue, Character; (Required) PERMANOVA result.
#'
#' @return
#' a scatterplot
#'
#' @usage Ord_Plot(Score=score, Axis_name=explains, Pvalue=label)
#' @examples
#'
#' data(ExprSet_species)
#'
#' PCA_res <- Ord_PCA(dataset=ExprSet_species, Normalization=TRUE, Group_info="Group")
#' Ord_Plot(Score=PCA_res$PCA, Axis_name=PCA_res$epn, Pvalue=PCA_res$PER)
#'
#' res_tsne <- Ord_tsne(dataset=ExprSet_species, Perplexity=20, Group_info="Group")
#' Ord_Plot(Score=res_tsne$tsne, Axis_name=res_tsne$epn, Pvalue=res_tsne$PER)
#'
#'
Ord_Plot <- function(Score=score,
                     Axis_name=name,
                     Pvalue=permanova){

  Score$Group <- factor(Score$Group)
  # main plot
  pmain <- ggplot(data=Score, aes(x=Axis1, y=Axis2))+
    geom_point(aes(fill=Group), size=2, shape=21, stroke=.8, color="black")+
    stat_ellipse(aes(color=Group), level=0.95, linetype=1, size=1.5)+
    annotate("text",
             x=max(Score$Axis1),
             y=min(Score$Axis1),
             label=Pvalue,
             size=6)+
    guides(color="none")+
    theme_bw()+
    theme(axis.title=element_text(size=10, color="black", face="bold"),
          axis.text=element_text(size=9, color="black"),
          text=element_text(size=8, color="black", family="serif"),
          strip.text=element_text(size=9, color="black", face="bold"),
          panel.grid=element_blank(),
          legend.position=c(0, 0),
          legend.justification=c(0, 0),
          legend.title=element_text(size=11, color="black", family="serif"),
          legend.text=element_text(size=10, color="black", family="serif"),
          legend.background=element_rect(color="black", fill="white", size=0.5))


  pmain <- pmain + labs(x=Axis_name[1], y=Axis_name[2])

  # Marginal bolxplot along x axis
  xbp <- ggplot(data=Score, aes(x=Group, y=Axis1, fill=Group))+
    stat_boxplot(geom="errorbar", width=.12)+
    geom_boxplot(width=.2, outlier.shape=3, outlier.size=1)+
    labs(x="", y="")+
    guides(fill="none")+
    theme_void()+
    theme(axis.text=element_blank(),
          panel.grid=element_blank(),
          axis.ticks=element_blank())

  # Marginal bolxplot along y axis
  ybp <- ggplot(data=Score, aes(x=Group, y=Axis2, fill=Group))+
    stat_boxplot(geom="errorbar", width=.12)+
    geom_boxplot(width=.2, outlier.shape=3, outlier.size=1)+
    labs(x="", y="")+
    guides(fill="none")+
    coord_flip()+
    theme_void()+
    theme(axis.text=element_blank(),
          panel.grid=element_blank(),
          axis.ticks=element_blank())

  # merge three plot
  p1 <- cowplot::insert_xaxis_grob(pmain, ybp, grid::unit(.2, "null"), position = "top")
  p2 <- cowplot::insert_yaxis_grob(p1, xbp, grid::unit(.2, "null"), position = "right")
  pl <- cowplot::ggdraw(p2)
  return(pl)
}
