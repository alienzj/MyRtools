#' Independent T-Test
#'
#' @description
#' Tests for the difference between the same variable from different populations (e.g., comparing boys to girls)
#'
#' @details 07/18/2019
#' @author  Hua Zou
#'
#' @param phen phenotype with sampleID and group; sampleID connected to profile.
#' @param prof profile table rownames->taxonomy; colnames->sampleID
#' @param DNAID names of sampleID to connect phen and prof
#' @param GROUP names of group information, only contail two levels if grp1 or grp2 haven't been provided
#' @param grp1  one of groups to be converted into 0
#' @param grp2  one of groups to be converted into 1
#'
#' @usage result <- unpaired_ttest(phen, prof, "SampleID", "Stage", "BASE", "WASH")
#' @return Returns a result of Independent T-Test
#'   type:       kind of data
#'   Block:      group information
#'   Num:        number of group
#'   P-value:    P by Independent T-Test
#'   FDR:        adjusted by BH
#'   Enrichment: directory by median
#'               directory by mean
#'   Occurence:  occurence of two groups
#'   median:     both or each group
#'   mean:       each group
#'   FDR:        adjusted P value by BH
#'   Odds Ratio:     95% Confidence interval
#'
#' @export
#'
unpaired_ttest <- function(x, y, DNAID, GROUP,
                           grp1=NULL, grp2=NULL){

  # determine x with two cols and names are corret
  phe <- x %>% select(DNAID, GROUP)
  colnames(phe)[which(colnames(phe) == DNAID)] <- "SampleID"
  colnames(phe)[which(colnames(phe) == GROUP)] <- "Stage"
  if (length(which(colnames(phe)%in%c("SampleID","Stage"))) != 2){
    warning("x without 2 cols: DNAID, GROUP")
  }

  # select groups
  if(length(grp1)){
    phe.cln <- phe %>% filter(Stage%in%c(grp1, grp2)) %>%
      mutate(Stage=factor(Stage, levels = c(grp1, grp2)))
    pr <- c(grp1, grp2)
  } else{
    phe.cln <- phe %>% mutate(Stage=factor(Stage))
    pr <- levels(phe.cln$Stage)
  }

  if (length(levels(phe.cln$Stage)) > 2) {
    stop("The levels of `group` are more than 2")
  }

  # profile
  sid <- intersect(phe.cln$SampleID, colnames(y))
  prf <- y %>% select(sid) %>%
    rownames_to_column("tmp") %>%
    # occurrence of rows more than 0.1
    filter(apply(select(., -one_of("tmp")), 1, function(x){sum(x > 0)/length(x)}) > 0.3) %>%
    data.frame() %>% column_to_rownames("tmp") %>%
    t() %>% data.frame()

  # judge no row of profile filter
  if (ncol(prf) == 0) {
    stop("No row of profile to be choosed\n")
  }

  # merge phenotype and profile
  mdat <- inner_join(phe.cln %>% filter(SampleID%in%sid),
                     prf %>% rownames_to_column("SampleID"),
                     by = "SampleID")
  dat.phe <- mdat %>% select(c(1:2))
  dat.prf <- mdat %>% select(-2)

  res <- apply(dat.prf[, -1], 2, function(x, grp){
    dat <- as.numeric(x)
    p <- signif(t.test(dat ~ grp, paired=F)$p.value, 6)
    # median
    md <- signif(median(dat), 4)
    mdn <- signif(tapply(dat, grp, median), 4)
    if ( mdn[1] > mdn[2] & p < 0.05) {
      enrich1 <- pr[1]
    } else if (mdn[1] < mdn[2] & p < 0.05) {
      enrich1 <- pr[2]
    } else if (p > 0.05 | mdn[1] == mdn[2]){
      enrich1 <- "No significance"
    }

    # mean
    mn <- signif(tapply(dat, grp, mean), 4)
    if ( mn[1] > mn[2] & p < 0.05) {
      enrich2 <- pr[1]
    } else if (mn[1] < mn[2] & p < 0.05) {
      enrich2 <- pr[2]
    } else if (p > 0.05 | mn[1] == mn[2]){
      enrich2 <- "No significance"
    }
    occ <- signif(tapply(dat, grp, function(x){
      round(sum(x > 0)/length(x), 4)}), 4)

    res <- c(p,enrich1,enrich2,occ,md,mdn,mn)
    return(res)
  }, dat.phe$Stage) %>%
    t(.) %>% data.frame(.) %>%
    rownames_to_column("type") %>%
    varhandle::unfactor(.)

  colnames(res)[2:11] <- c("Pvalue", "Enrich_median", "Enrich_mean",
                           paste0(pr, "_occurence"), "median_all",
                           paste0(pr, "_median"), paste0(pr, "_mn"))
  res$Block <- paste0(pr[1], "_vs_", pr[2])
  number <- as.numeric(table(dat.phe$Stage))
  res$Num <- paste0(pr[1], number[1], "_vs_",
                    pr[2], number[2])
  res.cln <- res %>% select(c(1,12:13, 2:11)) %>%
    mutate(Pvalue=as.numeric(Pvalue)) %>%
    mutate(FDR=p.adjust(Pvalue, method = "BH")) %>%
    arrange(FDR, Pvalue)
  res2 <- res.cln[,c(1:4,14,5:13)]


  # scale profile
  dat.prf.cln <- prf[, -1]
  dat.phe.cln <- dat.phe %>% mutate(Group=ifelse(Stage==pr[1], 0, 1))
  idx <- which(colnames(dat.phe.cln) == "Group")

  # glm result for odd ratios 95%CI
  glmFun <- function(m, n){
    # calculate the glm between profile and group information
    #
    # Args:
    #   m:  result of group information which must to be numeric
    #   n:  taxonomy to be glm
    #
    # Returns:
    #   the glm result of between taxonomy group
    group <- m
    marker <- n
    model <- summary(glm(marker ~ group, family = binomial(link = "logit")))
    res <- signif(exp(model$coefficients["group",1]) +
                    qnorm(c(0.025,0.5,0.975)) * model$coefficients["group",1], 2)

    return(res)
  }

  glm_res <- t(apply(dat.prf.cln, 2, function(x, group){
    res <- glmFun(group, as.numeric(x))
    return(res)
  }, group = dat.phe.cln[, idx]))
  Odd <- glm_res %>% data.frame() %>%
    setNames(c("upper", "expected","lower")) %>%
    mutate("Odds Ratio (95% CI)" = paste0(expected, " (", lower, ";", upper, ")"))
  Odd$type <- rownames(glm_res)

  res_merge <- inner_join(res2,
                          Odd[, c(4:5)], by = "type")

  return(res_merge)
}


#' Paired T-Test
#'
#' @description
#' Tests for the difference between two variables from the same population (e.g., a pre- and posttest score)
#'
#' @details 07/18/2019
#' @author  Hua Zou
#'
#' @details 07/18/2019
#' @author  Hua Zou
#'
#' @param phen phenotype with sampleID, ID and group; sampleID connected to profile.
#' @param prof profile table rownames->taxonomy; colnames->sampleID
#' @param DNAID names of sampleID to connect phen and prof
#' @param PID   id for paired test
#' @param GROUP names of group information, only contail two levels if grp1 or grp2 haven't been provided
#' @param grp1  one of groups to be converted into 0
#' @param grp2  one of groups to be converted into 1
#'
#' @usage result <- paired_ttest(phen, prof, "SampleID", "ID","Stage", "BASE", "WASH")
#' @return Returns a result of Paired T-Test
#'   type:       kind of data
#'   Block:      group information
#'   Num:        number of group
#'   P-value:    P by Paired T-Test
#'   FDR:        adjusted by BH
#'   Enrichment: directory by median
#'               directory by mean
#'   Occurence:  occurence of two groups
#'   median:     both or each group
#'   mean:       each group
#'   FDR:        adjusted P value by BH
#'   Odds Ratio:     95% Confidence interval
#'
#' @export
#'
paired_ttest <- function(x,
                         y,
                         DNAID,
                         PID,
                         GROUP,
                         grp1=NULL,
                         grp2=NULL){

  # determine x with two cols and names are corret
  phe <- x %>% select(DNAID, PID, GROUP)
  colnames(phe)[which(colnames(phe) == DNAID)] <- "SampleID"
  colnames(phe)[which(colnames(phe) == PID)] <- "ID"
  colnames(phe)[which(colnames(phe) == GROUP)] <- "Stage"
  if (length(which(colnames(phe)%in%c("SampleID","ID","Stage"))) != 3){
    warning("x without 2 cols: DNAID, ID, GROUP")
  }

  # select groups
  if(length(grp1)){
    phe.cln <- phe %>% filter(Stage%in%c(grp1, grp2)) %>%
      mutate(Stage=factor(Stage, levels = c(grp1, grp2))) %>%
      arrange(ID, Stage)
    pr <- c(grp1, grp2)
  } else {
    phe.cln <- phe %>% mutate(Stage=factor(Stage)) %>%
      arrange(ID, Stage)
    pr <- levels(phe.cln$Stage)
  }

  if (length(levels(phe.cln$Stage)) > 2) {
    stop("The levels of `group` are more than 2")
  }

  # profile
  sid <- intersect(phe.cln$SampleID, colnames(y))
  prf <- y %>% select(sid) %>%
    rownames_to_column("tmp") %>%
    # occurrence of rows more than 0.1
    filter(apply(select(., -one_of("tmp")), 1, function(x){sum(x > 0)/length(x)}) > 0.3) %>%
    data.frame() %>% column_to_rownames("tmp") %>%
    t() %>% data.frame()

  # judge no row of profile filter
  if (ncol(prf) == 0) {
    stop("No row of profile to be choosed\n")
  }

  # determine the right order and group levels
  for(i in 1:nrow(prf)){
    if ((rownames(prf) != phe.cln$SampleID)[i]) {
      stop(paste0(i, " Wrong"))
    }
  }

  # merge phenotype and profile
  mdat <- inner_join(phe.cln %>% filter(SampleID%in%sid),
                     prf %>% rownames_to_column("SampleID"),
                     by = "SampleID")
  dat.phe <- mdat %>% select(c(1:3))
  dat.prf <- mdat %>% select(-c(2:3))

  res <- apply(dat.prf[, -1], 2, function(x, grp){
    dat <- as.numeric(x)
    p <- signif(t.test(dat ~ grp, paired=T)$p.value, 6)
    # median
    md <- signif(median(dat), 4)
    mdn <- signif(tapply(dat, grp, median), 4)
    if ( mdn[1] > mdn[2] & p < 0.05) {
      enrich1 <- pr[1]
    } else if (mdn[1] < mdn[2] & p < 0.05) {
      enrich1 <- pr[2]
    } else if (p > 0.05 | mdn[1] == mdn[2]){
      enrich1 <- "No significance"
    }

    # mean
    mn <- signif(tapply(dat, grp, mean), 4)
    if ( mn[1] > mn[2] & p < 0.05) {
      enrich2 <- pr[1]
    } else if (mn[1] < mn[2] & p < 0.05) {
      enrich2 <- pr[2]
    } else if (p > 0.05 | mn[1] == mn[2]){
      enrich2 <- "No significance"
    }
    occ <- signif(tapply(dat, grp, function(x){
      round(sum(x > 0)/length(x), 4)}), 4)

    res <- c(p,enrich1,enrich2,occ,md,mdn,mn)

    return(res)
  }, dat.phe$Stage) %>%
    t(.) %>% data.frame(.) %>%
    rownames_to_column("type") %>%
    varhandle::unfactor(.)

  colnames(res)[2:11] <- c("Pvalue", "Enrich_median", "Enrich_mean",
                           paste0(pr, "_occurence"), "median_all",
                           paste0(pr, "_median"), paste0(pr, "_mean"))
  res$Block <- paste0(pr[1], "_vs_", pr[2])
  number <- as.numeric(table(dat.phe$Stage))
  res$Num <- paste0(pr[1], number[1], "_vs_",
                    pr[2], number[2])
  res.cln <- res %>% select(c(1,12:13, 2:11)) %>%
    mutate(Pvalue=as.numeric(Pvalue)) %>%
    mutate(FDR=p.adjust(Pvalue, method = "BH")) %>%
    arrange(FDR, Pvalue)
  res2 <- res.cln[,c(1:4,14,5:13)]


  # scale profile
  dat.prf.cln <- prf[, -1]
  dat.phe.cln <- dat.phe %>% mutate(Group=ifelse(Stage==pr[1], 0, 1))
  idx <- which(colnames(dat.phe.cln) == "Group")

  # glm result for odd ratios 95%CI
  glmFun <- function(m, n){
    # calculate the glm between profile and group information
    #
    # Args:
    #   m:  result of group information which must to be numeric
    #   n:  taxonomy to be glm
    #
    # Returns:
    #   the glm result of between taxonomy group
    group <- m
    marker <- n
    model <- summary(glm(marker ~ group, family = binomial(link = "logit")))
    res <- signif(exp(model$coefficients["group",1]) +
                    qnorm(c(0.025,0.5,0.975)) * model$coefficients["group",1], 2)

    return(res)
  }

  glm_res <- t(apply(dat.prf.cln, 2, function(x, group){
    res <- glmFun(group, as.numeric(x))
    return(res)
  }, group = dat.phe.cln[, idx]))
  Odd <- glm_res %>% data.frame() %>%
    setNames(c("upper", "expected","lower")) %>%
    mutate("Odds Ratio (95% CI)" = paste0(expected, " (", lower, ";", upper, ")"))
  Odd$type <- rownames(glm_res)

  res_merge <- inner_join(res2,
                          Odd[, c(4:5)], by = "type")

  return(res_merge)
}


#' The one-way analysis of variance (ANOVA)
#'
#' @description
#' one-factor ANOVA, is an extension of independent two-samples t-test for comparing means in a situation
#' where there are more than two groups. In one-way ANOVA, the data is organized into several groups base on
#' one single grouping variable
#'
#' @details 07/18/2019
#' @author  Hua Zou
#'
#' @param phen phenotype with sampleID and group; sampleID connected to profile.
#' @param prof profile table rownames->taxonomy; colnames->sampleID
#' @param DNAID names of sampleID to connect phen and prof
#' @param GROUP names of group information
#'
#' @usage result <- ANOVA_one(phen, prof, "SampleID", "Stage")
#' @return Returns a result of the one-way analysis of variance (ANOVA)
#'   type:       kind of data
#'   Num:        number of group
#'   P-value:    P by one-way analysis of variance (ANOVA)
#'               post test pvalue
#'   FDR:        adjusted by BH
#'   mean+/-sd   each group
#'   median:     each group
#'
#' @export
#'
ANOVA_one <- function(x, y, DNAID, GROUP,
                       grp1=NULL, grp2=NULL,grp3=NULL){

  # determine x with two cols and names are corret
  phe <- x %>% select(DNAID, GROUP)
  colnames(phe)[which(colnames(phe) == DNAID)] <- "SampleID"
  colnames(phe)[which(colnames(phe) == GROUP)] <- "Stage"
  if (length(which(colnames(phe)%in%c("SampleID","Stage"))) != 2){
    warning("x without 2 cols: DNAID, GROUP")
  }

  # select groups
  if(length(grp1)){
    phe.cln <- phe %>% filter(Stage%in%c(grp1, grp2, grp3)) %>%
      mutate(Stage=factor(Stage, levels = c(grp1, grp2, grp3)))
    pr <- c(grp1, grp2, grp3)
  } else {
    phe.cln <- phe %>% mutate(Stage=factor(Stage))
    pr <- levels(phe.cln$Stage)
  }

  if (length(levels(phe.cln$Stage)) < 2) {
    stop("The levels of `GROUP` no more than 2")
  }

  # profile
  sid <- intersect(phe.cln$SampleID, colnames(y))
  prf <- y %>% select(sid) %>%
    rownames_to_column("tmp") %>%
    # occurrence of rows more than 0.1
    filter(apply(select(., -one_of("tmp")), 1, function(x){sum(x > 0)/length(x)}) > 0.3) %>%
    data.frame() %>% column_to_rownames("tmp") %>%
    t() %>% data.frame()

  # judge no row of profile filter
  if (ncol(prf) == 0) {
    stop("No row of profile to be choosed\n")
  }

  # determine the right order and group levels
  for(i in 1:nrow(prf)){
    if ((rownames(prf) != phe.cln$SampleID)[i]) {
      stop(paste0(i, " Wrong"))
    }
  }

  # merge phenotype and profile
  mdat <- inner_join(phe.cln %>% filter(SampleID%in%sid),
                     prf %>% rownames_to_column("SampleID"),
                     by = "SampleID")
  dat.phe <- mdat %>% select(c(1:2))
  dat.prf <- mdat %>% select(-c(2))

  res <- apply(dat.prf[, -1], 2, function(x, grp){
    dat <- data.frame(y=as.numeric(x), grp)

    # anova; mean±sd
    res.aov <- aov(y ~ Stage, data = dat)
    res.aov.sum <- summary(res.aov)
    pvalue <- signif(res.aov.sum[[1]]$`Pr(>F)`[1], 4)

    # Tukey multiple pairwise-comparisons
    mtest <- TukeyHSD(res.aov)
    p.post <- mtest$Stage[, 4]

    mn <- tapply(dat$y, dat$Stage, function(x){
      num <- paste(signif(mean(x[!is.na(x)]), 4),
                   signif(sd(x[!is.na(x)]), 4),
                   sep= "+/-")})
    md <- tapply(dat$y, dat$Stage, function(x){
      signif(median(x[!is.na(x)]), 4)})
    nu <- tapply(dat$y, dat$Stage,function(x){
      length(x)})
    Num <-  paste(paste(pr, nu, sep="_"), collapse = " vs ")

    res <- c(pvalue, mn, md, Num, p.post)

    return(res)
  }, dat.phe[, 2, F]) %>% t(.) %>% data.frame(.) %>%
    rownames_to_column("tmp") %>% varhandle::unfactor(.)

  mean_sd <- unlist(lapply(pr, function(x){paste0(x, "\nMean+/-Sd")}))
  MD <- unlist(lapply(pr, function(x){paste0(x, "_Median")}))
  post_name <- NULL
  for(i in 1:(length(pr)-1)){
    for(j in (i+1):length(pr)){
      post_name <- c(post_name, paste("P.value\n", pr[i], "_vs_", pr[j], sep = ""))
    }
  }
  colnames(res)[2:ncol(res)] <- c("P.value", mean_sd, MD, "Number", post_name)
  res$FDR <- with(res, p.adjust(as.numeric(P.value), method = "BH"))

  return(res[, c(1,9,2,13,10:12,3:8)])
}


#' Two-Way ANOVA Test
#'
#' @description
#' Two-way ANOVA test is used to evaluate simultaneously the effect of
#' two grouping variables (A and B) on a response variable
#'
#' @details 07/18/2019
#' @author  Hua Zou
#'
#' @param phen phenotype with sampleID and group; sampleID connected to profile.
#' @param prof profile table rownames->taxonomy; colnames->sampleID
#' @param DNAID names of sampleID to connect phen and prof
#' @param GROUP1 names of group information: factor1
#' @param GROUP2 names of group information: factor2
#'
#' @usage result <- ANOVA_one(phen, prof, "SampleID", "Stage", "Group")
#' @return Returns a result of Two-Way ANOVA Test
#'   type:       kind of data
#'   Num:        number of group
#'   P-value:    P by Two-Way ANOVA Test
#'               post test pvalue
#'   FDR:        adjusted by BH
#'   mean+/-sd   each group
#'   median:     each group
#'
#' @export
#'
ANOVA_two <- function(x, y, DNAID, GROUP1, GROUP2){

  # determine x with two cols and names are corret
  phe <- x %>% select(DNAID, GROUP1, GROUP2)

  colnames(phe)[which(colnames(phe) == DNAID)] <- "SampleID"
  colnames(phe)[which(colnames(phe) == GROUP1)] <- "Group1"
  colnames(phe)[which(colnames(phe) == GROUP2)] <- "Group2"
  if (length(which(colnames(phe)%in%c("SampleID","Group1", "Group2"))) != 3){
    warning("x without 3 cols: DNAID, GROUP1, GROUP2")
  }

  # profile
  sid <- intersect(phe$SampleID, colnames(y))
  prf <- y %>% select(sid) %>%
    rownames_to_column("tmp") %>%
    # occurrence of rows more than 0.1
    filter(apply(select(., -one_of("tmp")), 1, function(x){sum(x > 0)/length(x)}) > 0.3) %>%
    data.frame() %>% column_to_rownames("tmp") %>%
    t() %>% data.frame()
  phe.cln <- phe %>% filter(SampleID%in%sid) %>%
    mutate(Group1=factor(Group1),
           Group2=factor(Group2))

  # judge no row of profile filter
  if (ncol(prf) == 0) {
    stop("No row of profile to be choosed\n")
  }

  # determine the right order and group levels
  for(i in 1:nrow(prf)){
    if ((rownames(prf) != phe$SampleID)[i]) {
      stop(paste0(i, " Wrong"))
    }
  }

  # merge phenotype and profile
  mdat <- inner_join(phe.cln, prf %>% rownames_to_column("SampleID"),
                     by = "SampleID")
  dat.phe <- mdat %>% select(c(1:3))
  dat.prf <- mdat %>% select(-c(2:3))

  res <- apply(dat.prf[, -1], 2, function(x, grp){

    #dat <- data.frame(y=as.numeric(dat.prf$Bacteroides), dat.phe[, c(2:3)])

    dat <- data.frame(y=as.numeric(x), grp)

    # anova; mean±sd
    res.aov <- aov(y ~ Group1 + Group2, data = dat)
    res.aov.sum <- summary(res.aov)
    pvalue <- signif(res.aov.sum[[1]]$`Pr(>F)`[1:2], 4)

    # Compute some summary statistics
    stat <- dat %>% group_by(Group1, Group2) %>%
      summarise(
        count = n(),
        mean = signif(mean(y, na.rm = TRUE), 3),
        sd = signif(sd(y, na.rm = TRUE), 3),
        ms = paste(mean, sd, sep="+/-")
      )
    count <- stat$count
    msd <- stat$ms

    # Tukey multiple pairwise-comparisons
    mtest <- TukeyHSD(res.aov)
    p.grp1 <- signif(mtest$Group1[, 4], 3)
    p.grp2 <- signif(mtest$Group2[, 4], 3)

    res <- c(pvalue, p.grp1, p.grp2, count, msd)

    return(res)
  }, dat.phe[, c(2:3)]) %>% t(.) %>% data.frame(.) %>%
    rownames_to_column("tmp") %>% varhandle::unfactor(.)

  # colnames
  pr1 <- levels(dat.phe$Group1)
  pr2 <- levels(dat.phe$Group2)
  ## aov p value
  aov.name <- paste0("P.value_", c(GROUP1, GROUP2))
  ## group1 p value
  grp1.p.name <- NULL
  for(i in 1:(length(pr1)-1)){
    for(j in (i+1):length(pr1)){
      grp1.p.name <- c(grp1.p.name, paste("P.value\n", pr1[i], "_vs_", pr1[j], sep = ""))
    }
  }
  ## group2 p value
  grp2.p.name <- NULL
  for(i in 1:(length(pr2)-1)){
    for(j in (i+1):length(pr2)){
      grp2.p.name <- c(grp2.p.name, paste("P.value_post\n", pr2[i], "_vs_", pr2[j], sep = ""))
    }
  }
  ## count mean sd
  count.name <- NULL
  ms.name <- NULL
  for(i in 1:length(pr1)){
    for(j in 1:length(pr2)){
      count.name <- c(count.name, paste(pr1[i], "_", pr2[j], "\ncount", sep = ""))
      ms.name <- c(ms.name, paste(pr1[i], "_", pr2[j], "\nMean+/-sd", sep = ""))
    }
  }

  colnames(res)[2:ncol(res)] <- c(aov.name, grp1.p.name, grp2.p.name,
                                  count.name, ms.name)
  return(res)
}
