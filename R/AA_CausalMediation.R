#' @title Causal Mediation Analysis (CMA)
#'
#' @description
#'  CMA examines the intermediate process(Mediation variable:M) by which the independent variable(X) affects the dependent variable(Y).
#'  Firstly, building Mediator Model(model.m): f(M | X + AdjVar); and then Outcome Model: f(Y = X + M + AdjVar).
#'  Secondly, mediation_out <- mediate{Y, M, X, AdjVar, model.y, model.m} to perform mediation analysis.
#'  Finally, extracting the results of mediation_out.
#'
#'
#' @details 12/5/2021 Guangzhou China
#' @author  Hua Zou
#'
#' @param dataset, Matrix; (Required) Metadata.
#' @param XVar, Character; independent variable.
#' @param YVar, Character; dependent variable.
#' @param MVar, Character; mediator.
#' @param AdjVar, Character; adjust variable(default: AdjVar=NULL).
#' @param Package, Character; (Required) package for CMA (default: Package="mediation").
#'
#' @return
#'  a list of results
#'   CMA model
#'   Sensitivity analysis of CMA model
#'   Table of CMA model
#'
#' @export
#'
#' @importFrom dplyr %>% select all_of
#' @importFrom stats formula na.omit glm
#' @importFrom mediation mediate
#' @importFrom intmed mediate
#'
#' @usage AA_Mediation(dataset=Metadata, XVar="Lysine", YVar="Tryptophan", AdjVar=c("Age", "Gender"), Package="mediation")
#' @examples
#'
#' \donttest{
#' data(ExprSet_species)
#' library(Biobase)
#' Metadata <- pData(ExprSet_species)
#'
#' # mediation package & one mediation
#' Mediation_res1 <- AA_Mediation(dataset=Metadata, XVar="Lysine", YVar="Tryptophan", AdjVar=c("Age", "Gender"), MVar="BMI", Package="mediation")
#'
#' # intmed package or two mediation
#' Mediation_res2 <- AA_Mediation(dataset=Metadata, XVar="Lysine", YVar="Tryptophan", AdjVar=c("Age", "Gender"), MVar=c("BMI", "Age"), Package="intmed")
#'
#'}
#'
AA_Mediation <- function(dataset=Metadata,
                         XVar="Lysine",
                         YVar="Tryptophan",
                         MVar="BMI",
                         AdjVar=c("Age", "Gender"),
                         Package="mediation"){

  dat <- dataset %>% dplyr::select(dplyr::all_of(c(XVar, YVar, MVar, AdjVar))) %>%
    stats::na.omit()
  if(nrow(dat) == 0){
    stop("Data with too many missing values couldn't be used please check your data")
  }


  # mediation package
  if(all(is.element(Package, "mediation"), length(MVar) == 1)){

    # formula
    if(is.null(AdjVar)){
      m_formula <- stats::formula(paste0(MVar, "~", XVar))
      y_formula <- stats::formula(paste0(YVar, "~", XVar, "+" , MVar))
    }else{
      m_formula <- stats::formula(paste0(MVar, "~", XVar, "+", paste(AdjVar, collapse = "+")))
      y_formula <- stats::formula(paste0(YVar, "~", XVar, "+" , MVar, "+", paste(AdjVar, collapse = "+")))
    }

    # extract the results of mediation::mediate(https://stackoverflow.com/questions/41582486/how-to-convert-r-mediation-summary-to-data-frame)
    # ACM: Average Causal Mediated Effect
    # Average Direct Effect
    extract_mediation_summary <- function(x){

      clp <- 100 * x$conf.level
      isLinear.y <- ((class(x$model.y)[1] %in% c("lm", "rq")) ||
                       (inherits(x$model.y, "glm") && x$model.y$family$family ==
                          "gaussian" && x$model.y$family$link == "identity") ||
                       (inherits(x$model.y, "survreg") && x$model.y$dist ==
                          "gaussian"))

      printone <- !x$INT && isLinear.y

      if(printone){

        smat <- c(x$d1, x$d1.ci, x$d1.p)
        smat <- rbind(smat, c(x$z0, x$z0.ci, x$z0.p))
        smat <- rbind(smat, c(x$tau.coef, x$tau.ci, x$tau.p))
        smat <- rbind(smat, c(x$n0, x$n0.ci, x$n0.p))

        rownames(smat) <- c("ACME", "ADE", "Total Effect", "Prop. Mediated")

      }else{

        smat <- c(x$d0, x$d0.ci, x$d0.p)
        smat <- rbind(smat, c(x$d1, x$d1.ci, x$d1.p))
        smat <- rbind(smat, c(x$z0, x$z0.ci, x$z0.p))
        smat <- rbind(smat, c(x$z1, x$z1.ci, x$z1.p))
        smat <- rbind(smat, c(x$tau.coef, x$tau.ci, x$tau.p))
        smat <- rbind(smat, c(x$n0, x$n0.ci, x$n0.p))
        smat <- rbind(smat, c(x$n1, x$n1.ci, x$n1.p))
        smat <- rbind(smat, c(x$d.avg, x$d.avg.ci, x$d.avg.p))
        smat <- rbind(smat, c(x$z.avg, x$z.avg.ci, x$z.avg.p))
        smat <- rbind(smat, c(x$n.avg, x$n.avg.ci, x$n.avg.p))

        rownames(smat) <- c("ACME (control)", "ACME (treated)",
                            "ADE (control)", "ADE (treated)", "Total Effect",
                            "Prop. Mediated (control)", "Prop. Mediated (treated)",
                            "ACME (average)", "ADE (average)", "Prop. Mediated (average)")

      }

      colnames(smat) <- c("Estimate", paste(clp, "% CI Lower", sep = ""),
                          paste(clp, "% CI Upper", sep = ""), "p-value")
      return(smat)

    }

    # fitting model use glm for either Continuous or Binary variables
    model_m <- stats::glm(m_formula, data = dat)
    model_y <- stats::glm(y_formula, data = dat)
    model_med <- mediation::mediate(model_m,
                                    model_y,
                                    treat = XVar,
                                    mediator = MVar,
                                    sims = 1000)
    if(class(model_m)[1] == "lm"){
      model_med_sens <- mediation::medsens(model_med, rho.by = 0.05)
    }else{
      model_med_sens <- NULL
    }

    model_med_df <- extract_mediation_summary(model_med)
  }else if(any(is.element(Package, "intmed"), length(MVar) > 1)){

    ModelType <- function(x, d=dat){
      if(is.numeric(d[, x])){
        model_type <- "regression"
      }else if(integer(d[, x])){
        model_type <- "poisson regression"
      }else{
        model_type <- "logistic regression"
      }

      return(model_type)
    }
    Ymodel <- as.character(sapply(YVar, ModelType))
    Mmodel <- as.character(sapply(MVar, ModelType))

    model_med <- intmed::mediate(y = YVar,
                                 med = MVar,
                                 treat = XVar,
                                 c = AdjVar,
                                 ymodel = Ymodel,
                                 mmodel = Mmodel,
                                 treat_lv = 1,
                                 control_lv = 0,
                                 incint = FALSE,
                                 inc_mmint = FALSE,
                                 conf.level = 0.95,
                                 data = dat,
                                 sim = 1000,
                                 digits = 3)
    model_med_sens <- NULL
    model_med_df <- NULL
  }

  res <- list(fit=model_med,
              sens=model_med_sens,
              df=model_med_df)

  return(res)
}
