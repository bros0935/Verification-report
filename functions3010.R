if (!require("pacman")) install.packages("pacman")
pacman::p_load(haven, psych, lm.beta, lmhelprs, olsrr, afex, car, emmeans)

descript_hist <- function(x, skew=FALSE, freq=TRUE, mar=c(2,2,1,1), dlty=0, global=FALSE, ...) {
  out1 <- psych::describe(x, skew=skew)
  out2 <- psych::multi.hist(Filter(is.numeric,x), freq=freq, mar=mar, dlty=dlty, global=global)
  return(out1)
  return(out2)
}


cor_lm <- function(fit) {
  dv <- names(fit$model)[1]
  dv_data <- fit$model[, dv]
  ivs <- names(fit$model)[-1]
  iv_data <- fit$model[, ivs, drop=FALSE]
  x <- fit$model
  x_omit <- lapply(ivs, function(X) x[, c(dv, setdiff(ivs, X)), drop = FALSE])
  names(x_omit) <- ivs
  lapply(x_omit, head)
  fits_omit <- lapply(x_omit, function(X) lm(as.formula(paste(dv, "~ .")),data = X))
  resid_omit <- sapply(fits_omit, resid)
  iv_omit <- lapply(ivs, function(X) lm(as.formula(paste(X, "~ .")), data = iv_data))
  resid_iv_omit <- sapply(iv_omit, resid)
  results <- sapply(seq(ivs), function(i) c(zeroorder = cor(iv_data[, i], dv_data), 
                                            partial = cor(resid_iv_omit[, i], resid_omit[, i]), 
                                            part = cor(resid_iv_omit[, i], dv_data)))    
  results <- data.frame(results)
  names(results) <- ivs
  results <- data.frame(t(results))
  results
}

col_diag <- function(x) {
  if(length(names(x$model)[-1]) < 2) {
    data.frame(Variables=names(x$model)[-1],Tolerance=1,VIF=1)
  } else {
    olsrr::ols_vif_tol(x) 
  } 
}

regression_output <- function(..., ci95=FALSE, part=FALSE, coll=FALSE, dp=.5, cf2=FALSE) {
  lm_outs <- list(...)  
  if (length(lm_outs) > 1) {
    lm_sorted <- do.call(lmhelprs::hierarchical, lm_outs)
    if (isTRUE(is.na(lm_sorted))) {
      stop("The models are not hierarchical.")
    }    
  } else {
    lm_sorted <- lm_outs
  }
  data.final <- model.frame(unlist(tail(lm_sorted,n=1)))
  m.form <- lapply(lm_sorted, function(x) deparse(formula(x)))
  if (length(lm_outs)==1) {
    R.sq <- summary(lm_outs[[1]])$r.squared
    bigR <- sqrt(R.sq)
    adj.R.sq <- summary(lm_outs[[1]])$adj.r.squared
    R.sq.change <- R.sq
    fstats <- summary(lm_outs[[1]])$fstatistic
    F.change <- fstats[1]
    p.change <- pf(fstats[1], fstats[2], fstats[3], lower.tail=FALSE)
    df1 <- fstats[2]
    cohens.f2.change <- R.sq.change / (1 - R.sq)
  } else {
    hier_output <- lmhelprs::hierarchical_lm(...)
    R.sq <- hier_output$R.sq
    bigR <- sqrt(R.sq)
    adj.R.sq <- hier_output$adj.R.sq
    R.sq.change <- hier_output$R.sq.change
    F.change <- hier_output$F
    p.change <- hier_output$"Pr(>F)"
    df1 <- hier_output$Df
    cohens.f2.change <- R.sq.change / (1 - R.sq)
  }
  anova_all <- lapply(lm_sorted, function(x) anova(update(x, . ~ 1), x))
  anova_clean <- lapply(anova_all, function(x) {
    ssp <- x[2,"Sum of Sq"]
    ssr <- x[2,"RSS"]
    sst <- ssp + ssr
    dfp <- x[2,"Df"]
    dfr <- x[2,"Res.Df"]
    dft <- dfp + dfr
    msp <- ssp/dfp
    msr <- ssr/dfr
    f <- x[2,"F"]
    p <- x[2,"Pr(>F)"]
    cohens.f2 <- (ssp/(ssp+ssr)) / (1 - (ssp/(ssp+ssr)))
    out <- cbind(c(ssp,ssr,sst), c(dfp,dfr,dft), c(msp,msr,NA), c(f,NA,NA), c(p,NA,NA), c(cohens.f2,NA,NA))
    colnames(out) <- c("SS","df","MS","F","Sig.","cohens.f2")
    rownames(out) <- c("Regression","Residual","Total")
    return(out)
  })
  see <- sqrt(sapply(anova_all, function(x) x[2,"RSS"]/x[2,"Res.Df"]))
  df2 <- sapply(anova_all, function(x) x[2,"Res.Df"])
  R.sq.change[1] <- R.sq[1]
  F.change[1] <- anova_all[[1]]$F[2]
  p.change[1] <- anova_all[[1]]$"Pr(>F)"[2]
  df1[1] <- anova_all[[1]]$Df[2]
  cohens.f2.change[1] <- R.sq.change[1] / (1 - R.sq[1])
  if (length(lm_outs)==1) {
    model.summary <- data.frame(R = bigR, R.sq = R.sq, adj.R.sq = adj.R.sq, Std.Error.Est = see) 
  } else {
    model.summary <- data.frame(R = bigR, R.sq = R.sq, adj.R.sq = adj.R.sq, Std.Error.Est = see,
                                R.sq.change = R.sq.change, F.change = F.change, df1 = df1, df2 = df2, sig.F.change = p.change, cohens.f2.change = cohens.f2.change)
  }
  rownames(model.summary) <- paste0("Model ", seq_along(lm_sorted))

  coef_all <- lapply(lm_sorted, function(x) {
    cis <- confint(x)
    colnames(cis) <- c("Lower","Upper")
    out <- cbind(coef(summary(lm.beta::lm.beta(x))),
            cis,
            rbind(NA,cor_lm(x)),
            rbind(NA,col_diag(x)[,-1]))
    as.matrix(out)
  })
  
  if (dp%%1==0) {
    if (length(lm_outs)==1) {
      model.summary <- format(round(model.summary, dp), nsmall=dp)
    } else {
      model.summary <- format(round(model.summary[,c(1:6,9:10)], dp), nsmall=dp)
      model.summary$df1 = df1
      model.summary$df2 = df2
      model.summary <- model.summary[,c(1:6,9:10,7:8)]
    }
    anova_clean <- lapply(anova_all, function(x) {
      ssp <- format(round(x[2,"Sum of Sq"],dp), nsmall=dp)
      ssr <- format(round(x[2,"RSS"],dp), nsmall=dp)
      sst <- format(round(x[2,"Sum of Sq"] + x[2,"RSS"],dp), nsmall=dp)
      dfp <- x[2,"Df"]
      dfr <- x[2,"Res.Df"]
      dft <- dfp + dfr
      msp <- format(round(x[2,"Sum of Sq"]/dfp,dp), nsmall=dp)
      msr <- format(round(x[2,"RSS"]/dfr,dp), nsmall=dp)
      f <- format(round(x[2,"F"],dp), nsmall=dp)
      p <- format(round(x[2,"Pr(>F)"],dp), nsmall=dp)
      # cohens.f2 <- (ssp/(ssp+ssr)) / (1 - (ssp/(ssp+ssr)))
      cohens.f2 <- format(round((x[2,"Sum of Sq"]/(x[2,"Sum of Sq"]+x[2,"RSS"])) / (1 - (x[2,"Sum of Sq"]/(x[2,"Sum of Sq"]+x[2,"RSS"]))),dp), nsmall=dp)
      out <- data.frame(SS = c(ssp,ssr,sst), df = c(dfp,dfr,dft), MS = c(msp,msr,NA), F = c(f,NA,NA), sig = c(p,NA,NA), cohens.f2 = c(cohens.f2,NA,NA))
      rownames(out) <- c("Regression","Residual","Total")
      return(out)
    })
    coef_all <- lapply(lm_sorted, function(x) {
      cis <- confint(x)
      colnames(cis) <- c("Lower","Upper")
      cbind(coef(summary(lm.beta::lm.beta(x))),
                    cis,
                    rbind(NA,cor_lm(x)),
                    rbind(NA,col_diag(x)[,-1]))
      })
  }
  if (ci95==FALSE) {
    coef_all <- lapply(coef_all, function(x) x[,setdiff(colnames(x),c("Lower","Upper"))])
  }
  if (part==FALSE) {
    coef_all <- lapply(coef_all, function(x) x[,setdiff(colnames(x),c("zeroorder","partial","part"))])
  }
  if (coll==FALSE) {
    coef_all <- lapply(coef_all, function(x) x[,setdiff(colnames(x),c("Tolerance","VIF"))])
  }
  if (cf2==FALSE) {
    model.summary <- model.summary[,setdiff(colnames(model.summary),"cohens.f2.change")]
    anova_clean <- lapply(anova_clean, function(x) x[,setdiff(colnames(x),c("cohens.f2"))])
  }
  
  d <- psych::describe(data.final, skew=FALSE)
  r.output <- psych::corr.test(data.final)
  
  tables <- list(
    descriptives = d,
    correlations = list(
      r = r.output$r,
      n = r.output$n,
      p = r.output$p
    ),
    model_summary = model.summary,
    anova = anova_clean,
    coefficients = coef_all,
    dp = dp,
    mforms = m.form
  )
  
  class(tables) <- "regression3X10"
  return(tables)
}

print.regression3X10 <- function(x,...) {
  dp <- x$dp
  
  cat(strrep("-", 100), "\n")
  # PRINT DESCRIPTIVES
  cat("DESCRIPTIVE STATISTICS", "\n\n")
  # d <- psych::describe(data.final, skew=FALSE)
  if (dp%%1==0) {
    print(x$descriptives, digits=dp)
  } else {
    print(x$descriptives, digits=6)
  }
  cat(strrep("-", 100), "\n")
  
  # PRINT CORRELATIONS
  cat("CORRELATIONS", "\n\n")
  # r.output <- psych::corr.test(data.final)
  cat("Matrix\n")
  if (dp%%1==0) {
    print(round(x$correlations$r, dp))
  } else {
    print(x$correlations$r)
  }
  cat("\n"); cat("Sample Size\n")
  print(x$correlations$n)
  cat("\n"); cat("Sig.\n")
  if (dp%%1==0) {
    print(format(round(x$correlations$p, dp), nsmall = dp), quote=FALSE)
  } else {
    print(x$correlations$p)
  } 
  cat(strrep("-", 100), "\n")
  
  # PRINT MODEL SUMMARY TABLE
  cat("MODEL SUMMARY", "\n\n")
  for (i in seq_along(x$mforms)) {
    cat(sprintf("Model %d: %s\n", i, x$mforms[[i]]))
  }
  cat("\n")
  print(x$model_summary)
  cat(strrep("-", 100), "\n")
  
  # PRINT ANOVA TABLE
  cat("ANOVA", "\n\n")
  for (i in seq_along(x$anova)) {
    cat(sprintf("Model %d\n", i))
    print(x$anova[[i]], na.print = "")      
    if (i != length(seq_along(x$anova))) {
      cat("\n")
    }
  }
  cat(strrep("-", 100), "\n")    
  
  # PRINT COEFFICIENTS TABLE
  cat("COEFFICIENTS", "\n\n")
  for (i in seq_along(x$coefficients)) {
    cat(sprintf("Model %d\n", i))
    if (dp%%1==0) {
      format <- format(round(x$coefficients[[i]],dp),nsmall=dp)
      format[is.na(x$coefficients[[i]])] <- ""
      print(format, quote=FALSE)   
    } else {
      print(x$coefficients[[i]], na.print = "")      
    }
    if (i != length(seq_along(x$coefficients))) {
      cat("\n")
    }
  }
  cat(strrep("-", 100), "\n") 
}

mr_scatterplot <- function(x, y, data, type="simple", z=NULL, xlim=range(data[[x]]), ylim=range(data[[y]]), leg_pos="topleft", points=FALSE, point_transp=1, ...){
  if (type %in% c("interaction_categorical", "interaction_continuous") && is.null(z)) {
    stop(paste("Error: A 'z' variable must be provided for plot type '", type, "'.", sep = ""))
  }
  if (!is.null(z) && x == z) {
    stop('Error: The predictor variable "x" and the moderator variable "z" cannot be the same. Use type="quadratic" instead.')
  }
  plot(1, type = "n", xlim=xlim, ylim=ylim, xlab=x, ylab=y, ...)
  if (type=="simple"){
    if (points == TRUE) {
      points(data[[x]], data[[y]], col = adjustcolor("black", alpha.f = point_transp))
    }
    m <- lm(as.formula(paste(y, "~", x)), data = data)
    abline(m, col = "red", lwd = 3)
  }
  if (type=="interaction_categorical") {
    z_levels <- unique(data[[z]])
    for (i in seq_along(z_levels)) {
      d <- data[data[[z]] == z_levels[i], ]
      if(points==TRUE) {
        points(d[[x]], d[[y]], pch=i, col=adjustcolor(i, alpha.f = point_transp))
      }
      m <- lm(as.formula(paste(y, "~", x)), d)
      abline(m, col=i, lwd=2)
    }
    legend(leg_pos, legend = z_levels, col = seq_along(z_levels), pch = seq_along(z_levels), lwd = 2, title = z)
  }
  if (type=="interaction_continuous") {
    m <- lm(as.formula(paste(y, "~", x, "*", z)), data = data)
    mean_z <- mean(data[[z]], na.rm = TRUE)
    sd_z   <- sd(data[[z]], na.rm = TRUE)
    z_levels <- c(mean_z - sd_z, mean_z, mean_z + sd_z)
    z_labels <- c(paste("Low", z, "(-1 SD)"), 
                  paste("Mean", z), 
                  paste("High", z, "(+1 SD)"))
    cols <- c("blue", "red", "darkgreen")
    cols_points <- adjustcolor(cols, alpha.f = point_transp)
    shapes <- c(16, 17, 18)
    lines_lty <- c(2, 1, 4)
    if (points == TRUE) {
      cut_points <- c(-Inf, mean_z - 0.5 * sd_z, mean_z + 0.5 * sd_z, Inf)
      data$z_group <- cut(data[[z]], breaks = cut_points, labels = FALSE)
      points(data[[x]], data[[y]], 
             pch = shapes[data$z_group], 
             col = cols_points[data$z_group])
    }
    for (i in 1:3) {
      current_z <- z_levels[i]
      b <- coef(m)
      interaction_term_name <- names(b)[grepl(paste0(":", z), names(b)) | grepl(paste0(z, ":"), names(b))]
      simple_intercept <- b["(Intercept)"] + b[z] * current_z
      simple_slope     <- b[x] + b[interaction_term_name] * current_z
      abline(a = simple_intercept, b = simple_slope, col = cols[i], lty = lines_lty[i], lwd = 3)
    }
    legend(leg_pos, legend = z_labels, col = cols,  lty = lines_lty, lwd = 3)
  }
  if (type=="quadratic") {
    if (points == TRUE) {
      points(data[[x]], data[[y]], col = adjustcolor("black", alpha.f = point_transp))
    }
    fmla_string <- paste(y, "~", x, "+ I(", x, "^2)")
    m <- lm(as.formula(fmla_string), data = data)
    x_curve <- seq(min(data[[x]], na.rm = TRUE), max(data[[x]], na.rm = TRUE), length.out = 200)
    newdata_for_pred <- setNames(data.frame(x_curve), x)
    y_curve <- predict(m, newdata = newdata_for_pred)
    lines(x_curve, y_curve, col = "red", lwd = 3)
  }
}

aov_barplot <- function(ff, data, ylim=c(0, ceiling(max(m+se)*1.5)), l = ifelse(length(all.vars(ff))>2, T, F), ...){
  m = xtabs(ff, aggregate(ff, data, mean))
  se = xtabs(ff, aggregate(ff, data, FUN=function(x) sd(x)/sqrt(length(x))))
  bp = barplot(m, beside=T, ylim=ylim, ylab=paste("mean", all.vars(ff)[1]), xlab=all.vars(ff)[3], legend=l, xpd=F, ...); box()
  arrows(bp, m-se, bp, m+se, angle=90, length=0.2, code=3)
}

contrast_output <- function(mod, emm, contrast_coefs, ...) {
  out <- summary(contrast(emm, contrast_coefs, ...))
  out <- as.data.frame(out)
  model_call <- emm@model.info$call
  if (!is.null(mod$data$long)) {
    model_data <- mod$data$long
    grouping_vars <- names(emm@grid)[sapply(emm@grid, function(x) is.factor(x) || is.character(x))]
    cell_n <- as.data.frame(table(model_data[, grouping_vars, drop = FALSE]))
  } else {
    model_env <- parent.frame()
    model_obj <- eval(model_call, envir = model_env)
    model_data <- model.frame(model_obj)
    grouping_vars <- names(emm@grid)[sapply(emm@grid, function(x) is.factor(x) || is.character(x))]
    cell_n <- as.data.frame(table(model_data[, grouping_vars, drop = FALSE]))
  }
  ss_denominator <- sapply(contrast_coefs, function(x) sum(x^2 / cell_n$Freq))
  out$SS <- (out$estimate ^ 2) / ss_denominator
  out$F.ratio <- out$t.ratio ^ 2    # observed F for the contrasts
  out$partial.eta.sq <- (out$t.ratio ^ 2) / (out$t.ratio ^ 2 + out$df)   # partial eta squared for contrasts
  return(out)
}

crit_F_bonferroni <- function(alpha, k, v2) {
  qf(alpha/k, 1, v2, lower.tail=F)
}

crit_F_scheffe <- function(alpha, v1, v2) {
  v1 * qf(alpha, v1, v2, lower.tail=F)
}


anova_output <- function(x, dp=.5, ...) {
  aov_df <- as.data.frame(x$anova_table)
  aov_df$`Sum Sq` <- aov_df$F * aov_df$MSE * aov_df$`num Df`
  aov_df$`Error SS` <- aov_df$MSE * aov_df$`den Df`
  aov_df <- aov_df[, c("Sum Sq", "num Df", "Error SS", "den Df", "F", "Pr(>F)", "pes")]
  attr(aov_df, "dp") <- dp
  class(aov_df) <- "anova3X10"
  return(aov_df)
}

print.anova3X10 <- function(x,...) {
  dp <- attr(x, "dp")
  if (dp%%1==0) {
    cols_to_round <- setdiff(names(x), c("num Df", "den Df"))
    x[cols_to_round] <- lapply(x[cols_to_round], function(col) {
      formatted <- round(as.numeric(col), digits = dp)
      format(formatted, nsmall = dp)})
  }
  print(structure(x, class = "data.frame"), right = TRUE, ...)
}

