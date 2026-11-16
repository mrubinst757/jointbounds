# A diagnostic only: never modifies the learner's data or basis.
l2_predictor_check <- function(x, context, sink = NULL) {
  x <- as.data.frame(x)
  row <- data.frame(context=context,n_train=nrow(x),n_predictors=ncol(x),
    design_columns=NA_integer_,rank=NA_integer_,constant_predictors="",
    duplicate_columns="",dependent_columns="",status="ok",detail="",
    stringsAsFactors=FALSE)
  constant <- vapply(x,function(z) length(unique(z[!is.na(z)]))<=1L,logical(1))
  row$constant_predictors <- paste(names(x)[constant],collapse=", ")
  # Factors with one observed level cannot be contrast-coded. Record them as
  # redundant with the intercept and exclude them only from this diagnostic.
  single_factor <- constant & vapply(x,function(z) is.factor(z)||is.character(z),logical(1))
  build <- tryCatch({
    xx <- x[,!single_factor,drop=FALSE]
    B <- if (!ncol(xx)) matrix(1,nrow(x),1,dimnames=list(NULL,"(Intercept)")) else
      stats::model.matrix(~ .,data=stats::model.frame(~ .,xx,na.action=stats::na.fail))
    if (!nrow(B) || any(!is.finite(B))) stop("empty or nonfinite design matrix")
    # Scale columns for the numerical rank check; preserve the intercept.
    scales <- apply(abs(B),2,max); scales[scales==0] <- 1
    q <- qr(sweep(B,2,scales,"/"),tol=1e-7,LAPACK=FALSE)
    dep <- if(q$rank<ncol(B)) colnames(B)[q$pivot[seq.int(q$rank+1L,ncol(B))]] else character()
    list(columns=ncol(B)+sum(single_factor),rank=q$rank,
      duplicate=colnames(B)[duplicated(t(B),MARGIN=1L)],
      dependent=c(dep,names(x)[single_factor]))
  },error=identity)
  if (inherits(build,"error")) {
    row$status <- "unchecked"; row$detail <- conditionMessage(build)
  } else {
    row$design_columns <- build$columns; row$rank <- build$rank
    row$duplicate_columns <- paste(build$duplicate,collapse=", ")
    row$dependent_columns <- paste(build$dependent,collapse=", ")
    if(length(build$dependent)||any(constant)) row$status <- "redundant"
  }
  if (!is.null(sink)) sink$rows[[length(sink$rows)+1L]] <- row
  if (row$status=="redundant") {
    items <- c(if(nzchar(row$constant_predictors)) paste0("constant predictors: ",row$constant_predictors),
      if(nzchar(row$duplicate_columns)) paste0("duplicate design columns: ",row$duplicate_columns),
      if(nzchar(row$dependent_columns)) paste0("dependent design columns: ",row$dependent_columns))
    warning(context,": predictor redundancy detected (design rank ",row$rank,"/",row$design_columns,
      "; n_train=",nrow(x),"). ",paste(items,collapse="; "),
      ". Predictors were not removed; GLM predictions may be non-estimable for some new rows.",call.=FALSE)
  } else if(row$status=="unchecked") warning(context,": predictor redundancy check could not be completed: ",
    row$detail,call.=FALSE)
  invisible(row)
}

l2_predictor_sink <- function() {
  z <- new.env(parent=emptyenv()); z$rows <- list(); z
}

l2_predictor_table <- function(sink) {
  if (!length(sink$rows)) return(data.frame())
  unique(do.call(rbind,sink$rows))
}

l2_add_predictor_table <- function(sink,table) {
  if (!is.null(table) && nrow(table)) sink$rows[[length(sink$rows)+1L]] <- table
  invisible(NULL)
}

# Check the actual screened designs used inside SL.glm's internal CV fits.
l2_checked_glm_wrapper <- function(original,context,sink) {
  force(original); force(context); force(sink); counter <- 0L
  function(Y,X,newX,...) {
    counter <<- counter+1L
    l2_predictor_check(X,paste0(context,"; SL.glm internal fit ",counter),sink)
    original(Y=Y,X=X,newX=newX,...)
  }
}
