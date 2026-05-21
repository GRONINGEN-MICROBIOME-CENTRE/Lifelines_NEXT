prep_glmmPen_data <- function(data,
                              model_variables,
                              responses = c("Factor1", "Factor2", "Factor3"),
                              time_vars = c("age_c", "age_c2", "age_c3"),
                              sample_id_var = "NG_ID",
                              individual_id_var = "NEXT_ID",
                              center_time = FALSE) {
  
  # one-sided formula: only RHS
  rhs <- paste(model_variables, collapse = " + ")
  ff  <- stats::as.formula(paste("~", rhs))
  
  # fixed-effects model matrix (includes intercept)
  X <- stats::model.matrix(ff, data = data)
  
  cn <- colnames(X)
  is_intercept <- cn == "(Intercept)"
  is_time_col  <- cn %in% time_vars
  
  # scale everything except intercept + time terms
  scale_these <- !(is_intercept | is_time_col)
  
  Xs <- X
  Xs[, scale_these] <- scale(Xs[, scale_these, drop = FALSE],
                             center = TRUE, scale = TRUE)
  
  # Optional: enforce centering time columns (ONLY if you want this function to do it)
  if (center_time) {
    Xs[, is_time_col] <- scale(Xs[, is_time_col, drop = FALSE],
                               center = TRUE, scale = FALSE)
  }
  
  # Build output:
  # - keep both IDs
  # - add response columns so the same df can be used for any response
  # - add standardized model matrix columns
  out <- data.frame(
    NG_ID   = data[[sample_id_var]],
    NEXT_ID = data[[individual_id_var]],
    data[, responses, drop = FALSE],
    Xs,
    check.names = FALSE
  )
  
  out1 <- out[, colnames(out) != "(Intercept)"]
  
  return(out1)
}

data_scaled <- prep_glmmPen_data(data = data_complete_cases, model_variables = model_variables, center_time = FALSE)


