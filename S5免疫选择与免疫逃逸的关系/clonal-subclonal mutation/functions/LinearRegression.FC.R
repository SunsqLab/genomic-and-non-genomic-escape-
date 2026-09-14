# Original adjusted linear-regression function and coefficient extraction.
LinearRegression.FC <- function(data,
                                outcome,
                                exposure,
                                adjust.vars = NULL,
                                conf.level = 0.95,
                                reference.levels = NULL,
                                return.vars = exposure,
                                outcome.log2 = TRUE,
                                digits = NULL) {
    quote_name <- function(x) {
        paste0("`", gsub("`", "``", x, fixed = TRUE), "`")
    }

    if (!is.data.frame(data)) {
        stop("data must be a data.frame.", call. = FALSE)
    }
    if (!is.character(outcome) || length(outcome) != 1) {
        stop("outcome must be a single column name.", call. = FALSE)
    }
    if (!is.character(exposure) || length(exposure) == 0) {
        stop("exposure must be a non-empty character vector.", call. = FALSE)
    }
    if (!is.null(adjust.vars) && !is.character(adjust.vars)) {
        stop("adjust.vars must be NULL or a character vector.", call. = FALSE)
    }
    if (!is.numeric(conf.level) || length(conf.level) != 1 || conf.level <= 0 || conf.level >= 1) {
        stop("conf.level must be a number between 0 and 1.", call. = FALSE)
    }
    if (!is.logical(outcome.log2) || length(outcome.log2) != 1) {
        stop("outcome.log2 must be TRUE or FALSE.", call. = FALSE)
    }
    if (!is.null(digits) && (!is.numeric(digits) || length(digits) != 1 || digits < 0)) {
        stop("digits must be NULL or a non-negative number.", call. = FALSE)
    }
    if (!is.null(reference.levels) && is.null(names(reference.levels))) {
        stop("reference.levels must be a named character vector or named list.", call. = FALSE)
    }

    predictors <- c(exposure, adjust.vars)
    if (!is.character(return.vars) || length(return.vars) == 0) {
        stop("return.vars must be a non-empty character vector or \"all\".", call. = FALSE)
    }
    if (!identical(return.vars, "all") && !all(return.vars %in% predictors)) {
        stop(sprintf("return.vars contains variables not included in the model: %s", paste(setdiff(return.vars, predictors), collapse = ", ")), call. = FALSE)
    }
    if (anyDuplicated(predictors)) {
        stop("exposure and adjust.vars contain duplicate column names; remove duplicates before fitting.", call. = FALSE)
    }
    if (outcome %in% predictors) {
        stop("outcome cannot also be included in exposure or adjust.vars.", call. = FALSE)
    }

    model.cols <- unique(c(outcome, predictors))
    missing.cols <- setdiff(model.cols, colnames(data))
    if (length(missing.cols) > 0) {
        stop(sprintf("Missing columns: %s", paste(missing.cols, collapse = ", ")), call. = FALSE)
    }

    model.data <- data[, model.cols, drop = FALSE]
    model.data <- model.data[stats::complete.cases(model.data), , drop = FALSE]
    if (nrow(model.data) == 0) {
        stop("No complete observations remain for model fitting.", call. = FALSE)
    }
    if (!is.numeric(model.data[[outcome]]) && !is.integer(model.data[[outcome]])) {
        stop("outcome must be a continuous numeric variable; use logistic regression for a binary outcome.", call. = FALSE)
    }
    if (length(unique(model.data[[outcome]])) < 2) {
        stop("outcome has fewer than two distinct values in complete observations.", call. = FALSE)
    }

    for (predictor in predictors) {
        if (is.character(model.data[[predictor]])) {
            model.data[[predictor]] <- factor(model.data[[predictor]])
        }
        if (is.logical(model.data[[predictor]])) {
            model.data[[predictor]] <- factor(model.data[[predictor]], levels = c(FALSE, TRUE))
        }
        if (is.factor(model.data[[predictor]])) {
            model.data[[predictor]] <- droplevels(model.data[[predictor]])
            if (!is.null(reference.levels) && predictor %in% names(reference.levels)) {
                ref.level <- as.character(reference.levels[[predictor]])[1]
                if (!ref.level %in% levels(model.data[[predictor]])) {
                    stop(sprintf("Reference level '%s' was not found in variable '%s'.", ref.level, predictor), call. = FALSE)
                }
                model.data[[predictor]] <- stats::relevel(model.data[[predictor]], ref = ref.level)
            }
            if (nlevels(model.data[[predictor]]) < 2) {
                stop(sprintf("Variable '%s' has fewer than two levels in complete observations.", predictor), call. = FALSE)
            }
        } else if (length(unique(model.data[[predictor]])) < 2) {
            stop(sprintf("Variable '%s' has fewer than two distinct values in complete observations.", predictor), call. = FALSE)
        }
    }

    model.formula <- stats::as.formula(
        paste(quote_name(outcome), "~", paste(vapply(predictors, quote_name, character(1)), collapse = " + "))
    )
    fit <- stats::lm(model.formula, data = model.data)
    if (fit$rank < ncol(stats::model.matrix(fit))) {
        warning("The model matrix is rank deficient; lm may have omitted collinear terms.", call. = FALSE)
    }

    coef.table <- stats::coef(summary(fit))
    coef.table <- coef.table[rownames(coef.table) != "(Intercept)", , drop = FALSE]
    if (nrow(coef.table) == 0) {
        stop("The model has no predictor coefficients to report.", call. = FALSE)
    }

    model.matrix <- stats::model.matrix(fit)
    assign.index <- attr(model.matrix, "assign")
    coef.index <- match(rownames(coef.table), colnames(model.matrix))
    variable.name <- predictors[assign.index[coef.index]]

    alpha <- 1 - conf.level
    t.value <- stats::qt(1 - alpha / 2, df = fit$df.residual)
    estimate <- coef.table[, "Estimate"]
    std.error <- coef.table[, "Std. Error"]
    ci.min <- estimate - t.value * std.error
    ci.max <- estimate + t.value * std.error

    global.p.value <- rep(NA_real_, length(predictors))
    names(global.p.value) <- predictors
    drop.table <- tryCatch(stats::drop1(fit, test = "F"), error = function(e) NULL)
    if (!is.null(drop.table) && "Pr(>F)" %in% colnames(drop.table)) {
        term.labels <- attr(stats::terms(fit), "term.labels")
        for (i in seq_along(predictors)) {
            possible.names <- unique(c(term.labels[i], predictors[i], quote_name(predictors[i])))
            matched.name <- intersect(possible.names, rownames(drop.table))
            if (length(matched.name) > 0) {
                global.p.value[predictors[i]] <- drop.table[matched.name[1], "Pr(>F)"]
            }
        }
    }

    result <- data.frame(
        Variable = unname(variable.name),
        Term = rownames(coef.table),
        Estimate = estimate,
        CI.min = ci.min,
        CI.max = ci.max,
        p.value = coef.table[, "Pr(>|t|)"],
        global.p.value = unname(global.p.value[variable.name]),
        FoldChange = if (outcome.log2) 2^estimate else NA_real_,
        FC.CI.min = if (outcome.log2) 2^ci.min else NA_real_,
        FC.CI.max = if (outcome.log2) 2^ci.max else NA_real_,
        N = nrow(model.data),
        Outcome = outcome,
        Reference = vapply(variable.name, function(x) {
            if (is.factor(model.data[[x]])) {
                levels(model.data[[x]])[1]
            } else {
                NA_character_
            }
        }, character(1)),
        stringsAsFactors = FALSE,
        row.names = NULL
    )

    if (!identical(return.vars, "all")) {
        result <- result[result$Variable %in% return.vars, , drop = FALSE]
    }

    if (!is.null(digits)) {
        numeric.cols <- c("Estimate", "CI.min", "CI.max", "p.value", "global.p.value", "FoldChange", "FC.CI.min", "FC.CI.max")
        result[numeric.cols] <- lapply(result[numeric.cols], round, digits = digits)
    }

    result
}
