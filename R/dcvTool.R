# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

dcvTool <- function(call, data = NULL, x = NULL, y, tuning, cost = rmspe, 
        outerFolds, innerFolds, names = NULL, predictArgs = list(), 
        costArgs = list(),  envir = parent.frame()) {
    # initializations
    if(is.null(names)) {
        names <- if(is.null(data)) c("x", "y") else "data"
    }
    if(inherits(innerFolds, "cvFolds")) {
        # use the same inner folds for each outer fold
        # applicable if outer folds all have the same number of observations
        if(is.null(data)) {
            cvExpr <- expression(
                lapply(subsetList, dcvXY, folds=innerFolds, call=call, x=x, 
                    y=y, tuning=tuning, names=names, predictArgs=predictArgs, 
                    envir=envir)
            )
        } else {
            cvExpr <- expression(
                lapply(subsetList, dcvData, folds=innerFolds, call=call, 
                    data=data, y=y, tuning=tuning, names=names, 
                    predictArgs=predictArgs, envir=envir)
            )
        }
    } else {
        # 'innerFolds' is a list containing the inner folds for each outer fold
        if(is.null(data)) {
            cvExpr <- expression(
                mapply(dcvXY, subsetList, innerFolds, 
                    MoreArgs=list(call=call, x=x, y=y, tuning=tuning, 
                        names=names, predictArgs=predictArgs, envir=envir), 
                    SIMPLIFY=FALSE, USE.NAMES=FALSE)
            )
        } else {
            cvExpr <- expression(
                mapply(dcvData, subsetList, innerFolds, 
                    MoreArgs=list(call=call, data=data, y=y, tuning=tuning, 
                        names=names, predictArgs=predictArgs, envir=envir),
                    SIMPLIFY=FALSE, USE.NAMES=FALSE)
            )
        }
    }
    R <- outerFolds$R
    # perform (repeated) cross-validation
    if(R == 1) {
        subsetList <- getSubsetList(outerFolds)
        yHat <- eval(cvExpr)  # perform cross-validation
        # instead of collecting the results from the folds in the original 
        # order of the observations, the response is re-ordered accordingly
        yHat <- combineData(yHat)  # combine predictions from the folds
        y <- dataSubset(y, unlist(subsetList))  # re-order response
        # compute cost function for predicted values
        if(is.null(dim(y)) && !is.null(dim(yHat))) {
            cv <- apply(yHat, 2, 
                function(yHat) doCall(cost, y, yHat, args=costArgs))
            if(is.list(cv)) {
                cv <- list(sapply(cv, function(x) x[[1]]), 
                    sapply(cv, function(x) x[[2]]))
            }
        } else cv <- doCall(cost, y, yHat, args=costArgs)
    } else {
        cv <- sapply(seq_len(R), 
            function(r) {
                subsetList <- getSubsetList(outerFolds, r)
                yHat <- eval(cvExpr)  # perform cross-validation
                # instead of collecting the results from the folds in the original 
                # order of the observations, the response is re-ordered accordingly
                yHat <- combineData(yHat)  # combine predictions from the folds
                y <- dataSubset(y, unlist(subsetList))  # re-order response
                # compute cost function for predicted values
                if(is.null(dim(y)) && !is.null(dim(yHat))) {
                    tmp <- apply(yHat, 2, 
                        function(yHat) doCall(cost, y, yHat, args=costArgs))
                    if(is.list(tmp)) tmp <- sapply(tmp, function(x) x[[1]])
                } else {
                    tmp <- doCall(cost, y, yHat, args=costArgs)
                    if(is.list(tmp)) tmp <- tmp[[1]]
                }
                tmp
            })
    }
    # prepare and return results
    if(is.list(cv)) {
        cv <- lapply(cv, 
            function(x) {
                if(is.null(names(x))) names(x) <- defaultCvNames(length(x))
                x
            })
    } else {
        cv <- if(is.null(dim(cv)) && R > 1) as.matrix(cv) else t(cv)
        if(is.null(colnames(cv))) colnames(cv) <- defaultCvNames(ncol(cv))
    }
    cv
}

dcvXY <- function(i, folds, call, x, y, tuning, names, predictArgs, envir) {
    # extract training and test data
    xTraining <- dataSubset(x, -i)
    yTraining <- dataSubset(y, -i)
    xTest <- dataSubset(x, i)
    # perform CV to determine the optimal tuning parameters
    cv <- cvTuning(call, x=xTraining, y=yTraining, tuning=tuning, folds=folds, 
        names=names, predictArgs=predictArgs, envir=envir)
    # plug training data into function call
    call[[names[1]]] <- xTraining
    call[[names[2]]] <- yTraining
    # evaluate function call with optimal tuning parameters and predict 
    # response for test data
    best <- cv$best
    tuningDF <- cv$tuning
    tuningNames <- names(tuningDF)
    if(length(best) == 1) {
        for(j in seq_along(tuningNames)) {
            call[[tuningNames[j]]] <- tuningDF[best, j]
        }
        fit <- eval(call, envir)
        doCall(predict, fit, xTest, args=predictArgs)
    } else {
        idx <- split(seq_along(best), best)
        tmp <- mapply(function(best, idx) {
                for(j in seq_along(tuningNames)) {
                    call[[tuningNames[j]]] <- tuningDF[best, j]
                }
                fit <- eval(call, envir)
                doCall(predict, fit, xTest, args=predictArgs)[, idx, drop=FALSE]
            }, sort(unique(best)), idx, SIMPLIFY=FALSE)
        do.call(cbind, tmp)[, order(unlist(idx))]
    }
}

dcvData <- function(i, folds, call, data, y, tuning, names, predictArgs, envir) {
    # extract training and test data
    dataTraining <- dataSubset(data, -i)
    yTraining <- dataSubset(y, -i)
    dataTest <- dataSubset(data, i)
    # perform CV to determine the optimal tuning parameters
    cv <- cvTuning(call, data=dataTraining, y=yTraining, tuning=tuning, 
        folds=folds, names=names, predictArgs=predictArgs, envir=envir)
    # plug training data and optimal tuning parameters into function call
    call[[names]] <- dataTraining
    # evaluate function call with optimal tuning parameters and predict 
    # response for test data
    best <- cv$best
    tuningDF <- cv$tuning
    tuningNames <- names(tuningDF)
    if(length(best) == 1) {
        for(j in seq_along(tuningNames)) {
            call[[tuningNames[j]]] <- tuningDF[best, j]
        }
        fit <- eval(call, envir)
        doCall(predict, fit, dataTest, args=predictArgs)
    } else {
        idx <- split(seq_along(best), best)
        tmp <- mapply(function(best, idx) {
                for(j in seq_along(tuningNames)) {
                    call[[tuningNames[j]]] <- tuningDF[best, j]
                }
                fit <- eval(call, envir)
                doCall(predict, fit, dataTest, args=predictArgs)[, idx, drop=FALSE]
            }, sort(unique(best)), idx, SIMPLIFY=FALSE)
        do.call(cbind, tmp)[, order(unlist(idx))]
    }
}
