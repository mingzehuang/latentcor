############################################################################################
# For multilinear interpolation approximation for bridge Inverse
############################################################################################

#' @importFrom chebpol ipol
#'
#'
NULL


############################################################################################
# Cutoff criteria based on the combination of variable types
############################################################################################
# For ml only
# cutoff_bc <- function(zratio1, zratio2 = NULL){rep(1, length(zratio1))}
# For mlbd
cutoff_bc <- function(zratio1, zratio2 = NULL){0.9 * 2 * zratio1 * (1 - zratio1)}
# For ml only
# cutoff_cb <- function(zratio1 = NULL, zratio2){rep(1, length(zratio2))}
# For mlbd
cutoff_cb <- function(zratio1 = NULL, zratio2){0.9 * 2 * zratio2 * (1 - zratio2)}
# For ml only
# cutoff_bb <- function(zratio1, zratio2){rep(1, length(zratio1))}
# For mlbd
cutoff_bb <- function(zratio1, zratio2){0.9 * 2 * pmin(zratio1, zratio2)*(1-pmax(zratio1, zratio2))}
# For ml only
# cutoff_tc <- function(zratio1, zratio2 = NULL){rep(1, length(zratio1))}
# For mlbd
cutoff_tc <- function(zratio1, zratio2 = NULL){0.9 * (1 - zratio1^2)}
# For ml only
# cutoff_ct <- function(zratio1 = NULL, zratio2){rep(1, length(zratio2))}
# For mlbd
cutoff_ct <- function(zratio1 = NULL, zratio2){0.9 * (1 - zratio2^2)}
# For ml only
# cutoff_tb <- function(zratio1, zratio2){rep(1, length(zratio1))}
# For mlbd
cutoff_tb <- function(zratio1, zratio2){0.9 * 2 * pmax(zratio2, 1 - zratio2) * (1 - pmax(zratio2, 1 - zratio2, zratio1))}
# For ml only
# cutoff_bt <- function(zratio1, zratio2){rep(1, length(zratio1))}
# For mlbd
cutoff_bt <- function(zratio1, zratio2){0.9 * 2 * pmax(zratio1, 1 - zratio1) * (1 - pmax(zratio1, 1 - zratio1, zratio2))}
# For ml only
# cutoff_tt <- function(zratio1, zratio2){rep(1, length(zratio1))}
# For mlbd
cutoff_tt <- function(zratio1, zratio2){0.9 * (1 - pmax(zratio1, zratio2)^2)}
# For ml only
# cutoff_nc <- function(zratio1, zratio2 = NULL){rep(1, nrow(zratio1))}
# For mlbd
cutoff_nc <- function(zratio1, zratio2 = NULL){0.9 * 2 * (zratio1[ , 1] * (zratio1[ , 2] - zratio1[ , 1]) + (1 - zratio1[ , 2]) * zratio1[ , 2])}
# For ml only
# cutoff_cn <- function(zratio1 = NULL, zratio2){rep(1, nrow(zratio2))}
# For mlbd
cutoff_cn <- function(zratio1 = NULL, zratio2){0.9 * 2 * (zratio2[ , 1] * (zratio2[ , 2] - zratio2[ , 1]) + (1 - zratio2[ , 2]) * zratio2[ , 2])}
# For ml only
# cutoff_nb <- function(zratio1, zratio2){rep(1, nrow(zratio1))}
# For mlbd
cutoff_nb <- function(zratio1, zratio2){0.9 * 2 * pmin(zratio1[ , 1] * (zratio1[ , 2] - zratio1[ , 1]) + (1 - zratio1[ , 2]) * zratio1[ , 2], zratio2 * (1 - zratio2))}
# For ml only
# cutoff_bn <- function(zratio1, zratio2){rep(1, nrow(zratio1))}
# For mlbd
cutoff_bn <- function(zratio1, zratio2){0.9 * 2 * pmin(zratio2[ , 1] * (zratio2[ , 2] - zratio2[ , 1]) + (1 - zratio2[ , 2]) * zratio2[ , 2], zratio1 * (1 - zratio1))}
# For ml only
# cutoff_nn <- function(zratio1, zratio2){rep(1, nrow(zratio1))}
# For mlbd
cutoff_nn <- function(zratio1, zratio2){0.9 * 2 * pmin(zratio1[ , 1] * (zratio1[ , 2] - zratio1[ , 1]) + (1 - zratio1[ , 2]) * zratio1[ , 2],
                                                       zratio2[ , 1] * (zratio2[ , 2] - zratio2[ , 1]) + (1 - zratio2[ , 2]) * zratio2[ , 2])}

cutoff_select <- function(type1, type2){
  if (type1 == "binary" & type2 == "binary") {
    cutoff_select <- cutoff_bb
  } else if (type1 == "trunc" & type2 == "trunc") {
    cutoff_select <- cutoff_tt
  } else if (type1 == "trunc" & type2 == "continuous") {
    cutoff_select <- cutoff_tc
  } else if (type1 == "continuous" & type2 == "trunc") {
    cutoff_select <- cutoff_ct
  } else if (type1 == "binary" & type2 == "continuous") {
    cutoff_select <- cutoff_bc
  } else if (type1 == "continuous" & type2 == "binary") {
    cutoff_select <- cutoff_cb
  } else if (type1 == "trunc" & type2 == "binary") {
    cutoff_select <- cutoff_tb
  } else if (type1 == "binary" & type2 == "trunc") {
    cutoff_select <- cutoff_bt
  } else if (type1 == "ternary" & type2 == "continuous") {
    cutoff_select <- cutoff_nc
  } else if (type1 == "continuous" & type2 == "ternary") {
    cutoff_select <- cutoff_cn
  } else if (type1 == "ternary" & type2 == "binary") {
    cutoff_select <- cutoff_nb
  } else if (type1 == "binary" & type2 == "ternary") {
    cutoff_select <- cutoff_bn
  } else if (type1 == "ternary" & type2 == "ternary") {
    cutoff_select <- cutoff_nn
  } else {
    stop("Unrecognized type of variables. Should be one of continuous, binary or trunc.")
  }
}



############################################################################################
# Select which bridge inverse function based on the combinatino of variable types
############################################################################################

bridgeInv_select <- function(type1, type2) {
  if (type1 == "binary" & type2 == "binary") { bridgeInv_select <- bridgeInv_bb
  } else if (type1 == "trunc" & type2 == "trunc") { bridgeInv_select <- bridgeInv_tt
  } else if (type1 == "trunc" & type2 == "continuous") { bridgeInv_select <- bridgeInv_tc
  } else if (type1 == "continuous" & type2 == "trunc") { bridgeInv_select <- bridgeInv_ct
  } else if (type1 == "binary" & type2 == "continuous") { bridgeInv_select <- bridgeInv_bc
  } else if (type1 == "continuous" & type2 == "binary") { bridgeInv_select <- bridgeInv_cb
  } else if (type1 == "trunc" & type2 == "binary") { bridgeInv_select <- bridgeInv_tb
  } else if (type1 == "binary" & type2 == "trunc") { bridgeInv_select <- bridgeInv_bt
  } else if (type1 == "ternary" & type2 == "binary") {bridgeInv_select <- bridgeInv_nb
  } else if (type1 == "binary" & type2 == "ternary") {bridgeInv_select <- bridgeInv_bn
  } else if (type1 == "ternary" & type2 == "continuous") {bridgeInv_select <- bridgeInv_nc
  } else if (type1 == "continuous" & type2 == "ternary") {bridgeInv_select <- bridgeInv_cn
  } else if (type1 == "ternary" & type2 == "ternary") {bridgeInv_select <- bridgeInv_nn
  } else {
    stop("Unrecognized type of variables. Should be one of continuous, binary or trunc.")
  }
}


# wrapper function for BC
bridgeInv_bc <- function(tau, zratio1, zratio2 = NULL){
  out <- BCipol(rbind(t(tau), t(zratio1))) / 10^7
  return(out)
}

# wrapper function for CB
bridgeInv_cb <- function(tau, zratio1 = NULL, zratio2){
  out <- BCipol(rbind(t(tau), t(zratio2))) / 10^7
  return(out)
}

# wrapper function
bridgeInv_bb <- function(tau, zratio1, zratio2){
  out <- BBipol(rbind(t(tau), t(zratio1), t(zratio2))) / 10^7
  return(out)
}

# wrapper functions
bridgeInv_tc <- function(tau, zratio1, zratio2 = NULL){
  out <- TCipol(rbind(t(tau), t(zratio1))) / 10^7
  return(out)
}

bridgeInv_ct <- function(tau, zratio1 = NULL, zratio2){
  out <- TCipol(rbind(t(tau), t(zratio2))) / 10^7
  return(out)
}

# wrapper functions
bridgeInv_tb <- function(tau, zratio1, zratio2){
  out <- TBipol(rbind(t(tau), t(zratio1), t(zratio2))) / 10^7
  return(out)
}

bridgeInv_bt <- function(tau, zratio1, zratio2){
  out <- TBipol(rbind(t(tau), t(zratio2), t(zratio1))) / 10^7
  return(out)
}

# wrapper function
bridgeInv_tt <- function(tau, zratio1, zratio2){
  out <- TTipol(rbind(t(tau), t(zratio1), t(zratio2))) / 10^7
  return(out)
}

# wrapper function
bridgeInv_nc <- function(tau, zratio1, zratio2 = NULL){
  out <- NCipol(rbind(t(tau), t(zratio1[ , 1] / zratio1[ , 2]), t(zratio1[ , 2]))) / 10^7
  return(out)
}

# wrapper function
bridgeInv_cn <- function(tau, zratio1 = NULL, zratio2){
  out <- NCipol(rbind(t(tau), t(zratio2[ , 1] / zratio2[ , 2]), t(zratio2[ , 2]))) / 10^7
  return(out)
}

# wrapper function
bridgeInv_nb <- function(tau, zratio1, zratio2){
  out <- NBipol(rbind(t(tau), t(zratio1[ , 1] / zratio1[ , 2]), t(zratio1[ , 2]), t(zratio2))) / 10^7
  return(out)
}

# wrapper function
bridgeInv_bn <- function(tau, zratio1, zratio2){
  out <- NBipol(rbind(t(tau), t(zratio2[ , 1] / zratio2[ , 2]), t(zratio2[ , 2]), t(zratio1))) / 10^7
  return(out)
}

# wrapper function
bridgeInv_nn <- function(tau, zratio1, zratio2){
  out <- NNipol(rbind(t(tau), t(zratio1[ , 1] / zratio1[ , 2]), t(zratio1[ , 2]), t(zratio2[ , 1] / zratio2[ , 2]), t(zratio2[ , 2]))) / 10^7
  return(out)
}

