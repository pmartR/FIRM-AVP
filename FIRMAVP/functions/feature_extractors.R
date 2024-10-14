#-------feature extractions------#

# amino acid composition

#' Amino Acid Composition Descriptor
#'
#' This function calculates the Amino Acid Composition descriptor (dim: 20).
#'
#' @param x A character vector, as the input protein sequence.
#'
#' @return A length 20 named vector
#'
#require(Biostrings)
require(DECIPHER)

extractAAC_revised <- function(x) {
  
  AADict <- c(
    "A", "R", "N", "D", "C", "E", "Q", "G", "H", "I",
    "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V" 
  )
  
  if (all(strsplit(x, split = "")[[1]] %in% AADict) == FALSE) {
    stop("x has unrecognized amino acid type")
  }
  
  # 20 Amino Acid Abbrevation Dictionary from
  # https://en.wikipedia.org/wiki/Amino_acid#Table_of_standard_amino_acid_abbreviations_and_properties
  
  
  
  AAC <- summary(
    factor(strsplit(x, split = "")[[1]], levels = AADict),
    maxsum = 21
  ) / nchar(x)
  
  return (AAC)
}


#' Dipeptide Composition Descriptor
#'
#' This function calculates the Dipeptide Composition descriptor (dim: 400).

extractDC_revised <- function(x) {
  
  protcheck <- function(x) {
    AADict <- c(
      "A", "R", "N", "D", "C", "E", "Q", "G", "H", "I",
      "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"
    )
    
    all(strsplit(x, split = "")[[1]] %in% AADict)
  }
  
  if (protcheck(x) == FALSE) {
    stop("x has unrecognized amino acid type")
  }
  
  # To achieve maximum performance, here we use dictionary directly
  # DCDict could also be generated with
  # AADict = c('A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I',
  #            'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')
  # as.vector((outer(AADict, AADict, paste, sep = '')))
  
  DCDict <- c(
    "AA", "RA", "NA", "DA", "CA", "EA", "QA", "GA", "HA", "IA",
    "LA", "KA", "MA", "FA", "PA", "SA", "TA", "WA", "YA", "VA",
    "AR", "RR", "NR", "DR", "CR", "ER", "QR", "GR", "HR", "IR",
    "LR", "KR", "MR", "FR", "PR", "SR", "TR", "WR", "YR", "VR",
    "AN", "RN", "NN", "DN", "CN", "EN", "QN", "GN", "HN", "IN",
    "LN", "KN", "MN", "FN", "PN", "SN", "TN", "WN", "YN", "VN",
    "AD", "RD", "ND", "DD", "CD", "ED", "QD", "GD", "HD", "ID",
    "LD", "KD", "MD", "FD", "PD", "SD", "TD", "WD", "YD", "VD",
    "AC", "RC", "NC", "DC", "CC", "EC", "QC", "GC", "HC", "IC",
    "LC", "KC", "MC", "FC", "PC", "SC", "TC", "WC", "YC", "VC",
    "AE", "RE", "NE", "DE", "CE", "EE", "QE", "GE", "HE", "IE",
    "LE", "KE", "ME", "FE", "PE", "SE", "TE", "WE", "YE", "VE",
    "AQ", "RQ", "NQ", "DQ", "CQ", "EQ", "QQ", "GQ", "HQ", "IQ",
    "LQ", "KQ", "MQ", "FQ", "PQ", "SQ", "TQ", "WQ", "YQ", "VQ",
    "AG", "RG", "NG", "DG", "CG", "EG", "QG", "GG", "HG", "IG",
    "LG", "KG", "MG", "FG", "PG", "SG", "TG", "WG", "YG", "VG",
    "AH", "RH", "NH", "DH", "CH", "EH", "QH", "GH", "HH", "IH",
    "LH", "KH", "MH", "FH", "PH", "SH", "TH", "WH", "YH", "VH",
    "AI", "RI", "NI", "DI", "CI", "EI", "QI", "GI", "HI", "II",
    "LI", "KI", "MI", "FI", "PI", "SI", "TI", "WI", "YI", "VI",
    "AL", "RL", "NL", "DL", "CL", "EL", "QL", "GL", "HL", "IL",
    "LL", "KL", "ML", "FL", "PL", "SL", "TL", "WL", "YL", "VL",
    "AK", "RK", "NK", "DK", "CK", "EK", "QK", "GK", "HK", "IK",
    "LK", "KK", "MK", "FK", "PK", "SK", "TK", "WK", "YK", "VK",
    "AM", "RM", "NM", "DM", "CM", "EM", "QM", "GM", "HM", "IM",
    "LM", "KM", "MM", "FM", "PM", "SM", "TM", "WM", "YM", "VM",
    "AF", "RF", "NF", "DF", "CF", "EF", "QF", "GF", "HF", "IF",
    "LF", "KF", "MF", "FF", "PF", "SF", "TF", "WF", "YF", "VF",
    "AP", "RP", "NP", "DP", "CP", "EP", "QP", "GP", "HP", "IP",
    "LP", "KP", "MP", "FP", "PP", "SP", "TP", "WP", "YP", "VP",
    "AS", "RS", "NS", "DS", "CS", "ES", "QS", "GS", "HS", "IS",
    "LS", "KS", "MS", "FS", "PS", "SS", "TS", "WS", "YS", "VS",
    "AT", "RT", "NT", "DT", "CT", "ET", "QT", "GT", "HT", "IT",
    "LT", "KT", "MT", "FT", "PT", "ST", "TT", "WT", "YT", "VT",
    "AW", "RW", "NW", "DW", "CW", "EW", "QW", "GW", "HW", "IW",
    "LW", "KW", "MW", "FW", "PW", "SW", "TW", "WW", "YW", "VW",
    "AY", "RY", "NY", "DY", "CY", "EY", "QY", "GY", "HY", "IY",
    "LY", "KY", "MY", "FY", "PY", "SY", "TY", "WY", "YY", "VY",
    "AV", "RV", "NV", "DV", "CV", "EV", "QV", "GV", "HV", "IV",
    "LV", "KV", "MV", "FV", "PV", "SV", "TV", "WV", "YV", "VV"
  )
  
  xSplitted <- strsplit(x, split = "")[[1]]
  n <- nchar(x)
  DC <- summary(factor(
    paste(xSplitted[-n], xSplitted[-1], sep = ""),
    levels = DCDict
  ), maxsum = 401) / (n - 1)
  
  return(DC)
}


#' Pseudo Amino Acid Composition (PseAAC) Descriptor
#'
#' This function calculates the Pseudo Amino Acid Composition (PseAAC)

extractPAAC_revised <- function(
  x, kayes, props = c("Hydrophobicity", "Hydrophilicity", "SideChainMass"),
  lambda = 5, w = 0.05, customprops = NULL) {
  
  protcheck <- function(x) {
    AADict <- c(
      "A", "R", "N", "D", "C", "E", "Q", "G", "H", "I",
      "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"
    )
    
    all(strsplit(x, split = "")[[1]] %in% AADict)
  }
  
  if (protcheck(x) == FALSE) {
    stop("x has unrecognized amino acid type")
  }
  
  if (nchar(x) <= lambda) {
    stop('Length of the protein sequence must be greater than "lambda"')
  }
  
  AAidx <- read.csv(kayes, header = TRUE)
  
  tmp <- data.frame(
    AccNo = c("Hydrophobicity", "Hydrophilicity", "SideChainMass"),
    A = c(0.62, -0.5, 15), R = c(-2.53, 3, 101),
    N = c(-0.78, 0.2, 58), D = c(-0.9, 3, 59),
    C = c(0.29, -1, 47), E = c(-0.74, 3, 73),
    Q = c(-0.85, 0.2, 72), G = c(0.48, 0, 1),
    H = c(-0.4, -0.5, 82), I = c(1.38, -1.8, 57),
    L = c(1.06, -1.8, 57), K = c(-1.5, 3, 73),
    M = c(0.64, -1.3, 75), F = c(1.19, -2.5, 91),
    P = c(0.12, 0, 42), S = c(-0.18, 0.3, 31),
    T = c(-0.05, -0.4, 45), W = c(0.81, -3.4, 130),
    Y = c(0.26, -2.3, 107), V = c(1.08, -1.5, 43)
  )
  
  AAidx <- rbind(AAidx, tmp)
  
  if (!is.null(customprops)) AAidx <- rbind(AAidx, customprops)
  
  aaidx <- AAidx[, -1]
  row.names(aaidx) <- AAidx[, 1]
  
  n <- length(props)
  
  # Standardize H0 to H
  
  H0 <- as.matrix(aaidx[props, ])
  
  H <- matrix(ncol = 20, nrow = n)
  for (i in 1:n) {
    H[i, ] <-
      (H0[i, ] - mean(H0[i, ])) / (sqrt(sum((H0[i, ] - mean(H0[i, ]))^2) / 20))
  }
  AADict <- c(
    "A", "R", "N", "D", "C", "E", "Q", "G", "H", "I",
    "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"
  )
  dimnames(H) <- list(props, AADict)
  
  # Compute (big) Theta
  
  Theta <- vector("list", lambda)
  
  xSplitted <- strsplit(x, split = "")[[1]]
  
  N <- length(xSplitted)
  
  for (i in 1:lambda) {
    for (j in 1:(N - i)) {
      Theta[[i]][j] <- mean((H[, xSplitted[j]] - H[, xSplitted[j + i]])^2)
    }
  }
  
  # Compute (small) theta
  
  theta <- sapply(Theta, mean)
  
  # Compute first 20 features
  
  fc <- summary(factor(xSplitted, levels = AADict), maxsum = 21)
  Xc1 <- fc / (1 + (w * sum(theta)))
  names(Xc1) <- paste("Xc1.", names(Xc1), sep = "")
  
  # Compute last lambda features
  
  Xc2 <- (w * theta) / (1 + (w * sum(theta)))
  names(Xc2) <- paste("Xc2.lambda.", 1:lambda, sep = "")
  
  # Combine (20 + lambda) features
  
  Xc <- c(Xc1, Xc2)
  
  return (Xc)
}

#' Amphiphilic Pseudo Amino Acid Composition (APseAAC) Descriptor
#'
#' This function calculates the Amphiphilic Pseudo Amino Acid

extractAPAAC_revised <- function(
  x, kayes, props = c("Hydrophobicity", "Hydrophilicity"),
  lambda = 5, w = 0.05, customprops = NULL) {
  
  protcheck <- function(x) {
    AADict <- c(
      "A", "R", "N", "D", "C", "E", "Q", "G", "H", "I",
      "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"
    )
    
    all(strsplit(x, split = "")[[1]] %in% AADict)
  }
  
  if (protcheck(x) == FALSE) {
    stop("x has unrecognized amino acid type")
  }
  
  if (nchar(x) <= lambda) {
    stop('Length of the protein sequence must be greater than "lambda"')
  }
  
  AAidx <- read.csv(kayes, header = TRUE)
  
  tmp <- data.frame(
    AccNo = c("Hydrophobicity", "Hydrophilicity", "SideChainMass"),
    A = c(0.62, -0.5, 15), R = c(-2.53, 3, 101),
    N = c(-0.78, 0.2, 58), D = c(-0.9, 3, 59),
    C = c(0.29, -1, 47), E = c(-0.74, 3, 73),
    Q = c(-0.85, 0.2, 72), G = c(0.48, 0, 1),
    H = c(-0.4, -0.5, 82), I = c(1.38, -1.8, 57),
    L = c(1.06, -1.8, 57), K = c(-1.5, 3, 73),
    M = c(0.64, -1.3, 75), F = c(1.19, -2.5, 91),
    P = c(0.12, 0, 42), S = c(-0.18, 0.3, 31),
    T = c(-0.05, -0.4, 45), W = c(0.81, -3.4, 130),
    Y = c(0.26, -2.3, 107), V = c(1.08, -1.5, 43)
  )
  
  AAidx <- rbind(AAidx, tmp)
  
  if (!is.null(customprops)) AAidx <- rbind(AAidx, customprops)
  
  aaidx <- AAidx[, -1]
  row.names(aaidx) <- AAidx[, 1]
  
  n <- length(props)
  
  # Standardize H0 to H
  
  H0 <- as.matrix(aaidx[props, ])
  
  H <- matrix(ncol = 20, nrow = n)
  for (i in 1:n) {
    H[i, ] <-
      (H0[i, ] - mean(H0[i, ])) / (sqrt(sum((H0[i, ] - mean(H0[i, ]))^2) / 20))
  }
  AADict <- c(
    "A", "R", "N", "D", "C", "E", "Q", "G", "H", "I",
    "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"
  )
  dimnames(H) <- list(props, AADict)
  
  # Compute H^{1,2, ..., n}_{i, j}
  
  Theta <- vector("list", lambda)
  for (i in 1:lambda) Theta[[i]] <- vector("list", n)
  
  xSplitted <- strsplit(x, split = "")[[1]]
  
  N <- length(xSplitted)
  
  for (i in 1:lambda) {
    for (j in 1:n) {
      for (k in 1:(N - i)) {
        Theta[[i]][[j]][k] <-
          H[props[j], xSplitted[k]] * H[props[j], xSplitted[k + i]]
      }
    }
  }
  
  # Compute tau
  
  tau <- sapply(unlist(Theta, recursive = FALSE), mean)
  
  # Compute first 20 features
  
  fc <- summary(factor(xSplitted, levels = AADict), maxsum = 21)
  Pc1 <- fc / (1 + (w * sum(tau)))
  names(Pc1) <- paste("Pc1.", names(Pc1), sep = "")
  
  # Compute last n * lambda features
  
  Pc2 <- (w * tau) / (1 + (w * sum(tau)))
  names(Pc2) <- paste(
    "Pc2", as.vector(outer(props, 1:lambda, paste, sep = ".")),
    sep = "."
  )
  
  # Combine 20 + (n * lambda) features
  
  Pc <- c(Pc1, Pc2)
  
  return (Pc)
}


#' CTD Descriptors - Composition
#'
#' This function calculates the Composition descriptor of the
#' CTD descriptors (dim: 24).
#'
#' @param x A character vector, as the input protein sequence.
#'
#' @return A length 24 named vector
#'


extractCTDC_revised <- function(x) {
  
  AADict <- c(
    "A", "C", "D", "E", "F", "G", "H", "I",
    "L", "K", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y" 
  )
  
  if (all(strsplit(x, split = "")[[1]] %in% AADict) == FALSE) {
    stop("x has unrecognized amino acid type")
  }
  
  group1 <- list(
    "hydrophobicity" = c("R", "K", "E", "D", "Q", "N"),
    "normwaalsvolume" = c("G", "A", "S", "T", "P", "D", "C"),
    "polarity" = c("L", "I", "F", "W", "C", "M", "V", "Y"),
    "polarizability" = c("G", "A", "S", "D", "T"),
    "charge" = c("K", "R"),
    "secondarystruct" = c("E", "A", "L", "M", "Q", "K", "R", "H"),
    "solventaccess" = c("A", "L", "F", "C", "G", "I", "V", "W"),
    "surfacetension" = c("G", "Q", "D", "N", "A", "H", "R")
  )
  
  group2 <- list(
    "hydrophobicity" = c("G", "A", "S", "T", "P", "H", "Y"),
    "normwaalsvolume" = c("N", "V", "E", "Q", "I", "L"),
    "polarity" = c("P", "A", "T", "G", "S"),
    "polarizability" = c("C", "P", "N", "V", "E", "Q", "I", "L"),
    "charge" = c(
      "A", "N", "C", "Q", "G", "H", "I", "L",
      "M", "F", "P", "S", "T", "W", "Y", "V"
    ),
    "secondarystruct" = c("V", "I", "Y", "C", "W", "F", "T"),
    "solventaccess" = c("R", "K", "Q", "E", "N", "D"),
    "surfacetension" = c("K", "T", "S", "E", "C")
  )
  
  group3 <- list(
    "hydrophobicity" = c("C", "L", "V", "I", "M", "F", "W"),
    "normwaalsvolume" = c("M", "H", "K", "F", "R", "Y", "W"),
    "polarity" = c("H", "Q", "R", "K", "N", "E", "D"),
    "polarizability" = c("K", "M", "H", "F", "R", "Y", "W"),
    "charge" = c("D", "E"),
    "secondarystruct" = c("G", "N", "P", "S", "D"),
    "solventaccess" = c("M", "S", "P", "T", "H", "Y"),
    "surfacetension" = c("I", "L", "M", "F", "P", "W", "Y", "V")
  )
  
  
  xSplitted <- strsplit(x, split = "")[[1]]
  n <- nchar(x)
  
  # Get groups for each property & each amino acid
  
  g1 <- lapply(group1, function(g) length(which(xSplitted %in% g)))
  names(g1) <- paste(names(g1), "Group1", sep = ".")
  g2 <- lapply(group2, function(g) length(which(xSplitted %in% g)))
  names(g2) <- paste(names(g2), "Group2", sep = ".")
  g3 <- lapply(group3, function(g) length(which(xSplitted %in% g)))
  names(g3) <- paste(names(g3), "Group3", sep = ".")
  CTDC <- unlist(c(g1, g2, g3)) *100 / n
  # This reorders the groups from
  # p1.g1, p2.g1, p3.g1 ...
  # to
  # p1.g1, p1.g2, p1.g3 ...
  # This reordering is not really important
  # just to make the results look pretty
  ids <- unlist(lapply(1:8, function(x) x + c(0, 8, 16)))
  
  #  CTDC[ids]
  
  val<-c()
  for(i in 1:24){
    val<-c(val, CTDC[[i]])
    
  }
  # return (val)
  
  return (CTDC[ids])
}

#extractCTDC1(x)



#' CTD Descriptors - Transition
#'
#' This function calculates the Transition descriptor of the
#' CTD descriptors (dim: 24).
#'
#' @param x A character vector, as the input protein sequence.
#'
#' @return A length 24 named vector
#'



extractCTDT_revised <- function(x) {
  
  AADict <- c(
    "A", "C", "D", "E", "F", "G", "H", "I",
    "L", "K", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y" 
  )
  
  if (all(strsplit(x, split = "")[[1]] %in% AADict) == FALSE) {
    stop("x has unrecognized amino acid type")
  }
  
  group1 <- list(
    "hydrophobicity" = c("R", "K", "E", "D", "Q", "N"),
    "normwaalsvolume" = c("G", "A", "S", "T", "P", "D", "C"),
    "polarity" = c("L", "I", "F", "W", "C", "M", "V", "Y"),
    "polarizability" = c("G", "A", "S", "D", "T"),
    "charge" = c("K", "R"),
    "secondarystruct" = c("E", "A", "L", "M", "Q", "K", "R", "H"),
    "solventaccess" = c("A", "L", "F", "C", "G", "I", "V", "W"),
    "surfacetension" = c("G", "Q", "D", "N", "A", "H", "R")
  )
  
  group2 <- list(
    "hydrophobicity" = c("G", "A", "S", "T", "P", "H", "Y"),
    "normwaalsvolume" = c("N", "V", "E", "Q", "I", "L"),
    "polarity" = c("P", "A", "T", "G", "S"),
    "polarizability" = c("C", "P", "N", "V", "E", "Q", "I", "L"),
    "charge" = c(
      "A", "N", "C", "Q", "G", "H", "I", "L",
      "M", "F", "P", "S", "T", "W", "Y", "V"
    ),
    "secondarystruct" = c("V", "I", "Y", "C", "W", "F", "T"),
    "solventaccess" = c("R", "K", "Q", "E", "N", "D"),
    "surfacetension" = c("K", "T", "S", "E", "C")
  )
  
  group3 <- list(
    "hydrophobicity" = c("C", "L", "V", "I", "M", "F", "W"),
    "normwaalsvolume" = c("M", "H", "K", "F", "R", "Y", "W"),
    "polarity" = c("H", "Q", "R", "K", "N", "E", "D"),
    "polarizability" = c("K", "M", "H", "F", "R", "Y", "W"),
    "charge" = c("D", "E"),
    "secondarystruct" = c("G", "N", "P", "S", "D"),
    "solventaccess" = c("M", "S", "P", "T", "H", "Y"),
    "surfacetension" = c("I", "L", "M", "F", "P", "W", "Y", "V")
  )
  
  
  xSplitted <- strsplit(x, split = "")[[1]]
  n <- nchar(x)
  
  G <- vector("list", 8)
  for (i in 1:7) G[[i]] <- rep(NA, n)
  
  # Get groups for each property & each amino acid
  
  for (i in 1:8) {
    try(G[[i]][which(xSplitted %in% group1[[i]])] <- "G1")
    try(G[[i]][which(xSplitted %in% group2[[i]])] <- "G2")
    try(G[[i]][which(xSplitted %in% group3[[i]])] <- "G3")
  }
  
  # Combine single amino acids by a 2-length step
  
  for (i in 1:8) G[[i]] <- paste(G[[i]][-n], G[[i]][-1], sep = "")
  G <- lapply(G, function(x)
    factor(x, levels = c(
      "G1G2", "G2G1", "G1G3",
      "G3G1", "G2G3", "G3G2",
      "G1G1", "G2G2", "G3G3"
    )))
  
  GSummary <- lapply(G, summary)
  
  # Compute (n_rs + n_sr) / (N - 1)
  
  CTDT <- vector("list", 8)
  
  for (i in 1:8) {
    CTDT[[i]][1] <- sum(GSummary[[i]][c("G1G2", "G2G1")]) / (n - 1)
    #CTDT[[i]][1]<-CTDT[[i]][1]*100
    CTDT[[i]][2] <- sum(GSummary[[i]][c("G1G3", "G3G1")]) / (n - 1)
    #CTDT[[i]][2]<-CTDT[[i]][2]*100
    CTDT[[i]][3] <- sum(GSummary[[i]][c("G2G3", "G3G2")]) / (n - 1)
    #CTDT[[i]][3]<-CTDT[[i]][3]*100
  }
  
  CTDT <- unlist(CTDT)
  
  names(CTDT) <- paste(
    "prop", rep(1:8, each = 3), ".",
    rep(c("Tr1221", "Tr1331", "Tr2332"), times = 8),
    sep = ""
  )
  
  return(CTDT*100)
}

#extractCTDT1(x)


#' CTD Descriptors - Distribution
#'
#' This function calculates the Distribution descriptor of the
#' CTD descriptors (dim: 120).
#'
#' @param x A character vector, as the input protein sequence.
#'
#' @return A length 120 named vector
#'



extractCTDD_revised <- function(x) {
  
  AADict <- c(
    "A", "C", "D", "E", "F", "G", "H", "I",
    "L", "K", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y" 
  )
  
  if (all(strsplit(x, split = "")[[1]] %in% AADict) == FALSE) {
    stop("x has unrecognized amino acid type")
  }
  
  group1 <- list(
    "hydrophobicity" = c("R", "K", "E", "D", "Q", "N"),
    "normwaalsvolume" = c("G", "A", "S", "T", "P", "D", "C"),
    "polarity" = c("L", "I", "F", "W", "C", "M", "V", "Y"),
    "polarizability" = c("G", "A", "S", "D", "T"),
    "charge" = c("K", "R"),
    "secondarystruct" = c("E", "A", "L", "M", "Q", "K", "R", "H"),
    "solventaccess" = c("A", "L", "F", "C", "G", "I", "V", "W"),
    "surfacetension" = c("G", "Q", "D", "N", "A", "H", "R")
  )
  
  group2 <- list(
    "hydrophobicity" = c("G", "A", "S", "T", "P", "H", "Y"),
    "normwaalsvolume" = c("N", "V", "E", "Q", "I", "L"),
    "polarity" = c("P", "A", "T", "G", "S"),
    "polarizability" = c("C", "P", "N", "V", "E", "Q", "I", "L"),
    "charge" = c(
      "A", "N", "C", "Q", "G", "H", "I", "L",
      "M", "F", "P", "S", "T", "W", "Y", "V"
    ),
    "secondarystruct" = c("V", "I", "Y", "C", "W", "F", "T"),
    "solventaccess" = c("R", "K", "Q", "E", "N", "D"),
    "surfacetension" = c("K", "T", "S", "E", "C")
  )
  
  group3 <- list(
    "hydrophobicity" = c("C", "L", "V", "I", "M", "F", "W"),
    "normwaalsvolume" = c("M", "H", "K", "F", "R", "Y", "W"),
    "polarity" = c("H", "Q", "R", "K", "N", "E", "D"),
    "polarizability" = c("K", "M", "H", "F", "R", "Y", "W"),
    "charge" = c("D", "E"),
    "secondarystruct" = c("G", "N", "P", "S", "D"),
    "solventaccess" = c("M", "S", "P", "T", "H", "Y"),
    "surfacetension" = c("I", "L", "M", "F", "P", "W", "Y", "V")
  )
  
  
  xSplitted <- strsplit(x, split = "")[[1]]
  n <- nchar(x)
  
  G <- vector("list", 8)
  for (i in 1:8) G[[i]] <- rep(NA, n)
  
  # Get groups for each property & each amino acid
  
  for (i in 1:8) {
    try(G[[i]][which(xSplitted %in% group1[[i]])] <- "G1")
    try(G[[i]][which(xSplitted %in% group2[[i]])] <- "G2")
    try(G[[i]][which(xSplitted %in% group3[[i]])] <- "G3")
  }
  
  # Compute Distribution
  
  D <- vector("list", 8)
  for (i in 1:8) D[[i]] <- matrix(ncol = 5, nrow = 3)
  
  for (i in 1:8) {
    for (j in 1:3) {
      inds <- which(G[[i]] == paste0("G", j))
      quartiles <- floor(length(inds) * c(0.25, 0.5, 0.75))
      
      quartiles[which(quartiles <= 0)] <- 1
      
      D[[i]][j, ] <- if (length(inds) > 0) {
        (inds[c(1, quartiles, length(inds))]) * 100 / n
      } else {
        0
      }
    }
  }
  
  D <- do.call(rbind, D)
  D <- as.vector(t(D))
  
  names(D) <- paste(
    rep(paste("prop", 1:8, sep = ""), each = 15),
    rep(rep(c(".G1", ".G2", ".G3"), each = 5), times = 8),
    rep(paste(".residue", c("0", "25", "50", "75", "100"), sep = ""), times = 24),
    sep = ""
  )
  
  flag<-matrix(0,120,1)
  lc<-c()
  
  i<-1
  while(i<=15){
    #print(i)
    
    if(flag[i,1]==0){
      z<-i
      count<-1
      while(count<=3){
        # print(z)
        lc<-c(lc,D[[z]])
        flag[z,1]<-1
        z<-z+5
        count<-count+1
      }
    }
    i<-i+1
    #print(i)
  }
  
  # lc 
  
  i<-16
  while(i<=30){
    #print(i)
    
    if(flag[i,1]==0){
      z<-i
      count<-1
      while(count<=3){
        # print(z)
        lc<-c(lc,D[[z]])
        flag[z,1]<-1
        z<-z+5
        count<-count+1
      }
    }
    i<-i+1
    #print(i)
  }
  
  # lc
  
  i<-31
  while(i<=45){
    #print(i)
    
    if(flag[i,1]==0){
      z<-i
      count<-1
      while(count<=3){
        # print(z)
        lc<-c(lc,D[[z]])
        flag[z,1]<-1
        z<-z+5
        count<-count+1
      }
    }
    i<-i+1
    #print(i)
  }
  
  # lc
  
  i<-46
  while(i<=60){
    #print(i)
    
    if(flag[i,1]==0){
      z<-i
      count<-1
      while(count<=3){
        # print(z)
        lc<-c(lc,D[[z]])
        flag[z,1]<-1
        z<-z+5
        count<-count+1
      }
    }
    i<-i+1
    #print(i)
  }
  
  #lc
  
  
  i<-61
  while(i<=75){
    #print(i)
    
    if(flag[i,1]==0){
      z<-i
      count<-1
      while(count<=3){
        # print(z)
        lc<-c(lc,D[[z]])
        flag[z,1]<-1
        z<-z+5
        count<-count+1
      }
    }
    i<-i+1
    #print(i)
  }
  
  #lc
  
  
  i<-76
  while(i<=90){
    #print(i)
    
    if(flag[i,1]==0){
      z<-i
      count<-1
      while(count<=3){
        # print(z)
        lc<-c(lc,D[[z]])
        flag[z,1]<-1
        z<-z+5
        count<-count+1
      }
    }
    i<-i+1
    #print(i)
  }
  
  #lc
  
  
  
  i<-91
  while(i<=105){
    #print(i)
    
    if(flag[i,1]==0){
      z<-i
      count<-1
      while(count<=3){
        # print(z)
        lc<-c(lc,D[[z]])
        flag[z,1]<-1
        z<-z+5
        count<-count+1
      }
    }
    i<-i+1
    #print(i)
  }
  
  #lc
  
  
  
  
  i<-106
  while(i<=120){
    #print(i)
    
    if(flag[i,1]==0){
      z<-i
      count<-1
      while(count<=3){
        # print(z)
        lc<-c(lc,D[[z]])
        flag[z,1]<-1
        z<-z+5
        count<-count+1
      }
    }
    i<-i+1
    #print(i)
  }
  
  #lc
  
  return (lc)
}



#x<-"DLGPPISLERLDVGTNLGNAIAKLEAKELLESSD"

extractSSF <- function(x) {
  
  
  
  #ncrna_na<-read.csv("/Applications/kayes/Antiviral\ peptide\ work/for\ X\ sequences/VEEV/VEEV_4/VEEV_4.csv",header = TRUE)
  
  #l<-nrow(ncrna_na)
  #l
  ncrna_ss<-x
  ncrna<-AAStringSet(ncrna_ss)
  #Secondary structure features
  
  #BiocManager::install(c("DECIPHER", "BiocGenerics","parallel"))
  
  #library(RSQLite)
  
  
  
  
  
  
  
  

  
  
  #fas <- system.file("extdata", "50S_ribosomal_protein_L2.fas", package="DECIPHER")
  
  #file = "C:\\Soil ML works\\Training\\AMR seq\\1561417752.fas.1"
  #file = "C:\\Soil ML works\\Training\\PATRIC_Non-AMR seq\\1561502488.fas.1_deleted _one_4290_no_seq.1"
  #file = "C:\\Soil ML works\\cleaned soil_all.fasta"
  
  #aa = readAAStringSet(file)
  #length(aa)
  
  
  
  
  
  #aa[[1]]
  #dna <- readDNAStringSet(fas)
  #aa <- translate(dna)
  #aa
  hec <- PredictHEC(ncrna)
  hec
  ln<-length (hec)
  #t<-hec[1]
  #bb<-toString(t)
  #bb[1]
  # ss<-strsplit(t, "")[[1]]
  # len<-length(ss)
  
  # hec[1459]
  
  
  l<-length(ncrna)
  
  #col_n<-"ss_1"
  
  #for(k in 2:6){
  #  col_n<-paste0(col_n," ","ss_",k)
  
  #}
  
  #ls<-c(col_n)
  ls<-c()
  
  for(i in 1:l){
    t<-hec[i]
    data_psi<-strsplit(t, "")[[1]]
    len<-length(data_psi)
    add_ch<-""
    add_ch_exclude<-""
    flag<-""
    freq<-0
    x<-0
    y<-0
    z<-0
    SH<-0
    SE<-0
    SC<-0
    
    
    
    for(v in 1:len){
      if(toString(data_psi[v])=="H"){
        SH<-SH+v
        SH
        
      }
      if(toString(data_psi[v])=="E"){
        SE<-SE+v
        SE
        
      }
      if(toString(data_psi[v])=="C"){
        SC<-SC+v
        SC
        
      }
      add_ch<-paste0(add_ch,toString(data_psi[v]))
    }
    
    rr <- rle(strsplit(add_ch,"")[[1]])
    rr
    #SH<-sum(rr$lengths[which(rr$values == "H")])
    MH<-max(rr$lengths[which(rr$values == "H")])
    if (MH==-Inf || MH==Inf){
      MH<-0
    }
    
    CMVH<-SH/(len*(len-1))
    NMH<-MH/len
    
    #SE<-sum(rr$lengths[which(rr$values == "E")])
    ME<-max(rr$lengths[which(rr$values == "E")])
    ME
    if (ME==-Inf || ME==Inf){
      ME<-0
    }
    
    CMVE<-SE/(len*(len-1))
    NME<-ME/len
    
    #SC<-sum(rr$lengths[which(rr$values == "C")])
    MC<-max(rr$lengths[which(rr$values == "C")])
    
    CMVC<-SC/(len*(len-1))
    #NMC<-MC/len
    rr_len<-length(rr$values)
    rr_len
    for(v in 1:rr_len){
      if(rr$values[v]!="C")  {
        if(rr$values[v]!= flag){
          add_ch_exclude<-paste0(add_ch_exclude, toString(rr$values[v]))
          flag<-toString(rr$values[v])
          freq<-freq+1
        }
      }
    }
    add_ch_exclude
    freq
    
    #count_EHE<-str_count(add_ch_exclude, fixed("EHE"))
    count_EHE<-gregexpr("(?=EHE)",add_ch_exclude,perl=TRUE)
    count_EHE1<-count_EHE[[1]]
    count_EHE1
    if(count_EHE1<0){
      count_EHE_len<-0
    }
    else{
      count_EHE_len<-length(count_EHE1)
    }
    
    count_EHE_len
    if(freq>2){
      f_EHE<-count_EHE_len/(freq-2)
    }
    else{
      f_EHE<-0
    }
    
    f_EHE
    #add_line<-paste0(toString(CMVH), " ", toString(CMVE), " ", toString(CMVC), " ", toString(NMH), " ", toString(NME), " ", toString(f_EHE))
    #add_line
    
    #ls<-c(ls,add_line)
    ls<-c(ls,CMVH, CMVE, CMVC, NMH, NME, f_EHE)
    
    
  }
  
  return (ls)
}

#--------- assemble extraced features -------------#
feature_extraction <- function(fasta_file_path, aa_idx_path = "./data/AAidx.csv", raw_sequences=NA){
  if(is.na(raw_sequences)){
    lines = readLines(fasta_file_path)
    sequences = c()
    sequences_name = c()
    for(line in lines){
      if(startsWith(line, ">")){
        sequences_name <- append(sequences_name, trimws(gsub(">", "",line)))
      }
      else if(line != ""){
        sequences <- append(sequences, trimws(line))
      }
    }
  } else {
    sequences <- raw_sequences
  }
  
  
  
  all_generated_features = data.frame()
  feat_names = paste0("aac_", 1:20, sep = "")
  feat_names = append(feat_names, paste0("dipep_", 1:400, sep = ""))
  feat_names = append(feat_names, paste0("pseudo_", 1:25, sep = ""))
  feat_names = append(feat_names, paste0("amphipseudo_", 1:30, sep = ""))
  feat_names = append(feat_names, paste0("comp_", 1:24, sep = ""))
  feat_names = append(feat_names, paste0("tran_", 1:24, sep = ""))
  feat_names = append(feat_names, paste0("dist_", 1:120, sep = ""))
  feat_names = append(feat_names, paste0("ss_", 1:6, sep = ""))
  
  for (k in feat_names)  all_generated_features [[k]] <- as.numeric()
  
  seqs_length = c()
  
  for(i in 1:length(sequences)){
    current_seq = sequences[i]
    seqs_length <- append(seqs_length, length(current_seq))
    aa_comp = extractAAC_revised(current_seq)
    dc_comp = extractDC_revised(current_seq)
    pseaac_comp=extractPAAC_revised(current_seq, aa_idx_path)
    apseaac_comp=extractAPAAC_revised(current_seq, aa_idx_path)
    ctd_comp = extractCTDC_revised(current_seq)
    ctd_trans = extractCTDT_revised(current_seq)
    ctd_distr = extractCTDD_revised(current_seq)
    ss_struct = extractSSF(current_seq)
    all_generated_features[i,] <- c(aa_comp, dc_comp, pseaac_comp, apseaac_comp, ctd_comp, ctd_trans, ctd_distr, ss_struct)
  }
  colnames(all_generated_features) <- feat_names
  
  features_selected=c("aac_1","aac_2","aac_3","aac_4","aac_6","aac_7","aac_8","aac_9","aac_10","aac_11","aac_12","aac_15","aac_16","aac_17","aac_18","aac_19","aac_20","dipep_32","dipep_51","dipep_111","dipep_211","dipep_220","dipep_340","pseudo_1","pseudo_2","pseudo_3","pseudo_4","pseudo_5","pseudo_10","pseudo_11","pseudo_12","pseudo_14","pseudo_16","pseudo_18","pseudo_20","pseudo_21","pseudo_22","pseudo_23","pseudo_24","pseudo_25","amphipseudo_21","amphipseudo_22","amphipseudo_23","amphipseudo_24","amphipseudo_25","amphipseudo_26","amphipseudo_27","amphipseudo_29","amphipseudo_30","comp_1","comp_2","comp_3","comp_4","comp_5","comp_6","comp_10","comp_11","comp_13","comp_14","comp_15","comp_16","comp_17","comp_18","comp_19","comp_21","comp_22","comp_23","comp_24","tran_1","tran_2","tran_3","tran_4","tran_5","tran_6","tran_11","tran_12","tran_13","tran_14","tran_16","tran_17","tran_18","tran_19","tran_20","tran_21","tran_22","tran_23","tran_24","dist_1","dist_2","dist_3","dist_4","dist_7","dist_8","dist_9","dist_10","dist_11","dist_12","dist_13","dist_14","dist_15","dist_16","dist_17","dist_18","dist_22","dist_23","dist_24","dist_25","dist_26","dist_27","dist_28","dist_29","dist_30","dist_32","dist_34","dist_38","dist_41","dist_46","dist_47","dist_50","dist_52","dist_53","dist_55","dist_56","dist_61","dist_62","dist_63","dist_65","dist_67","dist_68","dist_70","dist_71","dist_72","dist_73","dist_76","dist_77","dist_78","dist_79","dist_82","dist_83","dist_84","dist_85","dist_86","dist_87","dist_88","dist_89","dist_90","dist_91","dist_93","dist_94","dist_97","dist_99","dist_100","dist_102","dist_103","dist_105","dist_106","dist_107","dist_108","dist_109","dist_112","dist_113","dist_114","dist_115","dist_116","dist_117","dist_118","dist_119","dist_120","ss_1")


  
  selected_features = all_generated_features[, which(colnames(all_generated_features) %in% features_selected)]
  return(selected_features)
}