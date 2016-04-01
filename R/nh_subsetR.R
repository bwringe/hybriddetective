#' @name nh_subset
#' @title Subset NewHybrids format datasets
#' @description Function allows the user to subset NewHybrids format datasets to explore the efficacy of variously sized panels
#' @param NHData The file path to a NewHybrids formatted dataset to be subsetted
#' @param loci A vector of Loci names to be removed from the dataset
#' @importFrom dplyr filter summarise ungroup group_by
#' @importFrom grid arrow unit
#' @importFrom stringr str_extract str_split
#' @export



nh_subset <- function(NHData, loci = NULL){


  NHD <- read.table(NHData, header = FALSE, quote = "", sep = "\t", stringsAsFactors = FALSE)

  NHD_nometa <- NHD[-(1:5),]

  temp <- as.data.frame(do.call(rbind, strsplit(NHD_nometa, " ")))

  colnamevec <- unlist(strsplit(NHD[5,], " "))

  colnames(temp) <- colnamevec

  if(length(loci) > 0){ ### May add ability to subset populations later - have if statement here to allow

    loci.to.remove <- which(colnamevec %in% loci)

    temp_subsetted <- temp[, -loci.to.remove]
      }

  insertNumIndivs <- paste("NumIndivs", nrow(temp_subsetted))

  insertNumLoci <- paste("NumLoci", length(temp_subsetted[-1])) ## will probably have to be -1

  ### hard coded stuff
  insertDigits <- "Digits 3"
  insertFormat <- "Format Lumped"

  insertLociNames <- paste(colnames(temp_subsetted), collapse = " ")

  data.out <- do.call(paste, c(data.frame(temp_subsetted[,]), sep = " "))

  data.out <- c(insertNumIndivs, insertNumLoci, insertDigits, insertFormat, insertLociNames, data.out)


  NHsplit <- c(stringr::str_split(string = NHData, pattern = "/"))

  outNameHold <- stringr::str_extract(NHsplit, paste0("[:word:]{3,}", ".txt"))
  dat.name.remove <- which(unlist(NHsplit) == outNameHold)
  outNameHold <- gsub(x = outNameHold, pattern = ".txt", replacement = "_subsetted.txt")
  outNameOut <- paste0(paste(unlist(NHsplit)[-dat.name.remove], collapse = "/"), "/", outNameHold)

  write(x = data.out, file = outNameOut)

  }