#' @name freqbasedsim_GTFreq
#' @title Simulate Multi-Generational Hybrids - samples based on genotype frequency
#'
#' @description \code{freqbasedsim_GTFreq} generates simulated, centred Pure1, Pure2, F1, F2, BC1 and BC2 offspring based solely on the genotype frequencies of two ancestral populations provided. Allows the user to specify the number of individuals of each category to be simulated (including zero individuals should the user not wish to simulate a category).
#' @param NumSims an integer number of simulated datasets to be created. The default is 1
#' @param NumReps an integer number of replicates of each of the NumSims simulated dataset to be created. The default is 1
#' @param sample.sizePure an integer number of simulated Pure1 and Pure2 (centred ancestral populations) individuals to be created. The default is 200
#' @param sample.sizeF1 an integer number of simulated F1 individuals to be created. The default is 200
#' @param sample.sizeF2 an integer number of simulated F1 individuals to be created. The default is 200
#' @param sample.sizeBC an integer number of simulated Backross1 and Backcross2 (F1 X each of the Pure1 and Pure2) individuals to be created. The default is 200
#' @param outputName an optional character vector to be applied as the name of the output file(s). The default is NULL, in which case the output name is constructed from the name of the input file, with the suffix _SiRj_NH added. Where i is the number of simulations specified by NumSims, and j is the replicate number of the ith simulation, where j can take the values of 1:NumReps. NH refers to the fact that the output is in NewHybrids format
#' @param GenePopData file path to a GenePop formatted file containing genotypes from two (2) ancestral populations. This is the data from which the simulated hybrids will be constructed
#' @param pop.groups Optional character vector denoting how the individuals in the two ancestral populations should be named. The default is "PopA" and "PopB"
#' @export
#' @importFrom stringr str_extract str_extract_all str_split str_detect
#' @importFrom tidyr separate


freqbasedsim_GTFreq <- function(GenePopData, pop.groups = c("PopA", "PopB"), outputName = NULL, NumSims = 1, NumReps = 1, sample.sizePure = 200, sample.sizeF1 = 200, sample.sizeF2 = 200, sample.sizeBC = 200){


  max.ss <- max(sample.sizePure, sample.sizePure, sample.sizeF1, sample.sizeF2, sample.sizeBC, sample.sizeBC)

  if(max.ss<200){
    max.ss = 200
  }

  GenePop <- read.table(GenePopData, header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE)

  GPsplit <- c(stringr::str_split(string = GenePopData, pattern = "/"))

  outNameHold <- stringr::str_extract(GPsplit, paste0("[:word:]{3,}", ".txt"))
  outNameHold <- gsub(x = outputName, pattern = ".txt", replacement = "")

  NumIndivs <- (2*sample.sizePure) + sample.sizeF1 + sample.sizeF2 + (2*sample.sizeBC)

  stacks.version <- GenePop[1, ] # this could be blank or any other source. ## this was duplicated from another function - not sure if needed

  ## remove the first row which contains data normally ignored by GenePop, reformat data
  GenePop <- as.vector(GenePop)
  GenePop <- GenePop$V1[-1 ]
  GenePop <- data.frame(data=GenePop,ind=1:length(GenePop))
  GenePop$data <- as.character(GenePop$data)

  #ID the rows which flag the Populations
  Pops  <-  which(GenePop$data == "Pop" | GenePop$data == "pop" | GenePop$data == "POP")
  npops  <-  1:length(Pops)

  ## Seperate the data into the column headers (loci names) and the rest
  ColumnData <- GenePop[1:(Pops[1]-1),"data"]  ### SNP Names
  NumLoci <- length(ColumnData) ### NewHybrids Requires the number of LOCI be specified

  snpData <- GenePop[Pops[1]:NROW(GenePop),]  ### Genotypes - this is where the magic starts

  #Get a datafile with just the snp data no pops
  tempPops <- which(snpData$data == "Pop"| snpData$data == "pop" | snpData$data == "POP")
  snpData <- snpData[-tempPops, ]

  #Seperate the snpdata
  #First we pull out the population data which follows "TEXT ,  "
  temp <- tidyr::separate(snpData,data, into = c("Pops", "snps"), sep=",")
  temp$snps <- substring(temp$snps, 3) # delete the extra spaces at the beginning
  temp2 <- data.frame(do.call(rbind, stringr::str_extract_all(temp$snps, "[0-9]{3}")))

  ## Going to have to break the two alleles of the SNPS apart - this will thus double the number of columns
  ## SO <- will want to have SNP_A and SNP_A2
  ColumnData2 <- ColumnData ## Duplicatet the SNP names
  ColumnData2 <- paste(ColumnData2, "2", sep = ".")    ## add .2 to each duplicated name

  ## can't just append the duplicated names to the end of the original names - have to intersperse them
  places = rep(1:length(ColumnData)*2) ### creates a list of even numbers 2X as long as the number of columns i.e. the lenght of the original plus the duplicated names

  ## - this will also mark the position to insert the duplicates
  ColumnData.Dup = rep(NA, times = length(ColumnData)*2) ### make an object to feed names into
    for(i in 1:length(ColumnData)){ ### for loop to add original, then duplicated name
      a = places[i]-1 ## Original names go first, and are in the odd positions
      b = places[i] ### Duplicated names go second and are in the even posoitons
      Col.name.orig = ColumnData[i] ## Get the name in the ith position
      Col.name.plus2 = ColumnData2[i] ## get the name in the ith positon
      ColumnData.Dup[a] = Col.name.orig ## add the original name to the new vector
      ColumnData.Dup[b] = Col.name.plus2 ## add the duplicate name to the new vector
              } ## End of Loop


  #Contingency to see if R read in the top line as the "stacks version" -- modified to deal with the duplicated SNP names
  if (length(temp2) != length(ColumnData.Dup)){colnames(temp2) <- c(stacks.version, paste(stacks.version, "2", sep = "."), ColumnData.Dup)}
  if (length(temp2) == length(ColumnData.Dup)){colnames(temp2) <- ColumnData.Dup}

    ## Get the Alpha names
     NamePops=temp[,1] # Names of each

  if(length(pop.groups) == 0){ ### If unique grouping IDs ≠ number of "Pop" user must give vector of groupings equal to number of "Pop" or else the function will fail

    NameExtract=stringr::str_extract(NamePops, "[A-z]{3,}") ### if looking at higher order grouping (i.e. pops in  regions) can have more unique coding than "Pop" - will want to remove original names so can keep track of which unique groupings cross. i.e. Cross by "Pop", but remember ID of parents
          } ## End of IF statement


  # extract the text from the individuals names to denote population
  ## Now add the population tags using npops (number of populations and Pops for the inter differences)
    tPops <- c(Pops,NROW(GenePop))

    PopIDs <- NULL
      for (i in 2:length(tPops)){
        hold <- tPops[i]-tPops[i-1]-1

        if(i==length(tPops)){hold=hold+1}

        pophold <- rep(npops[i-1], hold)
        PopIDs <- c(PopIDs, pophold)
          } ## end of loop

    temp2$Pop <- PopIDs;


     if(length(pop.groups) != 0){
      hold.names=stringr::str_extract(NamePops, "[A-z]{3,}") ## This may need to be improved in published version
        for(i in 1:length(unique(PopIDs))){
          u.ID.no <- unique(PopIDs)[i]
          to <- min(which(PopIDs == u.ID.no))
          from <- max(which(PopIDs == u.ID.no))
          hold.names[to:from] = paste(pop.groups[i], hold.names[to:from], sep=".")
        } ## End I loop
      NameExtract <- hold.names
     } ## END IF

    ## get the nubmer of indivudals within each "Pop" grouping --- the Number of individuals in the two ancesntal populations need not be the same as the nubmer of individuals to be simulated
    PopLengths <- table(temp2$Pop)

    ## Need to be able to tell what row each individual is in, and what population it is
    ind.vector = c(1:nrow(temp)) ### make a vector that is the number of individuals
    ind.matrix = data.frame(temp2$Pop, ind.vector) ## add populatuions to that

    temp.split <- split(x = temp2, f = temp2$Pop)

    pop.recall <- NULL
      for(i in 1:length(temp.split)){
        popn <- paste(pop.groups[i], "pop", sep = "_")
        temp.split.hold = temp.split[[i]]
        temp.split.hold = temp.split.hold[-which(names(temp.split.hold) == "Pop")]
        assign(x = popn, value = temp.split.hold)
        pop.recall <- c(pop.recall, popn)
                } ## End of loop


    mat.name.recall <- NULL
      for(i in 1:length(pop.recall)){

        temp.mat <- data.frame(matrix(vector(), 2, length(temp2)/2))
        pop.get <- get(pop.recall[i])

          temp.mat.hold <- NULL
          for(k in 1:nrow(pop.get)){
            ind.hold <- pop.get[k, ]
            temp.mat[1,] <- t(t(ind.hold[c(T,F)]))
            temp.mat[2,] <-  t(t(ind.hold[c(F,T)]))
            temp.mat.hold <- rbind(temp.mat.hold, temp.mat)
                    } ## End of K loop

        mat.out.name <- paste(pop.recall[i],"matrix", sep = "_")
        assign(x = mat.out.name, value = temp.mat.hold)
        mat.name.recall <- c(mat.name.recall, mat.out.name)
                } ## End of I loop


                ##### Loooooooooooop tha Sims! #####

                for(sim in 1:NumSims){

                    ### MAKE PURE CROSS - centre the data -
                    pure.name.recall <- NULL
                      for(k in 1:length(pop.groups)){
                        pop1 <- get(mat.name.recall[k])
                        pop2 <- get(mat.name.recall[k])

                          off.interspersed.out <- NULL
                            for(i in 1:max.ss){

                              hold.off.pop1 <- apply(pop1, FUN = sample, 2, 1)
                              hold.off.pop2 <- apply(pop2, FUN = sample, 2, 1)
                              hold.off.interspersed <- data.frame(c(rbind(hold.off.pop1, hold.off.pop2)))

                              off.interspersed.out <- rbind(off.interspersed.out, t(hold.off.interspersed))
                                  } ## End of I loop

                          pure.name <- paste("Pure", pop.groups[k], sep = "_")
                          pure.name.recall <- c(pure.name.recall, pure.name)

                          assign(x = pure.name, value = off.interspersed.out)

                              } ## End of K loop


                    inv.pure.name.recall <- NULL
                    for(i in 1:length(pop.recall)){
                      temp.mat <- data.frame(matrix(vector(), 2, length(temp2)/2))
                      pop.get <- get(pure.name.recall[i])

                      temp.mat.hold <- NULL
                        for(k in 1:nrow(pop.get)){
                          ind.hold <- pop.get[k,]
                          temp.mat[1,] <- t(t(ind.hold[c(T,F)]))
                          temp.mat[2,] <-  t(t(ind.hold[c(F,T)]))
                          temp.mat.hold <- rbind(temp.mat.hold, temp.mat)
                            } ## End of K loop

                      inv.pure.out.name <- paste(pure.name.recall[i],"inv", sep = "_")
                      assign(x = inv.pure.out.name, value = temp.mat.hold)
                      inv.pure.name.recall <- c(inv.pure.name.recall, inv.pure.out.name)
                          } ## end of i loop

                    ### MAKE F1 CROSS ###

                    pop1 <- get(inv.pure.name.recall[1])
                    pop2 <- get(inv.pure.name.recall[2])

                    F1.out <- NULL
                    for(i in 1:max.ss){

                      hold.off.pop1 <- apply(pop1, FUN = sample, 2, 1)
                      hold.off.pop2 <- apply(pop2, FUN = sample, 2, 1)
                      hold.off.interspersed <- data.frame(c(rbind(hold.off.pop1, hold.off.pop2)))

                      F1.out <- rbind(F1.out, t(hold.off.interspersed))
                        } ## END of loop


                    temp.mat <- data.frame(matrix(vector(), 2, length(temp2)/2))
                    pop.get <- F1.out

                    inv.F1 <- NULL
                    for(k in 1:nrow(pop.get)){
                      ind.hold <- pop.get[k,]
                      temp.mat[1,] <- t(t(ind.hold[c(T,F)]))
                      temp.mat[2,] <-  t(t(ind.hold[c(F,T)]))
                      inv.F1 <- rbind(inv.F1, temp.mat)

                        } # end of k loop


                    ### MAKE F2 CROSS ###

                    pop1 <- inv.F1
                    pop2 <- inv.F1

                    F2.out <- NULL
                    for(i in 1:max.ss){

                      hold.off.pop1 <- apply(pop1, FUN = sample, 2, 1)
                      hold.off.pop2 <- apply(pop2, FUN = sample, 2, 1)
                      hold.off.interspersed <- data.frame(c(rbind(hold.off.pop1, hold.off.pop2)))

                      F2.out <- rbind(F2.out, t(hold.off.interspersed))
                      }


                    ### MAKE Back CROSS ###


                    BC.name.recall <- NULL
                    for(k in 1:length(pop.groups)){

                      pop1 <- get(inv.pure.name.recall[k])
                      pop2 <- inv.F1

                      off.interspersed.out <- NULL
                      for(i in 1:max.ss){

                        hold.off.pop1 <- apply(pop1, FUN = sample, 2, 1)
                        hold.off.pop2 <- apply(pop2, FUN = sample, 2, 1)
                        hold.off.interspersed <- data.frame(c(rbind(hold.off.pop1, hold.off.pop2)))

                        off.interspersed.out <- rbind(off.interspersed.out, t(hold.off.interspersed))

                          } # end of i loop

                      BC.name <- paste("BC", pop.groups[k], sep = "_")
                      BC.name.recall <- c(BC.name.recall, BC.name)

                      assign(x = BC.name, value = off.interspersed.out)

                        } ## end of k loop


                    for(i in 1:length(pure.name.recall)){
                      off.name <- paste(pure.name.recall[i], c(1:max.ss), sep="_")
                      hold.dat <- get(pure.name.recall[i])
                      hold.dat <- data.frame(off.name, hold.dat)
                      ColumnData.Dup.insert = c("ID", ColumnData.Dup)
                      colnames(hold.dat) = ColumnData.Dup.insert
                      assign(x = pure.name.recall[i], value = hold.dat)
                        } ## End I loop


                    f1.off.name <- paste("F1", c(1:max.ss), sep = "_")
                    F1.out <- data.frame(f1.off.name, F1.out)
                    colnames(F1.out) <- c("ID", ColumnData.Dup)

                    f2.off.name <- paste("F2", c(1:max.ss), sep = "_")
                    F2.out <-  data.frame(f2.off.name, F2.out)
                    colnames(F2.out) <- c("ID", ColumnData.Dup)


                    for(i in 1:length(BC.name.recall)){
                      off.name <- paste(BC.name.recall[i], c(1:max.ss), sep="_")
                      hold.dat <- get(BC.name.recall[i])
                      hold.dat <- data.frame(off.name, hold.dat)
                      ColumnData.Dup.insert = c("ID", ColumnData.Dup)
                      colnames(hold.dat) = ColumnData.Dup.insert
                      assign(x = BC.name.recall[i], value = hold.dat)
                        } ## End I loop


                    for(b in 1:length(pure.name.recall)){

                      fam.to.bind.name <- pure.name.recall[b]
                      fam.to.bind <- get(fam.to.bind.name)
                      indiv.hold <- fam.to.bind[,1]
                      loci.bind <- which(stringr::str_detect(string = colnames(fam.to.bind), pattern = "\\.2")==TRUE)

                      col.out <- NULL
                      for(k in 1:length(loci.bind)){
                        place.1 <- (loci.bind[k]-1)
                        place.2 <- loci.bind[k]
                        hold.col <- paste0(fam.to.bind[ ,place.1], fam.to.bind[ ,place.2])
                        col.out <- cbind(col.out, hold.col)
                        } ### End k loop

                      fam.reord <- cbind(indiv.hold, col.out)
                      colnames(fam.reord) <- c(colnames(fam.to.bind[1]), colnames(fam.to.bind[c((loci.bind-1))]))
                      assign(x = fam.to.bind.name, fam.reord)
                        } # End b loop

                    for(b in 1:length(pure.name.recall)){

                      fam.to.remove.untyped.name <- pure.name.recall[b]

                      fam.to.remove.untyped <- get(fam.to.remove.untyped.name)
                      fam.to.remove.untyped[which(stringr::str_detect(string = fam.to.remove.untyped, pattern = "000")==TRUE)] = "000000"
                      assign(x = fam.to.remove.untyped.name, value = fam.to.remove.untyped)
                    } # End b loop



                    ##### F1 #####

                    fam.to.bind.name <- "F1.out"

                    fam.to.bind <- get(fam.to.bind.name)
                    indiv.hold <- fam.to.bind[,1]
                    loci.bind <- which(stringr::str_detect(string = colnames(fam.to.bind), pattern = "\\.2")==TRUE)

                    col.out <- NULL
                    for(k in 1:length(loci.bind)){
                      place.1 <- (loci.bind[k]-1)
                      place.2 <- loci.bind[k]
                      hold.col <- paste0(fam.to.bind[,place.1], fam.to.bind[,place.2])
                      col.out <- cbind(col.out, hold.col)
                        } ## End K loop

                    fam.reord <- cbind(indiv.hold,col.out)
                    colnames(fam.reord) <- c(colnames(fam.to.bind[1]), colnames(fam.to.bind[c((loci.bind-1))]))
                    assign(x = fam.to.bind.name, fam.reord)


                    fam.to.remove.untyped.name <- "F1.out"

                    fam.to.remove.untyped <- get(fam.to.remove.untyped.name)
                    fam.to.remove.untyped[which(stringr::str_detect(string = fam.to.remove.untyped, pattern = "000")==TRUE)] = "000000"
                    assign(x = fam.to.remove.untyped.name, value = fam.to.remove.untyped)



                    fam.to.bind.name <- "F2.out"

                    fam.to.bind <- get(fam.to.bind.name)
                    indiv.hold <- fam.to.bind[,1]
                    loci.bind <- which(stringr::str_detect(string = colnames(fam.to.bind), pattern = "\\.2")==TRUE)

                    col.out <- NULL
                    for(k in 1:length(loci.bind)){
                      place.1 <- (loci.bind[k]-1)
                      place.2 <- loci.bind[k]
                      hold.col <- paste0(fam.to.bind[,place.1], fam.to.bind[,place.2])
                      col.out <- cbind(col.out, hold.col)

                    }

                    fam.reord <- cbind(indiv.hold,col.out)
                    colnames(fam.reord) <- c(colnames(fam.to.bind[1]), colnames(fam.to.bind[c((loci.bind-1))]))
                    assign(x = fam.to.bind.name, fam.reord)


                    fam.to.remove.untyped.name <- "F2.out"

                    fam.to.remove.untyped <- get(fam.to.remove.untyped.name)
                    fam.to.remove.untyped[which(stringr::str_detect(string = fam.to.remove.untyped, pattern = "000")==TRUE)] = "000000"
                    assign(x = fam.to.remove.untyped.name, value = fam.to.remove.untyped)


                    for(b in 1:length(BC.name.recall)){

                      fam.to.bind.name <- BC.name.recall[b]

                      fam.to.bind <- get(fam.to.bind.name)
                      indiv.hold <- fam.to.bind[,1]
                      loci.bind <- which(stringr::str_detect(string = colnames(fam.to.bind), pattern = "\\.2")==TRUE)

                      col.out <- NULL
                      for(s in 1:length(loci.bind)){
                        place.1 <- (loci.bind[s]-1)
                        place.2 <- loci.bind[s]
                        hold.col <- paste0(fam.to.bind[,place.1], fam.to.bind[,place.2])
                        col.out <- cbind(col.out, hold.col)
                          } ## End s loop

                      fam.reord <- cbind(indiv.hold,col.out)
                      colnames(fam.reord) <- c(colnames(fam.to.bind[1]), colnames(fam.to.bind[c((loci.bind-1))]))
                      assign(x = fam.to.bind.name, fam.reord)

                          } ## end b loop

                    for(b in 1:length(BC.name.recall)){

                      fam.to.remove.untyped.name <- BC.name.recall[b]

                      fam.to.remove.untyped <- get(fam.to.remove.untyped.name)
                      fam.to.remove.untyped[which(stringr::str_detect(string = fam.to.remove.untyped, pattern = "000")==TRUE)] = "000000"
                      assign(x = fam.to.remove.untyped.name, value = fam.to.remove.untyped)
                        } ## End b loop



                      #### Now recompile the NewHybrids #####


                      pop.names <- c(pure.name.recall, "F1.out", "F2.out", BC.name.recall)

                      no.sim.keep.vec <- c(sample.sizePure, sample.sizePure, sample.sizeF1, sample.sizeF2, sample.sizeBC, sample.sizeBC)

                      popvecout <- NULL
                      for(i in 1:length(pop.names)){

                        if(no.sim.keep.vec[i] > 0){
                          pvecmake <- paste0(pop.names[i], "_", c(1:no.sim.keep.vec[i]))
                          popvecout <- c(popvecout, pvecmake)
                          } # End IF statement
                        } # End loop



                      sim.out <- NULL
                      for(i in 1:length(pop.names)){
                        if(no.sim.keep.vec[i] > 0){
                          hold.pop <- get(pop.names[i])
                          hold.pop <- hold.pop[1:no.sim.keep.vec[i],] ### added to vary sample sizes
                          sim.out <- rbind(sim.out, hold.pop)
                        } # End IF statement
                      } # End loop


                      sim.out[ ,1] <- c(1:nrow(sim.out))
                      Loci.sim <- do.call(paste, c(data.frame(sim.out[,]), sep = " "))
                      cd2 <- paste(ColumnData, collapse = " ")

                      insertNumIndivs <- paste("NumIndivs", NumIndivs)
                      insertNumLoci <- paste("NumLoci", NumLoci)
                      insertYourDigits <- "Digits 3"
                      insertFormat <- "Format Lumped"
                      insertLociName <- paste("LocusNames", cd2)

                      Loci.out <- c(insertNumIndivs, insertNumLoci, insertYourDigits, insertFormat, insertLociName,   Loci.sim)

                      outNameGive <- gsub(x = GenePopData, pattern = ".txt", replacement = "")
                      outNameGive <- paste0(outNameGive, "_S", sim)

                      popvecout.fname <- gsub(x = GenePopData, pattern = ".txt", replacement = "_individuals.txt")
                      write(x = popvecout, file = popvecout.fname)

                      for(r in 1:NumReps){

                        outNameGiveOut <- paste0(outNameGive, "R", r, "_NH.txt")
                        write.table(x = Loci.out, file = outNameGiveOut, row.names = FALSE, col.names = FALSE, quote = FALSE)

                            } ## end r loop

                } ### End SIM Loop

} ## End function
