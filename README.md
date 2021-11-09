# hybriddetective

***
### **hybriddetective** is an R package designed to streamline, and where possible automate, the detection of hybrids using **NewHybrids** (Anderson and Thompson 2002) (available at <https://github.com/eriqande/newhybrids>) by moving the entire process into the R environment. This includes the development and testing of diagnostic panels of markers including measures of their accuracy and efficiency and the creation of multi-generation simulated hybrid datasets.
***


### Preamble  
**hybriddetective** is designed to work in conjunction with, and to link the functions contained in the R packages **parallelnewhybrid** (<https://github.com/bwringe/parallelnewhybrid>) and **genepopedit** (<https://github.com/rystanley/genepopedit>), to bring the entire hybrid detection workflow into the R environment.  

Please ensure that you have the correct version of **NewHybrids** installed. The source code and instructions for installation can be found at (<https://github.com/eriqande/newhybrids>)

By using **hybriddetective**, together with the packages [parallelnewhybrid](https://github.com/bwringe/parallelnewhybrid) and [genepopedit](https://github.com/rystanley/genepopedit) the user is able to complete the entire hybrid detection analysis starting with genotypic/genomic data through to identification of hybrid individuals, with quantifiable levels of certainty of assignment. These packages facilitate the manipulation of genotype files, the development of panels of markers not in linkage disequilibrium that have the highest loci-specific Fst, the simulation of multi-generation hybrids and subsequent quantification of the efficacy and accuracy (Vaha and Primmer 2006) of panels to identify hybrids, and the determination of hybrid category of experimental data.  
    
In addition, **hybriddetective** includes plotting functionality to generate both exploratory and publishable figures.  
    
**hybriddetective** and **parrallelnewhybrid** were designed primarily to work with the Bayesian model-based genetic hybrid detection program [NewHybrids](https://github.com/eriqande/newhybrids) (Anderson and Thompson 2002). However by using the format interconversion tools provided in **hybriddetective** and **genepopedit**, users are able to analyze simulated and experimental datasets in other programs such as *STRUCTURE* (Pritchard et al. 2000).  



***

### Overview of functions:  

**Function name** | **Synopsis** | **Main use**  
------------|------------------------------------------|------------------------------------------------------------------
**getTopLoc.R**  | Creates a panel comprised of the top *n* loci not in linkage disequilbrium. Exports a list of loci, and a dataset to be used in the simulation functions in [parallelnewhybrid](https://github.com/bwringe/parallelnewhybrid), or another program.| [_Data Preparation_](#gettoploc)
**freqbasedsim_GTFreq.R**  | Creates simulated hybrid datasets comprised of *two pure parental groups*, *F1*, *F2*, and *backcrosses to both ancestral groups*. Samples based solely on allele frequencies in the provided populations. Alleles are sampled without replacement from the two populations. Allows the number of individuals of each genotype frequency category output to be specified by the user.| [_Data Preparation_](#freqbasedADV)
**freqbasedsim_AlleleSample.R**  | Creates simulated hybrid datasets comprised of *two pure parental groups*, *F1*, *F2*, and *backcrosses to both ancestral groups*. Samples by taking one allele at each locus from each of the two parental populations in each generation. Alleles are sampled without replacement from the two populations.. Also allows the number of individuals of each genotype frequency category output to be specified by the user.| [_Data Preparation_](#freqbasedUB2)
**nh_analysis_generateR.R** |   Merges specified simulated reference datasets with the genotypes of unknown/experimental individuals, producing a file to be analyzed by [NewHybrids](https://github.com/eriqande/newhybrids) | [_Data Preparation_](#nhdatageneratR)
**nh_analysis_simulateR_generateR.R** |   Quickly creates simulated reference populations from user provided data, and merges these with the genotypes of unknown/experimental individuals, producing a file to be analyzed by *NewHybrids*| [_Data Preparation_](#nhanalysisdata)
**nh_subsetR.R** |  Subset loci from *NewHybrids* format datasets| [_Data Preparation_](#nhsubset)
**nh_analysis_data_generatoR.R** |   Quickly creates simulated reference populations from user provided data, and merges these with the genotypes of unknown/experimental individuals, producing a file to be analyzed by *NewHybrids*| [_Data Preparation_](#nhanalysisdata)
**nh_Zcore.R** |  Inserts known genotype category assignments to *NewHybrids* format files (see Anderson 2002)| [_Data Preparation_](#zcore)
**nh_preCheckR.R** | Checks all *NewHybrids* results within a given directory to ensure they have converged| [_Error Checking and Diagnostics_](#precheck)
**nh_plotR.R** |  Plots cumulative probability of assignments from *NewHybrids* results for each individual | [_Error Checking and Diagnostics_](#plotR)
**nh_multiplotR.R** |  Plots cumulative probability of assignments from *NewHybrids* results for each individual for all results within a given folder | [_Error Checking and Diagnostics_](#multiplotR)
**nh_accuracy_checkR.R** | Calculates and can print to screen the proportion of simulated individuals assigned to the correct hybrid category at three levels of stringency (PofZ >= 0.5, 0.75 and 0.9) | [_Error Checking and Diagnostics_](#accuracy)
**hybridPowerComp.R** |  Evaluates the accuracy with which *NewHybrids* assigns individuals of known hybrid class to the correct hybrid class in simulated datasets at varying levels of stringency (critical posterior probability thresholds). Will also write graphical and numerical results to a directory provided by the user. Also allows the efficacy of multiple panels of varying numbers of markers to be compared.| [_Quantification and Analysis_](#power)
**nh_panel_delta_plotR.R** |  Plots the change in individual genotype class probability of assignment for different panel sizes.| [_Quantification and Analysis_](#deltaplot)
**nh_build_example_results.R** |  Installs pre-analyzed example data to be used familiarize oneself with the functions *nh_preCheckR*, *nh_multiplotR*, and *hybridPowerComp* | [_Example Data_](#example)

***  
### Package Installation
**hybriddetective** requires the packages **genepopedit** (<https://github.com/rystanley/genepopedit>) and **parallelnewhybrid** (<https://github.com/bwringe/parallelnewhybrid>), and these must be installed from GitHub prior to the installation of **hybriddetective**.  

Installation and use of [NewHybrids](https://github.com/eriqande/newhybrids) to conduct analyses in parallel is also described in detail in the README for [parallelnewhybrid](https://github.com/bwringe/parallelnewhybrid) and the associated publication [Wringe et al. 2017](http://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12597/full)

**genepopedit**, **parallelnewhybrid**, and **hybriddetective** can be installed from GitHub using the R package *devtools* and the function *install_github*: 

```
devtools::install_github("bwringe/parallelnewhybrid") #required
devtools::install_github("rystanley/genepopedit") #required
devtools::install_github("bwringe/hybriddetective") #This package
```

**NOTE:** **hybriddetective** also relies on functions from the packages *dplyr*, *ggplot2*, *grid*, *parallel*, *plyr*, *stringr*, *scales*, *reshape2* and *tidyr*. The user should ensure these installed from CRAN prior to installing **hybriddetective**.  

**NOTE:** certain functions in **hybriddetective** rely on the programs **PGDSpider version 2.1.1.0** and **PLINK version 1.90 beta**. These can be downloaded from <http://www.cmpg.unibe.ch/software/PGDSpider/> and <https://www.cog-genomics.org/plink2> respectively, and must be installed prior to running **hybriddetective**.

***

***  
### Example Data:
Example datasets have been provided as R images (.rda files). These can be loaded into your workspace using the ``data()`` function  

**Example dataset** | **Contents**  
------------|---------------------------------------------------------------  
*PurePops* | A *GENEPOP* format file containing SNP genotypes of 150 individuals in each of two populations at 100 loci. To use this example data it must first be saved with the extension ".txt" to an empty folder on your hard drive.
*PurePops_Zeds* | A data frame of known genotype category assignments (Zeds) and the indviduals to which to match them in the *NewHybrids* format files created in the examples. Refer to Anderson 2002 and the function *nh_zcore* for more information.

```r

### PREPARATION OF EXAMPLE DATA

## Get the file path to the working directory, will be used to allow a universal example
path.hold <- getwd()

## Get the genotype file included along with the hybriddetective package and make it an object
sim_PurePops <- hybriddetective::PurePops

## Get the Zeds data file included along with the hybriddetective package and make it an object
sim_Zeds <- hybriddetective::PurePops_Zeds

## Save the genotype data to the working directory as a file called "SimPurePops.txt"
write.table(x = sim_PurePops, file = paste0(path.hold, "/SimPurePops.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)

## Save the Zeds data to the working directory as a file called "SimZeds.txt"
write.table(x = sim_Zeds, file = paste0(path.hold, "/SimZeds.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE)

## Create an empty folder within the working directory. This will be the directory for all examples.
dir.create(paste0(path.hold, "/hybriddetective example"))

## Copy the genotype file to the new folder
file.copy(from = paste0(path.hold, "/SimPurePops.txt"), to = paste0(path.hold, "/hybriddetective example"))

## Copy the Zeds data file to the new folder
file.copy(from = paste0(path.hold, "/SimZeds.txt"), to = paste0(path.hold, "/hybriddtective example"))


## Clean up the working directory by deleting the two files
file.remove(paste0(path.hold, "/SimPurePops.txt"))

file.remove(paste0(path.hold, "/SimZeds.txt"))

## NOTE: You must recall the file path to the hybriddetective example folder you created to use in the examples below!!!

your.path = paste0(path.hold, "/hybriddetective example/")

```


<a name = "example"/>
<h4 class="text-primary"> nh_build_example_results.R</h4>
This function will write and delete the provided pre-analyzed data that can be used to explore the use and results of the functions *nh_preCheckR*, *nh_multiplotR*, and *hybridPowerComp*. This function also demonstrates the manner/format in which [parallelnewhybrid](https://github.com/bwringe/parallelnewhybrid) exports the results of *NewHybrids* analyses. The function will write the example data to a folder/file path specified by the user.

```r
## INSTALL EXAMPLE DATA

## Get the file path to the working directory, will be used to allow a universal example
path.hold <- getwd()

## Note: the path provided must end in '/' 

nh_build_example_results(dir = paste0(path.hold, "/")) # this will write the example data to the working directory

## To remove the EXAMPLE DATA

nh_build_example_results(dir = paste0(path.hold, "/"), remove_example = TRUE)


```


### **hybriddetective** Contributors

* Dr. Brendan Wringe <https://github.com/bwringe>; bwringe(at)gmail.com; Corresponding author  
* Dr. Ian Bradbury <http://tinyurl.com/h3qrhkv>  
* Dr. Nick Jeffery <https://github.com/NickJeff13>  
* Dr. Ryan Stanley <https://github.com/rystanley> 

This package has been designed to work in conjunction with the packages [parallelnewhybrid](https://github.com/bwringe/parallelnewhybrid) [(Wringe et al. 2017)](http://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12597/abstract) and [genepopedit](https://github.com/rystanley/genepopedit) [(Stanley et al. 2017)](http://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12569/abstract) to offer a more efficient and repeatable method of detecting hybrid individuals based on genomic data. This is accomplished by moving the workflow into the R environment. The package is open, and interested individuals are encouraged to modify the code for their purposes if needed. I will endeavour to respond in a timely manner to any and all inquiries, and will try to add requested functions, or functionality where requested.

* If you don't understand anything, please contact me (bwringe[at]gmail.com)
* Any ideas on how to improve the functionality are very much appreciated.
* If you spot a typo, feel free to edit and send a pull request.  

Pull request how-to:

* Click the edit this page on the sidebar.
* Make the changes using GitHub’s in-page editor and save.
* Submit a pull request and include a brief description of your changes. (e.g. "spelling errors" or "indexing error").  

To cite the current version of **hybriddetective**, please refer to our zenodo DOI: 


*** 
***  

### Function descriptions

<a name = "gettoploc"/>
<h4 class="text-primary"> getTopLoc.R</h4>
This function is designed to produce panels of the *n* loci not in linkage disequilibrium with the highest loci-specific Fst, where *n* is a number specified by the user. The data provided to the function are broken into two data sets, one to be used to calculate the loci-specific Fsts, and the other to be used to created simulated hybrid datasets to test the efficacy of the panel. This is done to prevent high-grading bias where the same data is used to develop the panel(s) that is then used in their validation (Anderson 2010). The function will save three files to the folder that contains the data provided. The first file is a list of the top *n* loci, ordered by Fst, which can be used to subset data in with the function **genepop_subset** from the **genepopedit** package. The second file is a list of the individuals from each of the two populations that were subsetted out to form the testing panel. The final file are the genotypes of the individuals listed in the second file, at the top *n* loci not in linkage disequilibrium (listed in the first file). This dataset can be used in the **freqbasedsim** functions to create simulated hybrids for panel efficacy testing. 
**NOTE:** Data must be entered into the function as a two-population GENEPOP format file.  
**NOTE:** The programs **PLINK** and **PGDspider** are required to run **getTopLoc.R**. Refer to package installation for downloading instructions.  

**Parameter** | **Description**
------------|---------------------------------------------------------------  
*GPD*| A file path to the GENEPOP format file you wish to create your panel from
*LDpop*| A string which denotes which of the two populations you wish to calculate linkage disequilibrium in. The options are "Pop1" or "Pop2", or "Both" if the LD is to be calculated based on both populations.
*panel.size*| An integer number of loci to include in the panel
*FST.threshold*| The minimum Fst threshold required to retain a locus. This value is used where there are a large number of relatively uninformative (i.e. low Fst) loci to improve the speed of the function.
*r2.threshold*| The minimum r^2 threshold to consider a pair of loci to be in linkage disequilibrium.
*ld.window*| Number of adjacent SNPs to compare each SNP against for LD - default is NULL, which translates to a window size of 99999, which should allow for all pairwise comparisons.
*allocate.PGD.RAM*| An integer value in GB to specify the maximum amount of RAM to allocate to PGDspider. The default is 1 GB, which should be sufficient for most analyses. Note: This option cannot be changed on Windows (refer to PGDSpider help file for more information)
*where.PLINK*| A file path to the PLINK installation folder
*where.PGDspider*| A file path to the PGDspider installation folder


#### Usage

To create a panel consisting of the top 50 Loci from the example dataset *PurePops* (Note: There are only 100 loci total in the example data, so the requested panel size may not exceed 100)
```r
## There is no need with the example dataset to alter any function parameters apart from providing the required file paths and the desired panel size

## recall the filepath to the hybriddetective example folder
path.hold <- getwd()
your_example_data_path = paste0(path.hold, "/hybriddetective example/SimPurePops.txt")

### assume PGDSpider and PLINK are both installed in a folder called "software"
your_plink_path = "~/ ... /software/PLINK Folder/"
your_pgdspider_path = "~/ ... /software/PGDSpider Folder/"

getTopLoc(GPD = your_example_data_path, panel.size = 50, where.PLINK = your_plink_path, where.PGDspider = your_pgdspider_path)

### This will create a file containing the 50 loci with the greatest Fst values, that are not in LD above an r^2 of 0.2. 

### Note: The file "SimPurePops_50_Loci_Panel.txt" is the file to be used for the simulations

```

***

<a name = "freqbasedADV"/>
<h4 class="text-primary"> freqbasedsim_GTFreq.R</h4>
This function will simulate multigenerational hybrids based on the allele frequencies of the two reference populations provided. By so doing, the genotypes of the reference populations are centred, and linkage disequilibria or deviations from Hardy-Weinberg that may have been present in the sampled individuals are eliminated. **freqbasedsim_Advanced** allows users flexibility in specifying the number of individuals in each of the hybrid genotype frequency classes to be simulated.  

**Parameter** | **Description**
------------|---------------------------------------------------------------  
*GenePopData*| GenePop formatted file containing genotypes from two (2) ancestral populations. This is the data from which the simulated hybrids will be constructed  
*pop.groups*|An optional character vector denoting how the two ancestral populations should be named; default is PopA and PopB
*outputName*|An optioanal character vector to be applied as the name of the output. The default is NULL, in which case the output name is constructed from the name of the input, with the suffix *\_SiRj_NH* added where *i* is the number of *S*imulations corresponding to the output, and *j* is the number of *R*eplicates of the *ith* simulation. NH refers to the fact that the output is in *NewHybrids* format.  
*sample.sizePure*|An optional integer to specify the number of simulated  Pure1 and Pure2 individuals (centred ancestral populations) to be created. Default is 200  
*sample.sizeF1*|An optional integer to specify the number of simulated  F1 individuals (Pure1 X Pure2) to be created. Default is 200  
*sample.sizeF2*|An optional integer to specify the number of simulated  F2 dividuals (F1 X F1) to be created. Default is 200  
*sample.sizeBC*|An optional integer to specify the number of simulated  Back Cross1 and BackCross2 individuals (Pure1 X F1 and Pure2XF2) to be created. Default is 200  
*NumSims*|An integer number of simulated datasets to be created. Default is 1  
*NumReps*|An integer number of replicates to be created for each of the *i* simulated datasets specified in NumSims. Default is 1  

#### Usage
```r
## If the number of Pure, F1, F2 and Backcross individuals to be returned is not specified, this function will return 200 of each class. Because NewHybrids utilizes a Bayesian algorithm, results may be influenced by the proportions of individuals of different genotypic frequency categories provided, thus this function can be used to tailor the numbers of simulated individuals in each category to those expected of your study system.

### Suppose you expect that the majority of individuals in will be pure, with the numbers of Pure > F1 > F2 >= BC

### NOTE: "YourReferenceData.txt" is either data in that you wish to simulate in GenePop format OR
  ## if you going through the package example, it will be the file "SimPurePops_50_Loci_Panel.txt" created by getTopLoc

freqbasedsim_GTFreq(GenePopData = "~/ ... /YourReferenceData.txt", sample.sizePure = 100, sample.sizeF1 = 50, sample.sizeF2 = 10, sample.sizeBC = 10, NumSims = 3, NumReps = 3)

## This will create three independent simulations, each replicated three times.
```

***

<a name = "freqbasedUB2"/>
<h4 class="text-primary"> freqbasedsim_AlleleSample.R</h4>
This function will simulate multigenerational hybrids by drawing 2 alleles from each loci, without replacement from the two reference populations. The function requires that two reference populations must be provided as a single *GenePop* format text file, such as the file created by **getTopLoc**.  
The function first creates Pure1 and Pure2 individuals by drawing 2 alleles from each loci, without replacement from the two reference populations. By so doing, the genotypes of the reference populations are centred, and linkage disequilibria or deviations from Hardy-Weinberg that may have been present in the sampled individuals are eliminated. Subsequently, each F1 individual is created by drawing 1 allele per loci withouth replacement from the Pure1 individuals and 1 loci without replacement from the Pure2 individuals. F2 indivduals are created by drawning 2 alleles per loci without replacement from the F1 indivdiuals. Backcross individuals are created by drawing 1 allele per loci withouth replacement from the F1 individuals, and one without replacement from the appropriate Pure population individuals.

**Parameter** | **Description**
------------|---------------------------------------------------------------  
*GenePopData* | GenePop formatted file containing genotypes from two (2) ancestral populations. This is the data from which the simulated hybrids will be constructed  
*pop.groups* | An optional character vector denoting how the two ancestral populations should be named. Default is PopA and PopB  
*prop.sample* | The proportion of individuals in both ancestral PopA and PopB to sample to create the simulated hybrids (Pure1, Pure2, F1, F2, BC1, and BC2)  
*outputName* | An optional character vector to be applied as the name of the output. The default is NULL, in which case the output name is constructed from the name of the input, with the suffix *\_SiRj_NH* added where *i* is the number of *S*imulations corresponding to the output, and *j* is the number of *R*eplicates of the *ith* simulation. NH refers to the fact that the output is in *NewHybrids* format.  
*sample.sizePure1* | An optional integer to specify the number of simulated Pure1 individuals (centred ancestral PopA) to be output when the desired number is less than the number of individuals in Ancestral Population 1 * prop.sample. The default is NULL, where the number output = number of individuals in Ancestral Population 1 * prop.sample. If a number greater than number of individuals in Ancestral Population 1 * prop.sample is requested, the number of individuals in Ancestral Population 1 * prop.sample are output.  
*sample.sizePure2* | An optional integer to specify the number of simulated Pure2 individuals (centred ancestral PopB) to be output when the desired number is less than the number of individuals in Ancestral Population 2 * prop.sample. The default is NULL, where the number output = number of individuals in Ancestral Population 2 * prop.sample. If a number greater than number of individuals in Ancestral Population 2 * prop.sample is requested, the number of individuals in Ancestral Population 2 * prop.sample are output.  
*sample.sizeF1* | An optional integer to specify the number of F1 individuals to be output when the desired number of simulated F1 individuals is less than (number of individuals in Ancestral PopA + number of individuals in Ancestral PopB) * prop.sample. The default is NULL where the number returned = (number of individuals in Ancestral PopA + number of individuals in Ancestral PopB) * prop.sample. Unless sample.sizeF1 is explicitly stated, even when sample.sizePure1 and sample.sizePure2 are specified, the number of simulated F1 individuals returned will be equal to (number of individuals in Ancestral PopA + number of individuals in Ancestral PopB) * prop.sample.  
*sample.sizeF2* | An optional integer to specify the number of simulated F2 individuals to be output when the desired number of simulated F2 individuals is less than (number of individuals in Ancestral PopA + number of individuals in Ancestral PopB) * prop.sample. The default is NULL where the number returned = (number of individuals in Ancestral PopA + number of individuals in Ancestral PopB) * prop.sample. Unless sample.sizeF2 is explicitly stated, even when sample.sizePure1 and sample.sizePure2 are specified, the number of simulated F2 individuals returned will be equal to (number of individuals in Ancestral PopA + number of individuals in Ancestral PopB) * prop.sample.  
*sample.sizeBC1* | An optional integer to specify the number of simulated BC1 (PopA X F1) individuals to be output when the desired number of simulated BC1 individuals is less than the number of individuals in Ancestral PopA * prop.sample. The default is NULL where the number returned = number of individuals in Ancestral PopA * prop.sample. Unless sample.sizeBC1 is explicitly stated, even when sample.sizePure1 and sample.sizeF1 are specified, the number of simulated BC1 individuals returned will be equal to number of individuals in Ancestral PopA * prop.sample.  
*sample.sizeBC2* | An optional integer to specify the number of simulated BC2 (PopB X F1) individuals to be output when the desired number of simulated BC2 individuals is less than the number of individuals in Ancestral PopB * prop.sample. The default is NULL where the number returned = number of individuals in Ancestral PopB * prop.sample. Unless sample.sizeBC2 is explicitly stated, even when sample.sizePure2 and sample.sizeF1 are specified, the number of simulated BC2 individuals returned will be equal to number of individuals in Ancestral PopB * prop.sample.  
*NumSims* | An integer number of simulated datasets to be created. Default is 1  
*NumReps* | An integer number of replicates to be created for each of the *i* simulated datasets specified in NumSims. Default is 1  

#### Usage
```r
## This function will simulate multigenerational hybrids by drawing 2 alleles from each loci, without replacement from the two reference populations. The function requires that two reference populations must be provided as a single GenePop format text file, such as the file created by the function getTopLoc.  

### NOTE: "YourReferenceData.txt" is either data in that you wish to simulate in GenePop format OR if you going through the package example, it will be the file "SimPurePops_50_Loci_Panel.txt" created by getTopLoc

freqbasedsim_AlleleSample(GPD = "~/ ... /YourReferenceData.txt")
```

***

<a name = "nhdatageneratR"/>
<h4 class="text-primary"> nh_analysis_simulateR_generateR.R</h4>
This function quickly creates files to be used to detect hybridization within experimental data. It does so by creating simulated reference populations from user provided data, and then merging the simulated data with the genotypes of unknown/experimental individuals. The user is able to specify which hybrid categories to include in the simulated data, as well as the number of individuals per category to simulate. The function will save a *NewHybrids* format file for analysis, along with a dataframe containing the names of the individuals included (including those that were simulated) in the *NewHybrids* formatted file.  
**NOTE:** *nh_analysis_data_generatoR.R* creates a unique simulated data set each time it is run.   

**Parameter** | **Description**
------------|---------------------------------------------------------------  
*ReferencePopsData*| A file path to a GenePop formatted file containing genotypes from two (2) ancestral populations. This is the data from which the simulated hybrids will be constructed. 
*UnknownIndivs*| A dataframe of genotypes of unknown/experimental individuals to be analyzed for hybrid category generated using the function *genepop_flatten* from the package *genepopedit*. Individuals or loci can be subsetted from this flattened data set if desired, or manipulations can conducted before flattening using the appropriate functions from *genepopedit*. **NOTE:**: the number of loci in ReferencePopsData must equal the number in UnkownIndivs 
*outputName*| An optional character vector to be applied as the name of the output. The default is NULL, in which case the output name is constructed from the name of the input, with the *suffix _SiRj_NH* added where *i* is the number of simulations corresponding to the output, and *j* is the number of replicates of the *ith* simulation. NH refers to the fact that the output is in NewHybrids format 
*pop.groups*| An optional character vector denoting how the two ancestral populations should be named; default is "PopA" and "PopB" 
*sample.size*| An integer number of simulated individuals to be created for each of the six hybrid classes (e.g. 200 * each of Pure1, Pure2, F1, F2, BC1 and BC2 = 1200 total simulated individuals); default is 200 
*NumSims*| An integer number of simulated datasets to be created; default is 1 
*NumReps*| An integer number of replicates to be created for each of the *n* simulated datasets specified in NumSims; default is 1 
*cats.include*| Optional character vector list denoting which hybrid categories should be included in the output; default is to include all categories.

#### Usage

```r
### This function will create simulated reference data using the algorithm from "freqbasedsim_GTFreq" based on the genotypes of two reference populations provided by the user, and then combines these with unknown/experimental data

### NOTE: "TwoReferencePops.txt" data to be simulated and "GenotypesOfUnknownInds.txt" are genotypes of unknown/experimental individuals at the same loci
  ## if you going through the package example,"SimPurePops_50_Loci_Panel.txt" can be used for both "GenotypesOfUnknownInds.txt" and "TwoReferencePops.txt"

nh_analysis_simulateR_generateR(ReferencePopsData = "~/ ... /TwoReferencePops.txt", UnknownIndvs = "~/ ... /GenotypesOfUnknownInds.txt", sample.size = 50, NumSims = 3, NumReps = 3) ### create reference populations comprised of 50 individuals. Make 3 independent simulations, and 3 replicates of each.
```

***

<a name = "nhanalysisdata"/>
<h4 class="text-primary"> nh_analysis_generateR.R</h4>
This function combines pre-simulated datasets with experimental/unknown data. Both the simulated and experimental datasets can be in *GENEPOP* or *NewHybrids* formats. This function also allows the user to specify which hybrid categories to include from the simulated dataset. The function will save a *NewHybrids* format file for analysis, along with a dataframe containing the names of all individuals included in the *NewHybrids* formatted file.

**Parameter** | **Description**
------------|---------------------------------------------------------------  
*ReferencePopsData*| A file path to a either a *NewHybrids* or *GENEPOP* formatted file containing genotypes from the simulated ancestral populations. This can be the result of any of the *freqbasedsim* functions in this package, or a file created using the function *genepop_newhybrids* from the package *genepopedit*.  
*UnknownIndivs*| A file path to a file containing the genotypes of the individuals to be analyzed for possible hybrid ancestry. This can either be a *GENEPOP* or *NewHybrids* format file. **NOTE:**: the number of loci in ReferencePopsData must equal the number in UnkownIndivs 
*sim.pops.include*| An optional character vector list denoting which hybrid categories from the simulated data should be included in the output; default is Pure Population 1 and Pure Population 2. 
*outputName*| An optional character vector to be applied as the name of the output. The default is NULL, in which case the output name is constructed from the name of the input, with the *suffix _SiRj_NH* added where *i* is the number of simulations corresponding to the output, and *j* is the number of replicates of the *ith* simulation. NH refers to the fact that the output is in NewHybrids format 

#### Usage

```r
### This function will combine pre-created simulated reference with unknown/experimental data. 

## Note "PreSimulatedData.txt" is the results of any of the multi-generational simulation functions, and "GenotypesOfUnknownInds.txt" are genotypes of unknown/experimental individuals at the same loci
  ## if you going through the package example,"SimPurePops_50_Loci_Panel.txt" can be used for "GenotypesOfUnknownInds.txt" and the results of any of the multi-generational simulation functions can be used for "PreSimulatedData.txt"

nh_analysis_generateR(ReferencePopsData = "~/ ... /PreSimulatedData.txt", UnknownIndivs = "~/ ... /GenotypesOfUnknownInds.txt")
```

***

<a name = "nhsubset"/>
<h4 class="text-primary"> nh_subsetR.R</h4>
This function allows the user to subset markers out of *NewHybrids* formatted data files to test the efficacy of different panel sizes.  

**Parameter** | **Description**
------------|---------------------------------------------------------------  
*NHData*| The file path to the NewHybrids formatted dataset to be subsetted
*loci*| A vector of loci names to be retained in the dataset

#### Usage

```r
### Use to remove genotypes at specified loci to test panels of different sizes, or composition. Loci must be specified as the names of the loci.

### Find out what the names of the loci in your dataset are - assuming NewHybrids formatted file. 

loci.names <- scan("~/ ... /yourdata.txt", skip = 4, nlines = 1, what = "character")[-1] ### the first 4 lines are meta-data, read only the fifth line, remove the first entry which is "LocusNames"

### To keep the first 50 loci

nh_subsetR(NHData = "~/ ... /yourdata.txt", loci = loci.names[1:50]) ## function will save the new dataset in the same folder with the filename "yourdata_subsetted.txt"
```

***

<a name = "zcore"/>
<h4 class="text-primary"> nh_Zcore.R</h4>
This function will insert a column of known genotype categories <- Zvec into a NewHybrids formatted file to provide NewHybrids with prior information on the known genotype frequencies in the hybrid genotype frequency categories.The function matches a vector of categories (Zvec) to individuals, then output a NH format file, with a name indicating that the category vector has been added. The user must provide the Zvec(s) to the function in the form of a two column .csv file. The vector of genotype categories NEED NOT be the same length as the NewHybrids file - thus, time can be saved. 
NOTE Zvec file must have column names, and the first column must be named "Individual" and contain the names of individuals to which categories must be given that match those in the first column of the NewHybrids dataset. NOTE both multiapllyZvec and applyuniqueZvec cannot be used at the same time. 

**Parameter** | **Description**
------------|---------------------------------------------------------------  
*GetstheZdir*| A file path to the folder in which the *NewHybrids* formatted files to be appended reside
*multiapplyZvec*| **IF** *a single Zvec is to be applied to each **NewHybrids** format file*, use *multiapplyZvec*. Zvec will added in batch fashion to every *NewHybrids* format file within a working folder. *multiapplyZvec* must be specified as a file path + file name (e.g. "~/HoldenUniversity/CoolEthansComputer/HairDolls/Zvexxx.csv")
*applyuniqueZvec*| **IF** *each **NewHybrids** format file is to be given a UNIQUE vector of Zvecs* use *applyuniqueZvec*. **NOTE:** : to apply UNIQUE Zvecs, the Zvec files must all be placed in a single folder separate from the *NewHybrids* files *AND* they MUST follow the file name convention "NHFileName_Zvec.csv" where NHFileName is the same name as the file to which the Zvec is to be applied. Consequently, the number of *NewHybrids* files = number of Zvec files. *applyuniqueZvec* must be specified as a file path to the folder in which the Zvec files reside (e.g. "~/HarrisonUniversity/DEANGordonPritchard/LettersOfRecommendation/")  

#### Usage
```r
## to apply the same file of Zed scores to multiple files:

nh_zcore(GetstheZdir = "~/ ... /yourNHFiles/", multiapplyZvec = "~/ ... /yourZeds.txt")  ## NOTE: Single file

## to apply unique files of Zed scores to corresponding files:

nh_Zcore(GetstheZdir = "~/ ... /yourNHFiles/", applyuniqueZvec = "~/ ... /yourZedFiles/") ## Zed files must all be placed in a single folder separate from the NewHybrids files AND they MUST follow the file name convention "NHFileName_Zvec.csv" where NHFileName is the same name as the file to which the Zvec is to be applied
```

***

<a name = "precheck"/>
<h4 class="text-primary"> nh_preCheckR.R</h4>
When using *NewHybrids*, the Markov Chain Monte Carlo (MCMC) chains will occasionally fail to converge resulting, in *NewHybrids* spuriously returning the highest probability of assignment of nearly every individual to the F2 category. In such cases, the analysis should be re-run. *preCheckR.R* quickly assesses the results of many independent *NewHybrids* runs, and will flag instances in which the probability of assignment values indicate the chains appear to have failed to converge. Results flagged by *preCheckR.R* should be manually verified by the user.  
**NOTE:** several of the functions in **hybriddetective** will not function correctly, or will fail if they are provided with results of non-converged runs. Thus, it is strongly recommended that *preCheckR.R* be run on every *NewHybrids* result.  

**Parameter** | **Description**
------------|---------------------------------------------------------------  
*PreDir*| The directory in which the NewHybrids results (in individual folders as returned by parallelNH_XX) are located.
*propCutOff*| The proportion of individuals in either Pure Population 1 OR Pure Population 2 allowed to exceed the F2 assignment threshold (PofZCutOff). The default is 0.5. 
*PofZCutOff*| The threshold posterior probability of assignment (PofZ) to F2 above which Pure Population 1 or Pure Population 2 individuals are flagged to indicate possible non-convergence. The default is 0.1. 

#### Usage

```r
## precheckR will check for NewHybrids results in all folders within the folder specified by "PreDir". i.e. all results in the folders within the "NH.Results" folder created by parallelnh_xx from the package parallelnewhybrid

## to check results created by running parallelnh_xx

nh_preCheckR(PreDir = "~/ ... /NH.Results/", propCutOff = 0.5, PofZCutOff = 0.1) ## run with the default cut-off values; these could have been left blank.
### preCheckR will now flag any results where more than 50% (propCutOFF) of EITHER Pure1 or Pure2 individuals are given a probability of being an F2 greater than 10% (PofZCutOff)

```

***

<a name = "plotR"/>
<h4 class="text-primary"> nh_plotR.R</h4>
This function will plot the cumulative probability of assignment (PofZ) for each individual in a *NewHybrids* analysis. *nh_plotR.R* can be used with both simulated and experimental datasets. NEWHYBRIDS denotes one of the pure populations “pop1” and the other “pop2” (and thus backcrosses) at random (i.e. “pop1” can be equivalent to “pop2” in different runs). However, the ordering of individuals or samples remains constant and therefore, when examining the results of nh_multiplotR, it is important to bear in mind that the colouration of individuals in each plot can switch, yet the interpretation of the assignment(s) is the same. Because of this, *nh_plotR* contains an optional argument account for this where reversals have been noted.   

**Parameter** | **Description**
------------|---------------------------------------------------------------  
*NHResults*| A file path to the NewHybrids (PofZ file) result to be plotted 
*ReversePure* | A logical which can be used to reverse the Pure1/BC1, Pure2/BC2 colours if NEWHYBRIDS switches which population is called Population1 and Population2 between analyses of the same data (allows you to make your plots the same colours). Options are "1" or "2"
*ColourVector* | An optional vector to specify the six colours to be plotted as Pure1, Pure2, F1, F2, BC1, and BC2 respectively. The defaults are R's "red", "blue", "grey", "green", "black", and "yellow".

####Usage
```r
## Use plotR to visualize the cumulative probability of assignment for for each individual in the PofZ file specified by "NHResults"

nh_plotR(NHResults = "~/ ... /NH.Results/YourData_Results/YourData_PofZ.txt") ## plot results will be displayed by R 

```

***

<a name = "multiplotR"/>
<h4 class="text-primary"> nh_multiplotR.R</h4>
This function will plot the cumulative probability of assignment (PofZ) for each individual in a *NewHybrids* analysis, for all analyses within a given folder. *nh_plotR.R* can be used with both simulated and experimental datasets. This function is useful for visualizing results where **preCheckR.R** has indicated that *NewHybrids* may have failed to converge.  

**Parameter** | **Description**
------------|---------------------------------------------------------------  
*NHResults*| A file path to the directory in which the NewHybrids results (in individual folders as returned by parallelNH_XX) are located. 
*ColourVector* | An optional vector to specify the six colours to be plotted as Pure1, Pure2, F1, F2, BC1, and BC2 respectively. The defaults are R's "red", "blue", "grey", "green", "black", and "yellow".

#### Usage
```r
## Use multplotR to visualize the cumulative probability of assignment for each individual for multiple NewHybrids analyses. The function will plot the results in all folders within the directory specified by "NHResults". i.e. all results in the folders within the "NH.Results" folder created by parallelnh_xx

nh_multiplotR(NHResults = "~/ ... /NH.Results/") ## the function will create a folder called "NewHybrids Plots" within the directory specified by "NHResults", and the plots will be saved to that folder. 
```

***

<a name = "accuracy"/>
<h4 class="text-primary"> nh_accuracy_checkR.R</h4>
This function will quickly calculate, the proportion of individuals in a simulated dataset which were assigned to the correct hybrid category by *NewHybrids* at three probability of assignment thresholds (PofZ >= 0.5, 0.75 and 0.9). The proportion of individuals for each hybrid category correctly assigned will be returned as a dataframe with the name "NH.accuracy", and can also be printed to the screen if desired (print.results = TRUE). The user also has the option to check the proportion of individuals in all hybrid categories correctly identified as being of hybrid ancestry (i.e. proportion of F1, F2, BC1, and BC2 individuals given the greatest PofZ in an any of those categories as opposed to Pure1 or Pure2). *nh_accuracy_checkR.R* is intended primarily for exploratory analysis, and it's functionality is expanded upon in the function *hybridpowercomp*.  

**Parameter** | **Description**
------------|---------------------------------------------------------------  
*NHResultsDir*| A file path to the NewHybrids (PofZ file) result to be checked. 
*print.results*| A logical query for whether the function output should be printed to the screen in addition to exported as an object. Default is TRUE 
*all.hyb*| A logical query if the proportion of all individuals known to be hybrids were assigned to a hybrid category regardless if the category was the correct one. Default is FALSE


#### Usage

```r
## nh_accuracy_checkR can be used to quickly check the results of a NewHybrids analysis of simulated data to see the proportion of individuals of known genotype frequency class which have been assigned to the correct class at three different levels of PofZ stringency (i.e. 0.5, 0.75 and 0.9). 

nh_accuracy_checkR(NHResultsDir = "~/ ... /NH.Results/YourData_Results/", print.results = TRUE, all.hyb = FALSE) ## the function will return the results to the global environment as an object called "NH.accuracy"

```

***

<a name = "power"/>
<h4 class="text-primary"> hybridpowercomp</h4>
This function will extract data from *NewHybrids* runs and evaluate the correct assignment to a given hybrid classification (i.e. Pure1, Pure2, F1, F2, BC1, and BC2). From this data, hybrid power will generate summary plots which detail assignment success (among simulations and replicates) as a function of the classification threshold (minimum probability of assignment threshold required for an individual to be assigned a hybrid class, scaling from 0.50-1.00). Diagnostic plots are also provided which detail variation among simulations, as well as to which hybrid classes individuals are _missassigned_ to. Summary data and *ggplot* objects are saved to a .RData output automatically. Outputs from this function are saved to the folder "Figures and Data" which is generated within the directory (specified by *dir*) that contains the *NewHybrids* output folders.  Where panels of differrent sizes (number of markers) are tested, the function will partition the results according to the number of loci used in the assignment analysis. The plots will compare assignment power among different numbers of loci used.  
**NOTE:** If multiple panel sizes are not tested (i.e. only one panel size tested) *hybridpower* will return warnings related to missing factor levels from *ggplot2*. These warnings are expected, and will not affect the output or interpretation of the results.

**Parameter** | **Description**
------------|---------------------------------------------------------------  
*dir*| File path to the directory in which the NewHybrids results (in individual folders as returned by parallelNH_XX) are located. 
*filetag*| An optional character vector to be applied as the name of the outputs.
*Threshold*| A threshold which will be added to the plots showing the assignment success for different levels of probability of a given class estimated by NewHybrids. Default is (NULL) so if nothing is specified it will not add this to the output plots (success ~ threshold by class) 
*Thresholds*| A vector of thresholds which will be added to the plots showing the assignment success for different levels of probability of a given class estimated by NewHybrids. Default is (NULL) so if nothing is specified it will not add this to the output plots (success ~ threshold by class). 
*samplesize*| The number of individuals per NewHybrids class. By (default: NULL) this data will be extracted from the "*individuals.txt" output from parallelnewhybrids if present in the same folder as the PofZ file. This can also explicitly defined as a vector (6 values corresponding to # in P1,P2,F1,F2,BC1,BC2) or a path to a *_Individuals.txt. 
*CT*| The threshold posterior probability of assignment (PofZ) to F2 above which Pure Population 1 or Pure Population 2 individuals are flagged to indicate possible non-convergence. The default is 0.1.  
*CTI*| The proportion of individuals in either Pure Population 1 OR Pure Population 2 allowed to exceed the F2 assignment threshold (PofZCutOff). The default is 0.5.  

#### Usage

```r
## Use hybridpower to evaluate assignment of known individuals to the correct genotype frequency class. Function will create a folder called "Figures and Data" within the directory specified by "dir" and save results to this folder.
## NOTE: Unless an "_individuals" file from the "freqbasedsim" family of functions in the parallelnewhybrid package is present in the same folder as the results, the sample size must be specified.


hybridpowercomp(dir = "~/ ... /NH.Results/") ## Defaults appropriate for analysis
```

***


<a name = "deltaplot"/>
<h4 class="text-primary"> nh_panel_delta_plotR.R</h4>
This function is designed to illustrate the changes in individual probability of assignment to each of the six genotype frequency classes when analyzed using panels of differing numbers of markers. The user can specify if they would like the plot to be saved, or returned as a **ggplot2** object for further formatting.  

**Parameter** | **Description**
------------|---------------------------------------------------------------  
*GPD*| File path to the directory which holds the output from different runs through New Hybrids (e.g. 3 simulations with 3 replicate runs each through NH) note that this directory should only hold the output folders. 
*return.workspace*| A logical query of whether the ggplot2 object should be returned to the workspace for further editing or formatting. Can take values of "TRUE" or "FALSE". NOTE: if return.workspace = TRUE, the user must assign a variable to the function (i.e. mydata <- nh\_ind_panel_delta(...)) 
*save.plot*| A logical query of whether the plot should be saved as a file in the working directory. NOTE: The file will be saved "delta_plot.xxx", and so must be renamed if the function is to be run more than once to prevent overwriting. 
*plot.filetype*| If save.plot = TRUE, plot.filetype specifies the file type to save the plot as. The default is "png", with the option of "pdf". 

#### Usage

```r
### Use nh_panel_delta_plotR to visualize the change (if any) in individual assignment to genotype frequency class using panels comprised of different numbers of markers. The plot can returned to the global environment as an object called "delta.plot" by default, and can optionally be saved to the directory specified by GPD.

nh_panel_delta_plotR(GPD = "~/ ... /NH.Results/", return.workspace = TRUE, save.plot = FALSE) ### return the plot to the workspace, but do not save
```


***

### Citations

Anderson, E.C., and Thompson, E.A. 2002. A model-based method for identifying species hybrids using multilocus data. Genetics. 160(3): 1217-1229.

Anderson, E.C. 2010. Assessing the power of informative subsets of loci for population assignment: standard methods are upwardly biased. Molecular Ecology Resources. 10(4): 701-710.  

Pritchard, J.K., Stephens, M., and Donnelly, P. 2000. Inference of population structure using multilocus genotype data. Genetics 155(2): 945-959.

Stanley, R.R.E., Jeffery, N.W., Wringe, B.F., DiBacco, C., and Bradbury, I.R. 2017. genepopedit: a simple and flexible tool for manipulating multilocus molecular data in R. Molecular Ecology Resources. 17(1): 12-18.

Vaha, J.P. and Primmer, C.R. 2006. Efficiency of model-based Bayesian methods for detecting hybrid individuals under different hybridization scenarios and with different numbers of loci. Molecular Ecology. 15(1): 63-72.

Wringe B.F., Stanley, R.R.E., Jefferey, N.W., Anderson, E.C., and Bradbury, I.R. 2017. parallelnewhybrid: an R package for the parallelization of hybrid detection using newhybrids. Molecular Ecology Resources. 17(1): 91-95.

