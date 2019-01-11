#STEP 1: Protein database preparation

1. The following file was downloaded from: [http://gigaton.sigenae.org/ngspipelines/#!/NGSpipelines/Crassostrea%20gigas%20-%20GIGATON](http://gigaton.sigenae.org/ngspipelines/#!/NGSpipelines/Crassostrea%20gigas%20-%20GIGATON)

		contigs.fasta.transdecoder.pep	2017-02-10 11:30 	27M	

	for more details see Rhonda's jupyter notebook: [https://github.com/Ellior2/Fish-546-Bioinformatics/blob/master/analyses/DDA_2016/000-Remotelogin_filechange.ipynb](https://github.com/Ellior2/Fish-546-Bioinformatics/blob/master/analyses/DDA_2016/000-Remotelogin_filechange.ipynb)

2. The following common contaminant fasta files were obtained from: (?????)

		contam.bovin.fa 	2017-02-10 11:30 	3.9K 
		contam.human.fa 	2017-02-10 11:30 	26K
		contam.other.fa 	2017-02-10 11:30 	8.6K	 

3. Steven cleaned up the file 'contigs.fasta.transdecoder.pep' and converted it to a fasta file

		Cg_Gigaton_proteins.fa 	2017-02-10 11:30 	21M	 
	
	Then combined it with the common contaminant fastas and moved it here: [http://owl.fish.washington.edu/halfshell/bu-git-repos/nb-2017/C_gigas/data/](http://owl.fish.washington.edu/halfshell/bu-git-repos/nb-2017/C_gigas/data/)
	
	See Steven's jupyter notebook for more details: [https://github.com/sr320/nb-2017/blob/master/C_gigas/00-Protein-database.ipynb](https://github.com/sr320/nb-2017/blob/master/C_gigas/00-Protein-database.ipynb)

#STEP 2: Converting .raw file to .mzXML files

#### NOTE: This step is likely not necessary because Comet apparently does this for you according to Emma. See [issue #471](https://github.com/sr320/LabDocs/issues/471)

This step was done by Rhonda with Sam's help. See Rhonda's jupyter notebook: [https://github.com/Ellior2/Fish-546-Bioinformatics/blob/master/analyses/DDA_2016/000-Remotelogin_filechange.ipynb](https://github.com/Ellior2/Fish-546-Bioinformatics/blob/master/analyses/DDA_2016/000-Remotelogin_filechange.ipynb) 

This is also part of the [DDA data analysis wiki](https://github.com/sr320/LabDocs/wiki/DDA-data-Analyses#convert-mass-spec-ms-raw-files-mzxml-files-for-use-in-comet)

	for file in *.raw
	    do
	    no_path=${file##*/}
	    no_ext=${no_path%.raw}
	    WINEPREFIX=~/.wine32 ReAdW.2016010.msfilereader.exe "$file" "$no_ext".mzXML
	    done

#STEP 3:  running Comet and TPP and Abacus  
This step was done by Steven on Feb 13, 2017. See Steven's notebook entry: [https://github.com/sr320/sr320.github.io/blob/master/_posts/2017-02-13-Going-through-DDA.md](https://github.com/sr320/sr320.github.io/blob/master/_posts/2017-02-13-Going-through-DDA.md) and jupyter notebooks: [03-DDA-RE-converted.ipynb](https://github.com/sr320/nb-2017/blob/master/C_gigas/03-DDA-RE-converted.ipynb), [03.5-DDA-pipeline%3F.ipynb](https://github.com/sr320/nb-2017/blob/master/C_gigas/03.5-DDA-pipeline%3F.ipynb), and [04-Exploring-Abacus-out.ipynb](https://github.com/sr320/nb-2017/blob/master/C_gigas/04-Exploring-Abacus-out.ipynb)

	/home/shared/comet/comet.2016012.linux.exe \
	-Pcomet.params.high-low_de \
	-DCg-Giga_cont_AA.fa \
	*.mzXML \
	&>> output.error.comet.log

	find *xml | \
	xargs basename -s .pep.xml | \
	xargs -I {} /usr/tpp_install/tpp/bin/xinteract \
	-dDECOY_ \
	-N{} \
	{}.pep.xml \
	-p0.9 \
	-OAp \
	&>> output.error.xin.log

	/usr/tpp_install/tpp/bin/ProteinProphet \
	interact*.pep.xml \
	interact-COMBINED.prot.xml \
	&>> output.error.PP.log

	java -Xmx16g -jar /home/shared/abacus/abacus.jar -p \
	Abacus.params

This step was also done by Sean on March 1, 2017. See [https://github.com/sr320/LabDocs/issues/492](https://github.com/sr320/LabDocs/issues/492)

Steven's output files are here:  [https://github.com/sr320/nb-2017/tree/master/C_gigas/data](https://github.com/sr320/nb-2017/tree/master/C_gigas/data) and [http://owl.fish.washington.edu/halfshell/working-directory/17-02-14b/](http://owl.fish.washington.edu/halfshell/working-directory/17-02-14b/)

Sean's output files are here:  [http://owl.fish.washington.edu/scaphapoda/Sean/Rhonda-2016-Oyster-Intermediates/](http://owl.fish.washington.edu/scaphapoda/Sean/Rhonda-2016-Oyster-Intermediates/)


#STEP 4:  Extracting protein abundance values from Abacus output 

Based on [the DDA analysis wiki](https://github.com/sr320/LabDocs/wiki/DDA-data-Analyses#convert-mass-spec-ms-raw-files-mzxml-files-for-use-in-comet) and [Emma's JOP geoduck paper](https://pubs.acs.org/doi/abs/10.1021/acs.jproteome.7b00288):

- For NMDS analysis of technical replicates, **ADJNSAF values** should be used
	- IF technical replicates 'cluster closely together and show less variability than biological replicates,' ADJNSAF values can be averaged for downstream NMDS+ANOSIM analysis to determine differences between sample types. Differentially abundant proteins between sample types can then be identified by Fisher's exact test. 
	- For comparing biological replicates to identify differentially abundant proteins, you can sum **spectral counts** from technical replicates and perform a univariate statistical test 
	
- For QSPEC analysis of two individual samples, **NUMSPECSTOT values** should be used

In trying to run NMDS on technical replicate ADJNSAF data, I found discrepencies between the ADJNSAF values in Steven's [ABACUS_output021417NSAF.tsv](https://github.com/sr320/nb-2017/blob/master/C_gigas/data/ABACUS_output021417NSAF.tsv) and Sean's [Abacus_output.tsv](http://owl.fish.washington.edu/scaphapoda/Sean/Rhonda-2016-Oyster-Intermediates/ABACUS_output.tsv).  I compared Steven's [ABACUS_output021417.tsv](https://github.com/sr320/nb-2017/blob/master/C_gigas/data/ABACUS_output021417.tsv) file (from which he made ABACUS_output021417NSAF.tsv, see his jupyter notebook [https://github.com/sr320/nb-2017/blob/master/C_gigas/04-Exploring-Abacus-out.ipynb](https://github.com/sr320/nb-2017/blob/master/C_gigas/04-Exploring-Abacus-out.ipynb)) with Sean's [Abacus_output.tsv](http://owl.fish.washington.edu/scaphapoda/Sean/Rhonda-2016-Oyster-Intermediates/ABACUS_output.tsv) and found no difference: 

####R code for comparing files
	
	install.packages("arsenal")
	library(arsenal)
	#Compare 02/14/2017 data with Sean's march 1 data
	data_SR <- read.csv("~/Documents/GitHub/OysterSeedProject/raw_data/ABACUS_output021417.tsv", sep = "\t" , header=TRUE, stringsAsFactors = FALSE)
	data_SB <- read.csv("~/Documents/GitHub/OysterSeedProject/raw_data/ABACUS_output.tsv", sep = "\t" , header=TRUE, stringsAsFactors = FALSE)
	compare(data_SR,data_SB)
	#Output:
	  	#Compare Object
	  	#Function Call: 
	  	# 	compare.data.frame(x = data_SR, y = data_SB)
	  	#Shared: 457 variables and 8443 observations.
	 	#Not shared: 0 variables and 0 observations.
	  	
	  	#Differences found in 0/456 variables compared.
	  	#0 variables compared have non-identical attributes.
	  
	###SHOWS NO DIFFERENCES BETWEEN FILES

####confirmed by command line diff command
 	#D-10-18-212-233:Desktop Shelly$ diff ~/Documents/GitHub/OysterSeedProject/raw_data/ABACUS_outputMar1.tsv ~/Documents/GitHub/OysterSeedProject/raw_data/ABACUS_output021417.tsv 
  	#D-10-18-212-233:Desktop Shelly$ 
  	
The values in Steven's ABACUS_output021417NSAF.tsv are in fact **NUMSPECSADJ values**

####R code to determine what the values in ABACUS_output021417NSAF.tsv are
	data_SR_NSAF <- read.csv("~/Documents/GitHub/OysterSeedProject/raw_data/ABACUS_output021417NSAF.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
	data_SB_NUMSPECADJ <- data_SB[,c(1,grep("NUMSPECSADJ", colnames(data_SB)))]
	colnames(data_SB_NUMSPECADJ) <- gsub("NUMSPECSADJ","ADJNSAF", colnames(data_SB_NUMSPECADJ))
	compare(data_SR_NSAF,data_SB_NUMSPECADJ)
	#Output:
	#Compare Object
	#Function Call: 
  	#	compare.data.frame(x = data_SR_NSAF, y = data_SB_NUMSPECADJ)
	#Shared: 46 variables and 8443 observations.
	#Not shared: 0 variables and 0 observations.

	#Differences found in 0/45 variables compared.
	#0 variables compared have non-identical attributes.

	###SHOWS NO DIFFERENCES BETWEEN FILES SO VALUES IN 
	###ABACUS_output021417NSAF.tsv ARE ACTUALLY
	###NUMSPECADJ VALUES!!!!

**Determined values in ABACUS_output021417NSAF.tsv are in fact NUMSPECADJ values**

**Determined values in Kaitlyn's file [ABACUSdata_only.csv](https://github.com/kaitlynrm/OysterSeedProject/blob/master/data/ABACUSdata_only.csv) are in fact the averages of technical replicate NUMSPECADJ values**  
- see [markdown summary](https://github.com/shellytrigg/OysterSeedProject/blob/master/analysis/nmds_R/CompareAbacusOutputFiles.md) of [R analysis](https://github.com/shellytrigg/OysterSeedProject/blob/master/analysis/nmds_R/CompareAbacusOutputFiles.Rmd)


##Next steps:
1. NMDS analysis
	- extract ADJNSAF values from [ABACUS_output021417.tsv](https://github.com/sr320/nb-2017/blob/master/C_gigas/data/ABACUS_output021417.tsv) 
	- Find appropriate data transformation/normalization if necessary
		- Emma log transformed her NSAF values before doing NMDS
	- try NMDS again
	- determine if replicates can be pooled

2. Try downstream analyses with NUMSPECSTOT values from [ABACUS_output021417.tsv](https://github.com/sr320/nb-2017/blob/master/C_gigas/data/ABACUS_output021417.tsv)
	- if it makes sense to sum NUMSPECSTOT values for replicates, try that and then try running stats on those values

3. Determine what values make sense to use in Hierarchical clustering analysis and ASCA, then re-do those analyses

4. Look more closely at development over time
	- Try a fold-change analysis with each developmental time point relative to day 0