# STEP 1: Protein database preparation

1. The following file was downloaded from: [http://gigaton.sigenae.org/ngspipelines/#!/NGSpipelines/Crassostrea%20gigas%20-%20GIGATON](http://gigaton.sigenae.org/ngspipelines/#!/NGSpipelines/Crassostrea%20gigas%20-%20GIGATON)

		contigs.fasta.transdecoder.pep	2017-02-10 11:30 	27M	

	for more details see Rhonda's jupyter notebook: [https://github.com/Ellior2/Fish-546-Bioinformatics/blob/master/analyses/DDA_2016/000-Remotelogin_filechange.ipynb](https://github.com/Ellior2/Fish-546-Bioinformatics/blob/master/analyses/DDA_2016/000-Remotelogin_filechange.ipynb)

2. The following common contaminant fasta files were obtained from the [CRAPome](https://www.nature.com/articles/nmeth.2557)

		contam.bovin.fa 	2017-02-10 11:30 	3.9K 
		contam.human.fa 	2017-02-10 11:30 	26K
		contam.other.fa 	2017-02-10 11:30 	8.6K	 

3. Steven cleaned up the file 'contigs.fasta.transdecoder.pep' and converted it to a fasta file

		Cg_Gigaton_proteins.fa 	2017-02-10 11:30 	21M	 
	
	Then combined it with the common contaminant fastas and moved it here: [http://owl.fish.washington.edu/halfshell/bu-git-repos/nb-2017/C_gigas/data/](http://owl.fish.washington.edu/halfshell/bu-git-repos/nb-2017/C_gigas/data/)
	
	See Steven's jupyter notebook for more details: [https://github.com/sr320/nb-2017/blob/master/C_gigas/00-Protein-database.ipynb](https://github.com/sr320/nb-2017/blob/master/C_gigas/00-Protein-database.ipynb)

# STEP 2: Converting .raw file to .mzXML files

#### NOTE: This step is likely not necessary because Comet apparently does this for you according to Emma. See [issue #471](https://github.com/sr320/LabDocs/issues/471)

This step was done by Rhonda with Sam's help. See Rhonda's jupyter notebook: [https://github.com/Ellior2/Fish-546-Bioinformatics/blob/master/analyses/DDA_2016/000-Remotelogin_filechange.ipynb](https://github.com/Ellior2/Fish-546-Bioinformatics/blob/master/analyses/DDA_2016/000-Remotelogin_filechange.ipynb) 

This is also part of the [DDA data analysis wiki](https://github.com/sr320/LabDocs/wiki/DDA-data-Analyses#convert-mass-spec-ms-raw-files-mzxml-files-for-use-in-comet)

	for file in *.raw
	    do
	    no_path=${file##*/}
	    no_ext=${no_path%.raw}
	    WINEPREFIX=~/.wine32 ReAdW.2016010.msfilereader.exe "$file" "$no_ext".mzXML
	    done

# STEP 3:  running Comet and TPP and Abacus  
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


# STEP 4:  Extracting protein abundance values from Abacus output 

Based on [the DDA analysis wiki](https://github.com/sr320/LabDocs/wiki/DDA-data-Analyses#convert-mass-spec-ms-raw-files-mzxml-files-for-use-in-comet) and [Emma's JOP geoduck paper](https://pubs.acs.org/doi/abs/10.1021/acs.jproteome.7b00288):

- For NMDS analysis of technical replicates, **ADJNSAF values** should be used
	- IF technical replicates 'cluster closely together and show less variability than biological replicates,' ADJNSAF values can be averaged for downstream NMDS+ANOSIM analysis to determine differences between sample types. Differentially abundant proteins between sample types can then be identified by Fisher's exact test. 
	- For comparing biological replicates to identify differentially abundant proteins, you can sum **spectral counts** from technical replicates and perform a univariate statistical test 
	
- For QSPEC analysis of two individual samples, **NUMSPECSTOT values** should be used

