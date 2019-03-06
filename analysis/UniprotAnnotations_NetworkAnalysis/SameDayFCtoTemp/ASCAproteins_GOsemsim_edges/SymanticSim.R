###make a get_sim_grid to get protein relationships
###all proteins will be connected to each other 
###edges can be colored white if the similary is low and black if it's high
###edge weights can be small if similarity is low and heavier if it's high


#load libraries
library(reshape2)

###############################
####ATTEMPT with slim terms##
################################
ASCA_GOslims <- unique(par_ch_anc_slim[-grep("GO:0008150",par_ch_anc_slim$term),2])

beach <- list()
for(i in 1:length(ASCA_GOslims)){
  temp_beach <- try(go$id[[ASCA_GOslims[i]]], TRUE)
  if(isTRUE(class(temp_beach)=="try-error")) {next} else {beach[[i]] = temp_beach}
}

test <- get_term_info_content(go, term_sets = beach)

sim_matrix <- get_sim_grid(
  ontology=go, 
  information_content=test,
  term_sets=beach)

rownames(sim_matrix) <- ASCA_GOslims
colnames(sim_matrix) <- ASCA_GOslims
#convert lower triangle of matrix to NA val including the diagonal
sim_matrix[lower.tri(sim_matrix, diag = TRUE)] <- NA

term_term <- melt(sim_matrix)
#remove same ID-same ID instances and remove term-term with 0 similarity
term_term <- term_term[which(term_term$Var1 != term_term$Var2 & term_term$value > 0.7),]

#convert terms to names

longterm_term <- data.frame(matrix(0,nrow(term_term),3))
colnames(longterm_term) <- c("Var1","Var2", "value")

for(i in 1:nrow(term_term)){
  longterm_term$Var1[i] <- slims[which(slims$GOid == term_term$Var1[i]),"Term"]
  longterm_term$Var2[i] <- slims[which(slims$GOid == term_term$Var2[i]),"Term"]
  longterm_term$value[i] <- term_term$value[i]
}

longterm_term$type <- "term-term"
colnames(longterm_term)[3] <- "semsimValue"

#protein term data
prot_term <- unique(goslim_protein[,2:3])
colnames(prot_term)[2] <- "GOid"

prot_term <- merge(prot_term,slims[,c("Term","GOid")], by = "GOid")
prot_term$type <- "protein-term"
prot_term <- prot_term[,-1]
colnames(prot_term) <- c("Var1", "Var2", "type")

edge_attb <- rbind.fill(longterm_term,prot_term)

write.csv(edge_attb, "~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/SameDayFCtoTemp/ASCAproteins_GOsemsim_edges/edge_attb_semsim0.7.csv", row.names = FALSE, quote = FALSE)


##############################
####ATTEMPT with non-slimmed terms to see if clustering in cytoscape will overcome the complexity
################################

#terms not in ontologyX DB
exclude_terms <- c("GO:0103046", "GO:0103025", "GO:0062023", "GO:0101031", "GO:0061844", "GO:0106037")
ex_ASCA_sig_GOids <- subset(ASCA_sig_GOids, !(ASCA_sig_GOids %in% exclude_terms))

beach <- list()
for(i in 1:length(ex_ASCA_sig_GOids)){
  temp_beach <- try(go$id[[ex_ASCA_sig_GOids[i]]], TRUE)
  if(isTRUE(class(temp_beach)=="try-error")) {next} else {beach[[i]] = temp_beach}
}

test <- get_term_info_content(go, term_sets = beach)

sim_matrix <- get_sim_grid(
  ontology=go, 
  information_content=test,
  term_sets=beach)

rownames(sim_matrix) <- ex_ASCA_sig_GOids
colnames(sim_matrix) <- ex_ASCA_sig_GOids
#convert lower triangle of matrix to NA val including the diagonal
sim_matrix[lower.tri(sim_matrix, diag = TRUE)] <- NA

term_term <- melt(sim_matrix)
#remove same ID-same ID instances and remove term-term with 0 similarity
term_term <- term_term[which(term_term$Var1 != term_term$Var2 & term_term$value > 0.7),]

term_term$type <- "term-term"
write.csv(term_term, "~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/SameDayFCtoTemp/ASCAproteins_GOsemsim_edges/term_term.csv", row.names = FALSE, quote = FALSE)
