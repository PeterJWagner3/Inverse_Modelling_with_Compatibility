library(paleotree);		#install.packages("paleotree", dependencies=TRUE)
library(phangorn);		#install.packages("phangorn", dependencies=TRUE)
library(phytools);		#install.packages("ape", dependencies=TRUE)
library(ape);			#install.packages("ape", dependencies=TRUE)
library(gtools);		#install.packages("gtools", dependencies=TRUE)
library(BayesFactor);	#install.packages("BayesFactor", dependencies=TRUE)
library(phylobase);		#install.packages("phylobase", dependencies=TRUE)

# use the folder into which you put your general source R files
common_source_folder <- "~/Documents/R_Projects/Common_R_Source_Files/";	
source(paste(common_source_folder,"Wagner_R_Functions_for_Analytical_Class.r",sep=""));  #
UNKNOWN <- -11;
INAP <- -22;

desired_sample <- 25;
r <- 1.0	# expected samples per myr
p <- 1.0	# expected daughter taxa per myr
q <- 1.0	# expected extinctions per myr
pfr <- r/(r+q);	# expected proportion of taxa that we should sample
notu <- nTotalTaxa <- ceiling(desired_sample/pfr);
nTotalTaxa <- notu <- 25;
nchars <- 25;
states <- rep(2,nchars);	# let's just use binary characters for now
types <- rep(0,nchars);		# let's just unordered for now
ttl_steps <- 33;
# simulate a clade with sampling over time using parameters given above
fossil_record <- simFossilRecord(p,q,r,startTaxa=1,nTotalTaxa=nTotalTaxa);
simulated_history <- data.frame(fossilRecord2fossilTaxa(fossil_record),stringsAsFactors = F);	
sim_durations <- cbind(simulated_history$orig.time,simulated_history$ext.time);
sim_tree_info <- accio_vector_tree_and_ancestors_from_paleotree_output(fossil_record);
ctree_all <- ctree <- sim_tree_info$vector_tree;
simulated_ancestors <- sim_tree_info$node_ancestors;
durations_all <- durations <- rbind(sim_durations,sim_durations[simulated_ancestors,]);
chmatrix <- matrix(0,nrow=nTotalTaxa,ncol=nchars);
venn_tree_all <- accio_venn_tree_from_vector_tree(sim_tree_info$vector_tree);
branchings_all <- rep(1,length(sim_tree_info$vector_tree));
branchings_all[venn_tree_all[,1]] <- 0;
char_matrix_sim <- evolve_N_steps_on_tree(N=ttl_steps,init_chmatrix=chmatrix,venn_tree=venn_tree_all,branchings=branchings_all,nchars,states,types,hidden_reversals=F,UNKNOWN=-11,INAP=-22);
#compendium <- accio_fossil_record_summary_from_paleotree_output(fossil_record);
#strat_ranges <- compendium$strat_ranges
#write(sim_tree_info$vector_tree,file="Simulated_Tree.txt");
#strat_ranges$LAD[strat_ranges$FAD==0] <- (durations_all[strat_ranges$FAD==0,1]+durations_all[strat_ranges$FAD==0,2])/2;
#strat_ranges$FAD[strat_ranges$FAD==0] <- (durations_all[strat_ranges$FAD==0,1]+durations_all[strat_ranges$FAD==0,2])/2;
#rownames(char_matrix_sim) <- paste("taxon_",1:notu,sep="");
#write.table(char_matrix_sim,file="thirty_three.txt");
#latest <- round(min(strat_ranges$FAD),3);
#fossil_intervals <- data.frame(taxon=as.character(paste("taxon_",1:notu,sep="")),min=as.numeric(round(strat_ranges$FAD,3)-latest),max=as.numeric(round(strat_ranges$FAD,3)-latest),stringsAsFactors = F);
#write.table(fossil_intervals,file="simulated_taxon_intervals.tsv",row.names = F);
dissimilarities_sim <- pairwise_dissimilarity(chmatrix=char_matrix_sim,states=states,types=types,weight_ordered=F,polymorphs=F,UNKNOWN=-11,INAP=-22)
# first, let's get cumulative disparity
cumulative_disparity_sim <- c();
for (nn in 5:notu)
	cumulative_disparity_sim <- c(cumulative_disparity_sim,mean(extract_off_diagonal(orig_matrix=dissimilarities_sim[1:nn,1:nn])));
mxx <- notu;
mnx <- 5;
mxy <- 0.05*ceiling(max(cumulative_disparity_sim)/0.05);
mny <- 0.0;
#plot(5:notu,cumulative_disparity_sim,xlab="",ylab="",pch=21,bg="red")
plot_title <- "Simulated Cumulative Disparity at N Species";
plot(NA,type='n',axes=FALSE,main=plot_title,xlab="Total Species Evolved",ylab="Cumulative Disparity",xlim=c(mnx,mxx),ylim=c(mny,mxy));
specified_axis(axe=1,max_val=mxx,min_val=mnx,maj_break=10,med_break=5,min_break=1,orient=1);
specified_axis(axe=2,max_val=mxy,min_val=mny,maj_break=0.05,med_break=0.01,min_break=0.01,orient=2);
points(5:notu,cumulative_disparity_sim,pch=21,bg="red");
# again, that is just one: let's do it with a few.
#YOU DO THIS PART
ttl_reps <- 100;
sim_disparity <- matrix(0,nrow=ttl_reps,ncol=1+(notu-5));
for (rr in 1:ttl_reps)	{
	fossil_record <- simFossilRecord(p,q,r,startTaxa=1,nTotalTaxa=nTotalTaxa);
	sim_tree_info <- accio_vector_tree_and_ancestors_from_paleotree_output(fossil_record);
	chmatrix <- matrix(0,nrow=nTotalTaxa,ncol=nchars);
	venn_tree_all <- accio_venn_tree_from_vector_tree(sim_tree_info$vector_tree);
	branchings_all <- rep(1,length(sim_tree_info$vector_tree));
	branchings_all[venn_tree_all[,1]] <- 0;
	char_matrix_sim <- evolve_N_steps_on_tree(N=ttl_steps,init_chmatrix=chmatrix,venn_tree=venn_tree_all,branchings=branchings_all,nchars,states,types,hidden_reversals=F,UNKNOWN=-11,INAP=-22,repl=rr);
	dissimilarities_sim <- pairwise_dissimilarity(chmatrix=char_matrix_sim,states=states,types=types,weight_ordered=F,polymorphs=F,UNKNOWN=-11,INAP=-22)
	# first, let's get cumulative disparity
	cumulative_disparity_sim <- c();
	for (nn in 5:notu)
		cumulative_disparity_sim <- c(cumulative_disparity_sim,mean(extract_off_diagonal(orig_matrix=dissimilarities_sim[1:nn,1:nn])));
	sim_disparity[rr,] <- cumulative_disparity_sim;
	}
typical_disparity_sim <- colSums(sim_disparity)/ttl_reps;
plot_title <- "Simulated Cumulative Disparity at N Species (100 Runs)";
plot(NA,type='n',axes=FALSE,main=plot_title,xlab="Total Species Evolved",ylab="Cumulative Disparity",xlim=c(mnx,mxx),ylim=c(mny,mxy));
specified_axis(axe=1,max_val=mxx,min_val=mnx,maj_break=10,med_break=5,min_break=1,orient=1);
specified_axis(axe=2,max_val=mxy,min_val=mny,maj_break=0.05,med_break=0.01,min_break=0.01,orient=2);
points(5:notu,typical_disparity_sim,pch=22,bg="red");

# Now, let's do this with real data.
real_data <- accio_data_from_chosen_nexus_file();	# This lets you choose & read a nexus file.
#strat_ranges_file <- file.choose();
#strat_ranges <- read.csv(strat_ranges_file,header = T,stringsAsFactors = F);
strat_ranges <- read.csv(file.choose(),header = T,stringsAsFactors = F);
char_matrix <- real_data$Matrix;					# character matrix
nchars <- ncol(char_matrix);						# number of characters
nstates <- real_data$States;						# states per character
taxon_names <- real_data$OTUs;						# taxon names
init_types <- real_data$State_Order;
chtypes <- 1*vector(length=nchars);
chtypes[tolower(init_types)=="ordered"] <- 1;
#chtypes <- rep(0,nchars)
## remove any invariant characters
char_matrix <- char_matrix[,nstates>1];
chtypes <- chtypes[nstates>1];
nstates <- nstates[nstates>1];
nchars <- ncol(char_matrix);						# number of characters
notu <- nrow(char_matrix);
# to make our lives a little easier, let's turn all polymorphics into UNKNOWN
for (sp in 1:notu)	
	for (ch in 1:nchars)	
		if (char_matrix[sp,ch]<0 && char_matrix[sp,ch]!=INAP)	
			char_matrix[sp,ch] <- UNKNOWN;
dissimilarities_obs <- pairwise_dissimilarity(chmatrix=char_matrix,states=nstates,types=chtypes,weight_ordered=F,polymorphs=F,UNKNOWN=-11,INAP=-22)
oldest <- match(min(strat_ranges$FA),strat_ranges$FA)
dist_from_eldest <- dissimilarities_obs[1,];
strat_bins <- sort(unique(strat_ranges$FA));
appearance_order <- match(strat_ranges$FA,strat_bins);
sequence_order <- (1:notu)[order(appearance_order,dist_from_eldest)];
cumulative_disparity_obs <- c();
for (nn in 5:notu)	{
	cumul_richness <- sequence_order[1:nn];
	off_diag <- extract_off_diagonal(orig_matrix=dissimilarities_obs[cumul_richness,cumul_richness]);
	off_diag <- off_diag[!is.nan(off_diag)];
	cumulative_disparity_obs <- c(cumulative_disparity_obs,mean(off_diag));
	}
plot_title <- "Observed Cumulative Disparity at N Species";
mxy <- 0.05*ceiling(max(cumulative_disparity_obs)/0.05);
mxx <- notu;
plot(NA,type='n',axes=FALSE,main=plot_title,xlab="Total Species Evolved",ylab="Cumulative Disparity",xlim=c(mnx,mxx),ylim=c(mny,mxy));
specified_axis(axe=1,max_val=mxx,min_val=mnx,maj_break=5,med_break=5,min_break=1,orient=1);
specified_axis(axe=2,max_val=mxy,min_val=mny,maj_break=0.05,med_break=0.01,min_break=0.01,orient=2);
points((5:notu),cumulative_disparity_obs,pch=21,bg="orange");

### ttl_steps here shoud be you ml_steps from this mornings exercise!!!
# ml_steps=210
#char_matrix_sim <- evolve_N_steps_on_tree(N=ttl_steps,init_chmatrix=chmatrix,venn_tree=venn_tree_all,branchings=branchings_all,nchars,states,types,hidden_reversals=F,UNKNOWN=-11,INAP=-22,repl=rr);
ttl_reps <- 250;
ml_steps <- 210;	# Use the number from the prior exercise
inv_disparity <- matrix(0,nrow=ttl_reps,ncol=1+(notu-5));
for (rr in 1:ttl_reps)	{
	# notu should be the number of taxa in your matrix!
	fossil_record <- simFossilRecord(p,q,r,startTaxa=1,nTotalTaxa=notu);
	sim_tree_info_inv <- accio_vector_tree_and_ancestors_from_paleotree_output(fossil_record);
	chmatrix_inv <- char_matrix[sequence_order,];	# make sure that unknowns and inaps are in the right order!
	for (sp in 1:notu) chmatrix_inv[sp,chmatrix_inv[sp,]>0] <- 0;
	venn_tree_inv <- accio_venn_tree_from_vector_tree(sim_tree_info_inv$vector_tree);
	branchings_inv <- rep(1,length(sim_tree_info_inv$vector_tree));
	branchings_inv[venn_tree_inv[,1]] <- 0;
	char_matrix_inv <- evolve_N_steps_on_tree(N=ml_steps,init_chmatrix=chmatrix_inv,venn_tree=venn_tree_inv,branchings=branchings_inv,nchars,states=nstates,types=chtypes,hidden_reversals=F,UNKNOWN=-11,INAP=-22,repl=rr);
	dissimilarities_inv <- pairwise_dissimilarity(chmatrix=char_matrix_inv,states=nstates,types=chtypes,weight_ordered=F,polymorphs=F,UNKNOWN=-11,INAP=-22);
	# first, let's get cumulative disparity
	cumulative_disparity_inv <- c();
	for (nn in 5:notu)	{
		off_diag <- extract_off_diagonal(orig_matrix=dissimilarities_inv[1:nn,1:nn]);
		off_diag <- off_diag[!is.nan(off_diag)];
		cumulative_disparity_inv <- c(cumulative_disparity_inv,mean(off_diag));
		}
		
	inv_disparity[rr,] <- cumulative_disparity_inv;
	}
inv_cum_disparity_sorted <- inv_disparity;
inv_cum_disparity <- colSums(inv_disparity)/ttl_reps;
for (ii in 1:ncol(inv_disparity))	{
	inv_cum_disparity_sorted[,ii] <- sort(inv_disparity[,ii]);
	}
par(pin=c(4,3));
mxy <- 0.05*ceiling(max(inv_cum_disparity_sorted)/0.05);
mny <- 0;
mxx <- notu;
mnx <- 5;
plot(NA,type='n',axes=FALSE,main=plot_title,xlab="Total Species Evolved",ylab="Cumulative Disparity",xlim=c(mnx,mxx),ylim=c(mny,mxy));
specified_axis(axe=1,max_val=mxx,min_val=mnx,maj_break=5,med_break=5,min_break=1,orient=1);
specified_axis(axe=2,max_val=mxy,min_val=mny,maj_break=0.05,med_break=0.01,min_break=0.01,orient=2);
points((5:notu),inv_cum_disparity_sorted[round(0.5*nrow(inv_cum_disparity_sorted),0),],pch=22,bg="sienna");
lines((5:notu),inv_cum_disparity_sorted[round(0.975*nrow(inv_cum_disparity_sorted),0),]);
lines((5:notu),inv_cum_disparity_sorted[round(0.025*nrow(inv_cum_disparity_sorted),0),]);
points((5:notu),cumulative_disparity_obs,pch=21,bg="orange");
