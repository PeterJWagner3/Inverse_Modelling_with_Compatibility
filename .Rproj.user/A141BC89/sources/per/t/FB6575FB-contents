library(paleotree);		#install.packages("paleotree", dependencies=TRUE)
library(phangorn);		#install.packages("phangorn", dependencies=TRUE)
library(phytools);		#install.packages("ape", dependencies=TRUE)
library(ape);			#install.packages("ape", dependencies=TRUE)
library(gtools);		#install.packages("gtools", dependencies=TRUE)
library(BayesFactor);	#install.packages("BayesFactor", dependencies=TRUE)
library(phylobase);		#install.packages("phylobase", dependencies=TRUE)

# use the folder into which you put your general source R files
source("R/Compatibility.r");  #
source("R/Disparity.r");  #
source("R/General_Plot_Templates.r");  #
source("R/Nexus_File_Routines.r");  #
source("R/Paleo_Tree_extras.r");  #
source("R/Paleophylogeny_Routines.r");  #
source("R/Phylogeny_Simulations.r");  #
source("R/RevBayes_Setup.r");  #
source("R/Wagner_kluges.r");  #
source("R/Wagner_Stats_and_Probability_101.r");  #
UNKNOWN <- -11;
INAP <- -22;

desired_sample <- 25;
r <- 1.0	# expected samples per myr
p <- 1.0	# expected daughter taxa per myr
q <- 1.0	# expected extinctions per myr
pfr <- r/(r+q);	# expected proportion of taxa that we should sample
nTotalTaxa <- ceiling(desired_sample/pfr);
nchars <- 50;
states <- rep(2,nchars);	# let's just use binary characters for now
types <- rep(0,nchars);		# let's just unordered for now
max_changes <- 100;
# simulate a clade with sampling over time using parameters given above
fossil_record <- simFossilRecord(p,q,r,startTaxa=1,nTotalTaxa=nTotalTaxa);
sim_tree_info <- accersi_vector_tree_and_ancestors_from_paleotree_output(fossil_record);
chmatrix <- matrix(0,nrow=nTotalTaxa,ncol=nchars);
venn_tree_all <- transform_vector_tree_to_venn_tree(sim_tree_info$vector_tree);
branchings_all <- rep(1,length(sim_tree_info$vector_tree));
branchings_all[venn_tree_all[,1]] <- 0;
sim_compat <- evolve_compatibility_over_N_changes(N=max_changes,init_chmatrix=chmatrix,venn_tree=venn_tree_all,branchings=branchings_all,nchars,states,types,hidden_reversals=F,UNKNOWN=-11,INAP=-22);

par(pin=c(3,3));
mnx <- 5*floor(min(sim_compat[,1])/5);
mxx <- 5*ceiling(max(sim_compat[,1])/5);
mny <- 25*floor(min(sim_compat[,2])/25);
mxy <- 25*ceiling(max(sim_compat[,2])/25);
plot(NA,type='n',axes=FALSE,main="",xlab="Steps",ylab="Compatibility",xlim=c(mnx,mxx),ylim=c(mny,mxy));
specified_axis(axe=1,max_val=mxx,min_val=mnx,maj_break=10,med_break=5,min_break=1,orient=1);
specified_axis(axe=2,max_val=mxy,min_val=mny,maj_break=25,med_break=5,min_break=5,orient=2);
points(sim_compat[,1],sim_compat[,2],pch=22,bg="orange");

####
# That's nice, but it's just one run. Let's look at what 100 runs generates
ttl_reps <- 100;
sim_compat_inverse <- matrix(0,nrow=ttl_reps,ncol=1+(max_changes-nchars));
colnames(sim_compat_inverse) <- nchars:max_changes;
for (rr in 1:ttl_reps)	{
	fossil_record <- simFossilRecord(p,q,r,startTaxa=1,nTotalTaxa=nTotalTaxa);
	sim_tree_info <- accersi_vector_tree_and_ancestors_from_paleotree_output(fossil_record);
	chmatrix <- matrix(0,nrow=nTotalTaxa,ncol=nchars);
	venn_tree_all <- transform_vector_tree_to_venn_tree(sim_tree_info$vector_tree);
	branchings_all <- rep(1,length(sim_tree_info$vector_tree));
	branchings_all[venn_tree_all[,1]] <- 0;
	sim_compat <- evolve_compatibility_over_N_changes(N=max_changes,init_chmatrix=chmatrix,venn_tree=venn_tree_all,branchings=branchings_all,nchars,states,types,hidden_reversals=F,UNKNOWN=-11,INAP=-22,repl=rr);
	rsteps <- match(sim_compat[,1],as.numeric(colnames(sim_compat_inverse)));
	sim_compat_inverse[rr,rsteps] <- sim_compat[,2];
	}

cases <- colSums(sim_compat_inverse);
relv_steps <- as.numeric(names(cases[cases>0]));
cases <- cases[cases>0];
sim_compat_inverse <- sim_compat_inverse[,as.numeric(colnames(sim_compat_inverse)) %in% relv_steps];
#sim_compat_examples[,colSum(sim_compat_examples)>0]
mean_compatibility <- c();
for (cc in 1:ncol(sim_compat_inverse))
	mean_compatibility <- c(mean_compatibility,sum(sim_compat_inverse[sim_compat_inverse[,cc]>0,cc])/sum(sim_compat_inverse[,cc]>0));
changes <- min(as.numeric(colnames(sim_compat_inverse))):max_changes;
mnx <- 5*floor(min(changes)/5);
mxx <- 5*ceiling(max(changes)/5);
mny <- 25*floor(min(mean_compatibility)/25);
mxy <- 25*ceiling(max(mean_compatibility)/25);
plot(NA,type='n',axes=FALSE,main="",xlab="Steps",ylab="Compatibility",xlim=c(mnx,mxx),ylim=c(mny,mxy));
specified_axis(axe=1,max_val=mxx,min_val=mnx,maj_break=10,med_break=5,min_break=1,orient=1);
specified_axis(axe=2,max_val=mxy,min_val=mny,maj_break=25,med_break=5,min_break=5,orient=2);
points(changes,mean_compatibility,pch=22,bg="orange");

### OK, that's cool.  But what about comparing this to real data?
### First Step: get some real data!
real_data <- accersi_data_from_chosen_nexus_file();	# This lets you choose & read a nexus file.
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
compmat <- compatibility_matrix(chmatrix=char_matrix,states=nstates,types=chtypes,UNKNOWN,INAP);
char_comp <- rowSums(compmat)-1;
total_compat <- (sum(compmat)-nchars)/2;
taxa_scored <- count_scored_otu_per_character(chmatrix=char_matrix);
hist(char_comp,breaks=((min(char_comp)-1):max(char_comp)),main="Character Compatibilities",xlim=c(0,(nchars-1)),xlab="Compatibility",ylab="No. Characters");
#par(pin=c(3.5,3.5))
#plot((1:nchars),char_comp,xlab="Character No.",ylab="No. Compatibilities",ylim=c(0,nchars-1))
#points((1:nchars),nstates,pch=21,bg="red")
#points((1:nchars),taxa_scored,pch=22,bg="blue")
maxsteps_per_char <- accersi_maximum_parsimony_steps_per_character(chmatrix=char_matrix,states=nstates);
maxsteps <- round((sum(maxsteps_per_char)+sum(taxa_scored))/2,0);

# let's get a basic idea of how many steps we need first
fossil_record_inv <- simFossilRecord(p,q,r,startTaxa=1,nTotalTaxa=notu);	# fossil record for inverse modeling
sim_tree_info_inv <- accersi_vector_tree_and_ancestors_from_paleotree_output(fossil_record_inv);
chmatrix_inv <- char_matrix;
for (sp in 1:notu) chmatrix_inv[sp,chmatrix_inv[sp,]>0] <- 0;
venn_tree_inv <- transform_vector_tree_to_venn_tree(sim_tree_info_inv$vector_tree);
branchings_inv <- rep(1,length(sim_tree_info_inv$vector_tree));
branchings_inv[venn_tree_inv[,1]] <- 0;
sim_compat <- evolve_compatibility_over_N_changes(N=maxsteps,init_chmatrix=chmatrix_inv,venn_tree=venn_tree_inv,branchings=branchings_inv,nchars,states=nstates,types=chtypes,hidden_reversals=F,UNKNOWN=-11,INAP=-22);
needed_steps <- max(sim_compat[sim_compat[,2]>=total_compat,1]);
maxsteps <- ceiling(min(1.5*needed_steps,maxsteps));

# we probably greatlyy reduced maxsteps, which will speed up things.
# So, let's repeat the exercise, setting up matrices with the same dimensions as the real data.
ttl_reps <- 250;
sim_compat_inverse <- matrix(0,nrow=ttl_reps,ncol=1+(maxsteps-nchars));
colnames(sim_compat_inverse) <- nchars:maxsteps;
for (rr in 1:ttl_reps)	{
	fossil_record_inv <- simFossilRecord(p,q,r,startTaxa=1,nTotalTaxa=notu);
	sim_tree_info_inv <- accersi_vector_tree_and_ancestors_from_paleotree_output(fossil_record_inv);
	chmatrix_inv <- char_matrix;
	for (sp in 1:notu) chmatrix_inv[sp,chmatrix_inv[sp,]>0] <- 0;
	venn_tree_inv <- transform_vector_tree_to_venn_tree(sim_tree_info_inv$vector_tree);
	branchings_inv <- rep(1,length(sim_tree_info_inv$vector_tree));
	branchings_inv[venn_tree_inv[,1]] <- 0;
	sim_compat <- evolve_compatibility_over_N_changes(N=maxsteps,init_chmatrix=chmatrix_inv,venn_tree=venn_tree_inv,branchings=branchings_inv,nchars,states=nstates,types=chtypes,hidden_reversals=F,UNKNOWN=-11,INAP=-22,repl=rr);
	rsteps <- match(sim_compat[,1],as.numeric(colnames(sim_compat_inverse)));
	sim_compat_inverse[rr,rsteps] <- sim_compat[,2];
	}
# we have results!  So, what maximizes the probability of getting observed disparity?
sim_compat_inverse_init <- sim_compat_inverse;
sim_compat_inverse_relv <- sim_compat_inverse[,colSums(sim_compat_inverse)>0];
prob_observed <- mean_compat_inv <- var_compat_inv <-c();
for (nc in 1:ncol(sim_compat_inverse_relv))	{
	compats <- sim_compat_inverse_relv[sim_compat_inverse_relv[,nc]>0,nc];
	mean_compat_inv <- c(mean_compat_inv,mean(compats));
	var_compat_inv <- c(var_compat_inv,var(compats));
	dnorm(x=total_compat,mean=mean(compats),sd=sqrt(var(compats)))
	min_comp <- min(total_compat,compats);
	max_comp <- max(total_compat,compats);
	prob_observed <- c(prob_observed,dnorm(total_compat,mean=mean(compats),sd=sqrt(var(compats)))/sum(dnorm(min_comp:max_comp,mean=mean(compats),sd=sqrt(var(compats)))));
	}

mnx <- 25*floor(min(as.numeric(colnames(sim_compat_inverse_relv)))/25);
mxx <- 25*ceiling(max(as.numeric(colnames(sim_compat_inverse_relv)))/25);
mny <- 0.001*floor(min(prob_observed)/0.001);
mxy <- 0.001*ceiling(max(prob_observed)/0.001);
ml_steps <- as.numeric(colnames(sim_compat_inverse_relv))[prob_observed==max(prob_observed)];
plot_title <- paste("ML Steps = ",ml_steps,sep="");
plot(NA,type='n',axes=FALSE,main=plot_title,xlab="Steps",ylab="P[Obs. Compatibility | Steps]",xlim=c(mnx,mxx),ylim=c(mny,mxy));
specified_axis(axe=1,max_val=mxx,min_val=mnx,maj_break=50,med_break=25,min_break=51,orient=1);
specified_axis(axe=2,max_val=mxy,min_val=mny,maj_break=0.001,med_break=0.0005,min_break=0.0001,orient=2);
points(as.numeric(colnames(sim_compat_inverse_relv)),prob_observed,pch=22,bg="orange");
ml_alpha <- (ml_steps/sum(branchings_inv))/sum(nstates-1);

ml_steps/(nchars*notu)
ttt <- (-1:10)/20
alphas <- c(0.167,0.121,0.44,0.21,0.21,0.226,0.107,0.231,0.13,0.216,0.255,0.28,0.32)
hist(alphas,xlab="per-taxon change per-character",breaks=ttt);
c("bivalves","mammals","mammals","trilobites","gastropods","echinoderms","gastropods","echinoderms","dinosaurs","mammals","brachiopods","fish")

