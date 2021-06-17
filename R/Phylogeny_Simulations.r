
# ROUTINES RELATED INVERSE MODELLING WITH COMPATIBILITY
simulate_cladogram <- function(otu,tree_type=c("vector","matrix"))	{
htu <- otu+1
branch_durations <- vector(length=(2*otu)-1)
branch_durations[1:2] <- 1
if (tree_type=="matrix")	{
	cladogram <- c(1,2)
	for (s in 3:otu)	{
		htu <- htu+1
		split <- ceiling(runif(1)/(1/(s-1)))
		new_cl <- c(split,s)
		cladogram <- rbind(cladogram,new_cl)
		# find where the tree needs to be altered
		alter <- which(cladogram==split,arr.ind=TRUE)
		alter <- alter[order(alter[,1],decreasing=FALSE),]
		# replace split species with htu number
		cladogram[alter[1,1],alter[1,2]] <- htu
		cladogram[alter[1,1],] <- sort(cladogram[alter[1,1],])
		branch_durations[htu] <- branch_durations[split]
		branch_durations[split] <- 0
		branch_durations[array(cladogram)[array(cladogram)<=otu]] <- branch_durations[array(cladogram)[array(cladogram)<=otu]]+1
		}
	} else	{
	cladogram <- rep(-1,(2*otu)-1)
	cladogram[1:2] <- htu
	branch_durations <- vector(length=(2*otu)-1)
	branch_durations[1:2] <- 1
	for (s in 3:otu)	{
		htu <- htu+1
		split <- ceiling(runif(1)/(1/(s-1)))
		cladogram[htu] <- cladogram[split]
		cladogram[split] <- cladogram[s] <- htu
		branch_durations[htu] <- branch_durations[split]
		branch_durations[split] <- 0
		branch_durations[1:s] <- branch_durations[1:s]+1
#		s <- s+1
		}
	}
output <- list(cladogram,branch_durations)
names(output) <- c("Cladogram","Rel_Branch_Lengths")
return (output)
}

# routine to evolve S contemporaneous (e.g., extant) taxa
evolve_to_standing_richness_S <- function(S,lambda,mu,bifurcation=1,temp_prec=0.1)	{
# cumulative and standing diversity
spc <- cumulative <- standing <- 1
birth <- c(temp_prec)
life <- temp_prec*ceiling((-log(runif(1))/mu)/temp_prec)
death <- c(birth[spc]+life)
ancestral <- c(0)
remaining <- c(1)
dbr <- 1
while (max(standing) < S)	{
	spc <- remaining[1]
	life <- abs((death[spc]-birth[spc]))
	daughters <- rpois(1,(lambda*life))
	if (daughters > 0)	{
		if (bifurcation==1)	{
			splits <- birth[spc]+sort(life*runif(daughters))
			splits <- temp_prec*ceiling(splits/temp_prec)
			new_lines <- resid_lines <- c()
			resid_lines <- c(resid_lines,cumulative + (2*(1:daughters)-1))
			new_lines <- c(new_lines,cumulative + (2*(1:daughters)))
			prior_anc <- spc
#			ancestral <- c(ancestral,rep(spc,2*daughters))
#			birth <- c(birth,rep(0,2*daughters))
#			death <- c(death,rep(0,2*daughters))
			for (i in 1:daughters)	{
				rb <- resid_lines[i]
				lb <- new_lines[i]
				ancestral <- c(ancestral,prior_anc,prior_anc)
				birth <- c(birth,splits[i],splits[i])
				if (i < daughters)  {
					death <- c(death,splits[i+1]-temp_prec)
					}	else	{
					death <- c(death,death[spc])
					death[spc] <- splits[1]-temp_prec
					}
				life <- temp_prec*ceiling((-log(runif(1))/mu)/temp_prec)
				death <- c(death,(splits[i] + life))
				prior_anc <- rb	# the right branch leads to the next split
				}
			remaining <- c(remaining,new_lines)
			}	else	{
			# for budding model, just get appearance and death of each daughter
			birth <- c(birth,birth[spc]+sort(life*runif(daughters)))
			ancestral <- c(ancestral,rep(spc,daughters))
			for (d in 1:daughters)	{
				life <- temp_prec*ceiling((-log(runif(1))/mu)/temp_prec)
				death[d+cumulative] <- birth[d+cumulative]+life
				}
			remaining <- c(remaining,cumulative+(1:daughters))
			}
		cumulative <- max(remaining)
		ranges <- cbind(birth,death)
		standing <- tally_richness_from_continuous_ranges(ranges,temp_prec)
		}
	if (length(remaining) > 1)	{
		remaining <- remaining[2:length(remaining)]
		}	else if (max(standing) < S)	{
		dbr <- dbr+1
		# reboot..... :-<     
		spc <- cumulative <- standing <- 1
		birth <- c(0.1)
		life <- temp_prec*ceiling((-log(runif(1))/mu)/temp_prec)
		death <- c(birth[spc]+life)
		ancestral <- c(0)
		remaining <- c(1)
		}
	}

sim_species_present <- accersi_species_present_in_bins(ranges,temp_prec=0.1)
if (!is.na(match(S,standing)))	{
	present <- match(S,standing)
	}	else	{
	present <- match(max(standing),standing)
	}
sampled <- sim_species_present[present,]

if (length(sampled)>S)	{
	xx <- sampled[order(birth[sampled],decreasing=FALSE)]
	sampled <- sort(xx[1:S],decreasing=FALSE)
	}

venn_tree_info <- accersi_simulated_venn_tree_info(sampled,ancestral)
venn_tree <- venn_tree_info$Venn_Tree
branchings <- venn_tree_info$Branchings
branchings[S+1] <- 0
vtree <- transform_venn_tree_to_vector_tree(venn_tree)
divergence_dates <- birth[venn_tree_info$First_Taxon]
divergence_dates[S+1] <- max(divergence_dates)
divergence_dates <- divergence_dates-(min(divergence_dates))
divergence_dates[S+1] <- 0
#print(branchings[1:S])
#print(divergence_dates)
output <- list(vtree,venn_tree,branchings,divergence_dates)
names(output) <- c("Vector_Tree","Venn_Tree","Branchings","Divergence_Dates")
return(output)
}

# routine to evolve S contemporaneous (e.g., extant) taxa
MBL_style_to_standing_richness_S <- function(S,lambda,mu,bifurcation=1,temp_prec=0.1)	{
# cumulative and standing diversity
spc <- cumulative <- standing <- 1
birth <- c(temp_prec)
life <- temp_prec*ceiling((-log(runif(1))/mu)/temp_prec)
death <- c(birth[spc]+life)
ancestral <- c(0)
remaining <- c(1)
dbr <- 1
while (max(standing) < S)	{
	spc <- remaining[1]
	life <- abs((death[spc]-birth[spc]))
	daughters <- rpois(1,(lambda*life))
	if (daughters > 0)	{
		if (bifurcation==1)	{
			splits <- birth[spc]+sort(life*runif(daughters))
			splits <- temp_prec*ceiling(splits/temp_prec)
			new_lines <- resid_lines <- c()
			resid_lines <- c(resid_lines,cumulative + (2*(1:daughters)-1))
			new_lines <- c(new_lines,cumulative + (2*(1:daughters)))
			prior_anc <- spc
#			ancestral <- c(ancestral,rep(spc,2*daughters))
#			birth <- c(birth,rep(0,2*daughters))
#			death <- c(death,rep(0,2*daughters))
			for (i in 1:daughters)	{
				rb <- resid_lines[i]
				lb <- new_lines[i]
				ancestral <- c(ancestral,prior_anc,prior_anc)
				birth <- c(birth,splits[i],splits[i])
				if (i < daughters)  {
					death <- c(death,splits[i+1]-temp_prec)
					}	else	{
					death <- c(death,death[spc])
					death[spc] <- splits[1]-temp_prec
					}
				life <- temp_prec*ceiling((-log(runif(1))/mu)/temp_prec)
				death <- c(death,(splits[i] + life))
				prior_anc <- rb	# the right branch leads to the next split
				}
			remaining <- c(remaining,new_lines)
			}	else	{
			# for budding model, just get appearance and death of each daughter
			birth <- c(birth,birth[spc]+sort(life*runif(daughters)))
			ancestral <- c(ancestral,rep(spc,daughters))
			for (d in 1:daughters)	{
				life <- temp_prec*ceiling((-log(runif(1))/mu)/temp_prec)
				death[d+cumulative] <- birth[d+cumulative]+life
				}
			remaining <- c(remaining,cumulative+(1:daughters))
			}
		cumulative <- max(remaining)
		ranges <- cbind(birth,death)
		standing <- tally_richness_from_continuous_ranges(ranges,temp_prec)
		}
	if (length(remaining) > 1)	{
		remaining <- remaining[2:length(remaining)]
		}	else if (max(standing) < S)	{
		dbr <- dbr+1
		# reboot..... :-<     
		spc <- cumulative <- standing <- 1
		birth <- c(0.1)
		life <- temp_prec*ceiling((-log(runif(1))/mu)/temp_prec)
		death <- c(birth[spc]+life)
		ancestral <- c(0)
		remaining <- c(1)
		}
	}
output <- cbind(ancestral,birth,death)
#output <- list(vtree,venn_tree,branchings,divergence_dates)
#names(output) <- c("Vector_Tree","Venn_Tree","Branchings","Divergence_Dates")
return(output)
}

logistic_to_standing_richness_S <- function(R=0.75,mu=0.50,K=15,bifurcation=0,temp_prec = 0.1)	{
spc <- cumulative <- standing <- 1
birth <- c(temp_prec)
life <- temp_prec*ceiling((-log(runif(1))/mu)/temp_prec)
death <- c(birth[spc]+life)
ancestral <- c(0)
remaining <- c(1)
dbr <- 1
bin <- 1	# general stage number
started <- temp_prec
elapsed <- 2*temp_prec
while (max(standing) < K)	{
	lambda <- expected_origination_given_logistic_constant_extinction(mu,R,K,S=standing)[1]
	for (s in 1:standing)	{
		spc <- remaining[s]
		if (birth[spc]<=started)	{
			daughters <- rpois(1,(lambda*temp_prec))
			} else	{
			daughters <- rpois(1,(lambda*(elapsed-birth[spc])))
			}	# if taxon appears partway through time slice, then reduce chance of daughters
		if (daughters > 0)	{
			if (bifurcation==1)	{
				splits <- started+runif(1)*temp_prec
				new_lines <- resid_lines <- c()
				resid_lines <- c(resid_lines,cumulative + (2*(1:daughters)-1))
				new_lines <- c(new_lines,cumulative + (2*(1:daughters)))
				prior_anc <- spc
				for (i in 1:daughters)	{
					rb <- resid_lines[i]
					lb <- new_lines[i]
					ancestral <- c(ancestral,prior_anc,prior_anc)
					birth <- c(birth,splits[i],splits[i])
					if (i < daughters)  {
						death <- c(death,splits[i+1]-temp_prec)
						}	else	{
						death <- c(death,death[spc])
						death[spc] <- splits[1]-temp_prec
						}
					life <- temp_prec*ceiling((-log(runif(1))/mu)/temp_prec)
					death <- c(death,(splits[i] + life))
					prior_anc <- rb	# the right branch leads to the next split
					}
				remaining <- c(remaining,new_lines)
			}	else	{
			# for budding model, just get appearance and death of each daughter
				xx <- max(started,birth[spc])
				birth <- c(birth,xx+sort(temp_prec*runif(daughters)))
				ancestral <- c(ancestral,rep(spc,daughters))
				for (d in 1:daughters)	{
					life <- temp_prec*ceiling((-log(runif(1))/mu)/temp_prec)
					death[d+cumulative] <- birth[d+cumulative]+life
					}
				remaining <- c(remaining,cumulative+(1:daughters))
				}
			cumulative <- max(remaining)
			ranges <- cbind(birth,death)
			standing <- tally_richness_from_continuous_ranges(ranges,temp_prec)
			}
		}
	remaining <- remaining[death[remaining]>elapsed]
	if (length(remaining) == 0)	{
#		remaining <- remaining[2:length(remaining)]
#		}	else if (max(standing) < S)	{
		dbr <- dbr+1	# trials needed.
		# reboot..... :-<     
		spc <- cumulative <- standing <- 1
		birth <- c(0.1)
		life <- temp_prec*ceiling((-log(runif(1))/mu)/temp_prec)
		death <- c(birth[spc]+life)
		ancestral <- 0
		remaining <- 1
		started <- 0
		elapsed <- temp_prec
		} # end reboot
	standing <- length(remaining)
	started <- started+temp_prec
	elapsed <- elapsed+temp_prec
	}
output <- cbind(ancestral,birth,death)
return(output)
}

# routine to evolve S taxa sampled from throughout the clade's history
evolve_to_sampled_richness_S <- function(S,lambda,mu,freqrat,bifurcation=1,temp_prec=0.1)	{
# cumulative and standing diversity
maxfinds <- dbr <- spc <- cumulative <- standing <- 1
birth <- c(temp_prec)
life <- temp_prec*ceiling((-log(runif(1))/mu)/temp_prec)
death <- c(birth[spc]+life)
ancestral <- c(0)
remaining <- c(1)
found <- 0	# this will keep track of how many species we have sampled
sampled <- c()	# list of sampled species
fa <- la <- occurrences <- vector(length=1)	# occurrences per species
all_finds <- matrix(0,1,maxfinds)
keep_going <- TRUE
while (keep_going)	{
	spc <- remaining[1]
	life <- abs((death[spc]-birth[spc]))
	daughters <- rpois(1,(lambda*life))
	fossils <- rpois(1,freqrat*life)
	if (fossils > 0)	finds <-  sort(birth[spc]+runif(fossils)*life)
	if (daughters > 0)	{
		if (bifurcation==1)	{
			splits <- birth[spc]+sort(life*runif(daughters))
			splits <- temp_prec*ceiling(splits/temp_prec)
			new_lines <- resid_lines <- c()
			resid_lines <- c(resid_lines,cumulative + (2*(1:daughters)-1))
			new_lines <- c(new_lines,cumulative + (2*(1:daughters)))
			prior_anc <- spc
#			ancestral <- c(ancestral,rep(spc,2*daughters))
#			birth <- c(birth,rep(0,2*daughters))
#			death <- c(death,rep(0,2*daughters))
			for (i in 1:daughters)	{
				rb <- resid_lines[i]
				lb <- new_lines[i]
				ancestral <- c(ancestral,prior_anc,prior_anc)
				birth <- c(birth,splits[i],splits[i])
				if (i < daughters)  {
					death <- c(death,splits[i+1]-temp_prec)
					}	else	{
					death <- c(death,death[spc])
					death[spc] <- splits[1]-temp_prec
					}
				life <- temp_prec*ceiling((-log(runif(1))/mu)/temp_prec)
				death <- c(death,(splits[i] + life))
				prior_anc <- rb	# the right branch leads to the next split
				}
			# distribute any finds among ancestors
			if (fossils > 0)	{
				newsplits <- c(birth[spc],splits)
				anagenetics <- c(spc,resid_lines)
				occurrences <- c(occurrences,rep(0,daughters*2))
				fa <- c(fa,rep(0,daughters*2))
				la <- c(la,rep(0,daughters*2))
				all_finds <- rbind(all_finds,matrix(0,daughters*2,maxfinds))
				for (f in 1:fossils)	{
					ap <- anagenetics[sum(finds[f]>=newsplits)]
					occurrences[ap] <- occurrences[ap]+1	# species number
					if (maxfinds < occurrences[ap])	{
						maxfinds <- occurrences[ap]
						all_finds <- cbind(all_finds,rep(0,length(occurrences)))
						}
					if (occurrences[ap]==1)	{
						sampled <- c(sampled,ap)
						all_finds[ap,1] <- fa[ap] <- la[ap] <- finds[f]
						}	else	{
						all_finds[ap,occurrences[ap]] <- la[ap] <- finds[f]
						}
					}
				}
			### add something to divide samples here.
			remaining <- c(remaining,new_lines)
			}	else	{
			# for budding model, just get appearance and death of each daughter
#			birth <- c(birth,birth[spc]+sort(life*runif(daughters)))
			birth <- c(birth,birth[spc]+sort(temp_prec*ceiling((life*runif(daughters))/temp_prec)))
			ancestral <- c(ancestral,rep(spc,daughters))
			for (d in 1:daughters)	{
				life <- temp_prec*ceiling((-log(runif(1))/mu)/temp_prec)
				death[d+cumulative] <- birth[d+cumulative]+life
				}
			# update dimensions of information arrays
			remaining <- c(remaining,cumulative+(1:daughters))
			added <- rep(0,daughters)
			occurrences <- c(occurrences,added)
			fa <- c(fa,added)
			la <- c(la,added)
			all_finds <- rbind(all_finds,matrix(0,daughters,maxfinds))
			}
		cumulative <- max(remaining)
		ranges <- cbind(birth,death)
		standing <- tally_richness_from_continuous_ranges(ranges,temp_prec)
		}
	if ((daughters == 0 || bifurcation==0) && fossils > 0)	{
		occurrences[spc] <- fossils
		fa[spc] <- min(finds)
		la[spc] <- max(finds)
		sampled <- c(sampled,spc)
		if (maxfinds < occurrences[spc])	{
			added <- occurrences[spc] - maxfinds
			maxfinds <- occurrences[spc]
			all_finds <- cbind(all_finds,matrix(0,length(occurrences),added))
			all_finds[spc,] <- finds
			}	else	{
			all_finds[spc,1:fossils] <- finds	
			}
		}
	if (length(remaining) > 1)	{
		remaining <- remaining[2:length(remaining)]
		}	else if (max(standing) < S)	{
		dbr <- dbr+1
		# reboot..... :-<     
		spc <- cumulative <- standing <- 1
		birth <- c(0.1)
		life <- temp_prec*ceiling((-log(runif(1))/mu)/temp_prec)
		death <- c(birth[spc]+life)
		ancestral <- c(0)
		remaining <- c(1)
		# clear sampling vectors
		sampled <- c()	# list of sampled species
		fa <- la <- occurrences <- vector(length=1)	# occurrences per species
		maxfinds <- 1
		all_finds <- matrix(0,1,maxfinds)
		}
	
	found <- length(sampled)
	relv_fas <- sort(fa[sampled])
	# stop when it is not possible to find any species with older fossils after we have S
	if (found>=S)
		if (relv_fas[S]<min(birth[remaining]))
			keep_going <- FALSE
	}

# the simulations usually over-sample to avoid arbitrary cutoff wonkiness
#	So, reduce to just S species
if (length(sampled)>S)	{
	xx <- sampled[order(fa[sampled],decreasing=FALSE)]
	sampled <- sort(xx[1:S],decreasing=FALSE)
	}

venn_tree_info <- accersi_simulated_venn_tree_info(sampled,ancestral)
venn_tree <- venn_tree_info$Venn_Tree
branchings <- venn_tree_info$Branchings
branchings[S+1] <- 0
vtree <- transform_venn_tree_to_vector_tree(venn_tree)
divergence_dates <- birth[venn_tree_info$First_Taxon]
if (sampled[1]==1)	{
	divergence_dates <- divergence_dates-divergence_dates[1]
	}	else	{
	divergence_dates <- divergence_dates-(min(divergence_dates[divergence_dates>divergence_dates[S+1]])-1)
	}
divergence_dates[S+1] <- 0
durations <- cbind(birth[sampled],death[sampled])
strat_ranges <- cbind(temp_prec*floor(fa[sampled]/temp_prec),temp_prec*ceiling(la[sampled]/temp_prec))
#print(branchings[1:S])
#print(divergence_dates)
output <- list(vtree,venn_tree,branchings,divergence_dates,durations,strat_ranges,occurrences[sampled])
names(output) <- c("Vector_Tree","Venn_Tree","Branchings","Divergence_Dates","Durations","Stratigraphic_Ranges","Finds")
return(output)
}

# routine to convert simulation into phylogeny
accersi_simulated_venn_tree_info <- function(sampled,ancestral)	{
# first, get the venn tree for all simulated taxa leading to sampled
raw_venn_tree <- accersi_raw_venn_tree_from_ancestor_list(sampled,ancestral=ancestral)
#write.table(raw_venn_tree,"raw_tree.xls",sep="\t")
# now, reduce it to unique nodes, adding redundant nodes to branch lengths
#	for species, this will be venn_tree nodes in which just that specees appears
#	for nodes, this will duplicate vectors
otus <- length(raw_venn_tree[1,])
reduction <- remove_singletons_from_venn_tree(venn_tree=raw_venn_tree)
#write.table(reduction$Condensed_Venn_Tree,"raw_clades_only.xls",sep="\t")
condensed_raw_venn_tree <- unique(reduction$Condensed_Venn_Tree)
#write.table(condensed_raw_venn_tree,"reduced_raw_tree.xls",sep="\t")
#condensed_raw_venn_tree <- reduction$Condensed_Venn_Tree
unique_nodes <- dim(condensed_raw_venn_tree)[1]
branchings <- reduction$Branchings
first_taxon <- raw_venn_tree[1,]
node_ancestors <- as.numeric(rownames(condensed_raw_venn_tree))
for (n in 1:otus)	{
	u <- 1
	while (u < branchings[n])	{
		first_taxon[n] <- ancestral[first_taxon[n]]
		u <- u+1
		}
	}
first_taxon <- c(first_taxon,node_ancestors)
branchings <- c(branchings,rep(0,unique_nodes))
# now, take duplicate nodes and reduce them, adding to their branch lengths
for (u in 2:unique_nodes)	{
	htu <- otus + u
	row.is.a.match <- apply(raw_venn_tree, 1, identical, condensed_raw_venn_tree[u,]) 
	match.idx <- sort(which(row.is.a.match))
	branchings[htu] <- sum(row.is.a.match)
	}
#	first_taxon <- c(first_taxon,all_ancestors[match.idx[1]])
#	if (is.na(match(node_ancestors[htu],all_ancestors)))	{
#		first_taxon[htu] <- all_ancestors[match.idx[1]]
#		}	else	{
#		first_taxon[htu] <- condensed_raw_venn_tree[u,1]
#		branchings[condensed_raw_venn_tree[u,1]] <- 0
#		}
#	u <- u+1
#	first_taxon
#	branchings
	# add something here to track down original ancestor using all_ancestors vector
#	which(condensed_raw_venn_tree[u,1]==raw_venn_tree[,1],arr.ind=TRUE)
#	i<-1:unique_nodes
#	which(condensed_raw_venn_tree[u,],raw_venn_tree[i,])
#	}

# there might be ancestral species sampled: if so, then they have
#	branch durations of zero
all_ancestors <- accersi_all_ancestors_for_sampled_taxa(sampled,ancestral)
sampled_ancestors <- c()
#if (!is.na(match(otus,all_ancestors)))	{
## PICK UP HERE!!!!
for (s in 1:otus)	{
	if (!is.na(match(sampled[s],all_ancestors)))	{
		sampled_ancestors <- c(sampled_ancestors,sampled[s])
		htu <- otus+match(sampled[s],condensed_raw_venn_tree[,1])
		### ancestral branches already have at least one
		###		add any additional branches
		branchings[htu] <- branchings[htu] + (branchings[s] - 1)
		branchings[s] <- 0
		}
	}

# now, change all of the sampled taxon numbers to 1…S
venn_tree <- match(condensed_raw_venn_tree[1,],raw_venn_tree[1,])
for (n in 2:unique_nodes)	{
	xx <- match(condensed_raw_venn_tree[n,],raw_venn_tree[1,])
	xx[is.na(xx)] <- 0
	venn_tree <- rbind(venn_tree,xx)
	}
# name rows after simulated founder
#rownames(venn_tree) <- first_taxon[(otus+1):length(first_taxon)]
rownames(venn_tree) <- rownames(condensed_raw_venn_tree)

output <- list(venn_tree,branchings,first_taxon)
names(output) <- c("Venn_Tree","Branchings","First_Taxon")
return(output)
}

# routine to eliminate nodes with only one sampled descendant
remove_singletons_from_venn_tree <- function(venn_tree)	{
sampled <- sort(venn_tree[1,])
init_nodes <- dim(venn_tree)[1]
otus <- length(sampled)
branchings <- rep(1,otus)
fs <- c()
for (n in 1:init_nodes)
	if (sum(venn_tree[n,]>0)==1)
		fs <- c(fs,n)
# add to branch lengths of species with singleton representations
for (sn in 1:length(fs))	{
	n <- fs[sn]
	spc <- match(venn_tree[n,1],sampled)
	branchings[spc] <- branchings[spc]+1
	}
# list retained nodes
ret_nodes <- (1:init_nodes)[!(1:init_nodes) %in% fs]
red_venn_tree <- venn_tree[ret_nodes,]
output <- list(red_venn_tree,branchings)
names(output) <- c("Condensed_Venn_Tree","Branchings")
return(output)
}

#mine <- 1:6 
#table.combos <- matrix(data = 1:12, nrow = 10, ncol = 6, byrow=T) 
#row.is.a.match <- apply(table.combos, 1, identical, mine) 
#match.idx <- which(row.is.a.match) 
#total.matches <- sum(row.is.a.match) 

# routine to get basic venn tree from initial simulation
accersi_raw_venn_tree_from_ancestor_list <- function(sampled,ancestral)	{
notu <- length(sampled)
all_ancestors <- accersi_all_ancestors_for_sampled_taxa(sampled,ancestral)
all_ancestors_per_otu <- list_all_ancestors_for_all_sampled_taxa(sampled,ancestral)
raw_venn_tree <- matrix(0,length(all_ancestors),notu)
all_desc <-rep(0,length(all_ancestors))
for (i in 1:notu)	{
	spc <- sampled[i]
	nodes <- match(all_ancestors_per_otu[[i]],all_ancestors)
	if (!is.na(match(spc,all_ancestors)))	nodes <- c(nodes,match(spc,all_ancestors))
	all_desc[nodes] <- all_desc[nodes]+1
	for (n in 1:length(nodes))	{
#		print(c(n,nodes[n],all_desc[nodes[n]]))
		raw_venn_tree[nodes[n],all_desc[nodes[n]]] <- spc
		}
#	raw_venn_tree[,1:3]
#	i <- i+1
	}
rNodes <- dim(raw_venn_tree)[1]
rownames(raw_venn_tree) <- all_ancestors
while (identical(raw_venn_tree[1,],raw_venn_tree[2,]))	{
	raw_venn_tree <- raw_venn_tree[2:rNodes,]
	rNodes <- rNodes - 1
	}
return(raw_venn_tree)	# it is fine at this point: see what happens to it....
}

# routine to list all simulated ancestors for "sampled" simulated taxa
accersi_all_ancestors_for_sampled_taxa <- function(sampled,ancestral)	{
otu <- length(sampled)
added <- c()
for (i in 1:otu)	{
	spc <- sampled[i]
	added <- sort(unique(c(added,accersi_all_ancestors_for_taxon(spc,ancestral))),decreasing=FALSE)
	}
return(added)
}

# routine to list all simulated ancestors for each "sampled" simulated taxa
list_all_ancestors_for_all_sampled_taxa <- function(sampled,ancestral)	{
otu <- length(sampled)
all_anc <- list()
for (i in 1:otu)	{
	spc <- sampled[i]
	all_anc[[i]] <- sort(accersi_all_ancestors_for_taxon(spc,ancestral),decreasing=FALSE)
	}
return(all_anc )
}

# routine to list all simulated ancestors for a "sampled" simulated taxa
accersi_all_ancestors_for_taxon <- function(spc,ancestral)	{
anc <- ancestral[spc]
while (min(anc)>=1)	{
	anc <- c(anc,ancestral[min(anc)])
	}
return(anc)
}

# simulate character evolution up to N steps & tally compatibility
#		2017-02-22: use beta distribution to "smooth" P[compatibility|steps]!!!
#		2017-08-10: make sure initial compatibility is maximum or close to it.....
evolve_compatibility_over_N_changes <- function(N,init_chmatrix,venn_tree,branchings,nchars,states,types,hidden_reversals=TRUE,UNKNOWN,INAP,repl=1)	{
# N: maximum number of steps
# venn_tree: a tree giving all of the observed descendants (direct & indirect)
# branchings: length of each branch
# nchars: number of characters
# states: number of states for each character
# types: types for each character (0: unordered; 1: ordered)
# hidden_reversals: if true, then character can change 2+ times per sampled branch
# UNKNOWN: numeric code for "?"
# INAP: numeric code for "-"
# repl: replication number
notu <- dim(venn_tree)[2]	# number of observed taxa in tree
nodes <- dim(venn_tree)[1]
#simchmatrix <- matrix(0,notu,nchars)
simchmatrix <- init_chmatrix
simchmatrix[simchmatrix>0] <- 0	# added 2017-08-10
# get those branches where there can be change (sampled ancestors are zero, with nodal branch ≥ 1)
unique_ab <- (1:length(branchings))[branchings>0]
#ab <- vector(length=branches)		# available branches, with branches with 2+ species listed 2+ times

# each branch gets one additional representative per unsampled ancestor.
# 		With perfect sampling, all taxa are entered once.
ab <- c()
for (b in 1:length(unique_ab))	{
	br <- unique_ab[b]
	ab <- c(ab,rep(br,branchings[br]))
#	if (b==1)	{
#		# ab: availabe branches
#		ab <- rep(br,branchings[br])
#		}	else	{
#		ab <- c(ab,rep(br,branchings[br]))
#		}
	}
tb <- sum(branchings)
#tb <- length(ab)			# total branches; = branches above, so why both?
fb <- length(unique_ab)		# free branches, with each branch that can change listed only once
bcc <- matrix(0,nchars,fb)	# branch changes per character
st_rich <- matrix(0,nchars,max(states))	# richness of each state

mx_ch <- count_scored_otu_per_character(chmatrix=simchmatrix)

rich <- vector(length=nodes)
for (n in 1:nodes)  rich[n] <- sum(venn_tree[n,]>0)
# for (n in 1:nodes)  rich[n] <- length(subset(venn_tree[n,],venn_tree[n,]>0))
pabm <- get_possible_branches_for_all_characters(simchmatrix,branchings,ab,venn_tree) # possible available branches matrix
pcc <- vector(length=nchars)  # possible character changes
for (ch in 1:nchars)  {
	# scramble order in which branches are sampled
	if (hidden_reversals==TRUE)	{
		x <- permute(subset(pabm[ch,],pabm[ch,]>0))
		}	else	{
		x <- unique(permute(subset(pabm[ch,],pabm[ch,]>0)))
		}
	pcc[ch] <- length(x)
	pabm[ch,1:pcc[ch]] <- x
	}
	# this routine should be unnecessary, and it crashes sometimes.
#	if (length(x) < tb)	{
#	if (pcc[ch] < tb)	{
#		y <- 1+pcc[ch]
#		pabm[ch,y:tb] <- 0	# problem here sometimes...
#		pabm[ch,(pcc[ch]+1):tb] <- 0	# problem here sometimes...
#		}
#	if (pcc[ch]<length(ab)	for (c in (pcc[ch]+1):length(ab))   pabm[ch,c] <- 0
#	}

# first, make sure that all states appear
branch_changes <- vector(length=notu+nodes)
char_changes <- vector(length=nchars)
simcompat <- vector(length=N)
for (ch in 1:nchars)	{
	if (states[ch]>=2)	{
    	use_br <- sort(pabm[ch,1:(states[ch]-1)],decreasing=TRUE)   # do high nodes first, then low nodes, then otus
		for (c in 1:(states[ch]-1))	{
			br <- use_br[c]
			if (br<=notu)	{
				simchmatrix[br,ch] <- c
				} else {
				ndd <- br-notu
				for (r in 1:rich[ndd])	{
					s <- venn_tree[ndd,r]
					if (simchmatrix[s,ch]==0)	simchmatrix[s,ch] <- c
					}
				}
			char_changes[ch] <- 1+char_changes[ch]
			branch_changes[br] <- branch_changes[br]+1
			}	# add state c to matrix
		}	# only do variant characters
	}

simstates <- vector(length=nchars)
for (ch in 1:nchars)	simstates[ch] <- sum(unique(simchmatrix[,ch])>=0)

# makes sure that all branches have 1+ change
# branch_changes[unique(ab)]
#  PUT THIS ROUTINE IN OTHER MODULES!
if (sum(branch_changes[unique_ab]==0)>0)	{
	needy <- unique(ab)[branch_changes[unique(ab)]==0]
	nb <- length(needy)
	for (b in 1:nb)	{
		# get first branch in need of a change
		br <- needy[b]
		# target characters where it would change earliest
		candidates <- which(pabm==br,arr.ind=TRUE)
		cn <- 1
		ch <- candidates[cn,1]
		while (states[ch]<2
			   || char_changes[ch]>=mx_ch[ch]
			   || char_changes[ch]>=candidates[cn,2]
			   || is.na(match(br,pabm[ch,])))	{
			cn <- cn+1
			ch <- candidates[cn,1]
			}	# make sure that this is an appropriate character
		# determine shifts
		if (states[ch]==2)	{
			dstates <- c(1,0)
			} else if (types[ch]==0)	{
			dstates <- scramble_multistates(states[ch])	
			} 
		# routine for tips
		if (br<=notu)	{
			if (simchmatrix[br,ch]!=UNKNOWN && simchmatrix[br,ch]!=INAP)	{
				cc <- simchmatrix[br,ch]
				simchmatrix[br,ch] <- dstates[cc+1]
				} else {
				ndg <- br-notu
				for (d in 1:length(subset(venn_tree[ndg,],venn_tree[ndg,]>0)))	{
					sp <- venn_tree[ndg,d]
					if (simchmatrix[sp,ch]!=UNKNOWN && simchmatrix[sp,ch]!=INAP)	simchmatrix[sp,ch] <- dstates[1+simchmatrix[sp,ch]]
					}	
				}
			#change_character_on_branch(ch,branchings[b])
			}
		branch_changes[br] <- 1
		char_changes[ch] <- char_changes[ch]+1
		# now, flip over branches
		pabm[ch,candidates[cn,2]] <- pabm[ch,char_changes[ch]]
		pabm[ch,char_changes[ch]] <- br
		}	# end case where branch needed changes
	}
#third, tally compatibility at this point
delta <- sum(char_changes)
simcompmat <- compatibility_matrix(simchmatrix,states,types,UNKNOWN,INAP)
character_compats <- rowSums(simcompmat)-1
simcompat[delta] <- sum(character_compats)/2
#simcompat[delta] <- total_compatibility(simchmatrix,nchars,states,types,notu,UNKNOWN,INAP)
steps_v_compat <- c(delta,simcompat[delta])
if (repl>0)	print(c(repl,delta,simcompat[delta]))

counter <- 0
for (d in (delta+1):N)	{
	counter <- counter+1
	# add the br != 0 check!!!
	br <- 0
	while (br==0)	{
		ch <- ceiling(nchars*runif(1))
		while (ch>nchars || char_changes[ch]>=mx_ch[ch])	ch <- ceiling(nchars*runif(1))
		c <- char_changes[ch]+1
		br <- pabm[ch,c]
		}
	
	if (states[ch]==2)	{
		dstates <- c(1,0)       # this is a transition matrix: 0->1, 1->0
		} else if (types[ch]==0)	{
		dstates <- scramble_multistates(states[ch])	# this transition matrix will have either 0->1+1->2+2->0 or 0->2+1->0+2->1 for a 3 state character
		} 
#	prior_compat <- compatibility_of_a_character <- (ch1,simchmatrix,nchars,states,types,notu,UNKNOWN,INAP)
	prior_compat <- simcompmat[ch,]
	
	if (br<=notu)	{
		if (simchmatrix[br,ch]!=UNKNOWN && simchmatrix[br,ch]!=INAP)	{
			simchmatrix[br,ch] <- dstates[1+simchmatrix[br,ch]]
			} else {
			ndh <- br-notu
			for (d in 1:length(subset(venn_tree[ndh,],venn_tree[ndh,]>0)))	{
				sp <- venn_tree[ndh,d]
				if (simchmatrix[sp,ch]!=UNKNOWN && simchmatrix[sp,ch]!=INAP)	simchmatrix[sp,ch] <- dstates[1+simchmatrix[sp,ch]]
				}	
			}
		#change_character_on_branch(ch,branchings[b])
		}
	# new vector of compatibilities
	new_compat <- compatibility_of_a_character(ch,simchmatrix,states,types,UNKNOWN,INAP)

	# update compatibility matrix
	# tally new matrix compatibility
	#simcompat[d] <- (sum(simcompmat)-nchars)/2
	dcompat <- sum(new_compat - prior_compat)
	simcompat[d] <- simcompat[d-1]+dcompat
#	simcompat[d] <- simcompat[d-1]-(sum(prior_compat)-sum(new_compat))
	if (repl>0 && counter%%10==0)	print(c(repl,d,simcompat[d]))
	steps_v_compat <- rbind(steps_v_compat,c(d,simcompat[d]))
	branch_changes[br] <- branch_changes[br]+1
	char_changes[ch] <- char_changes[ch]+1
	if (sum(prior_compat==new_compat)!=nchars)	{
		simcompmat[ch,] <- new_compat
		simcompmat[,ch] <- new_compat
		for (c2 in 1:nchars)	{
			if (new_compat[c2]!=prior_compat[c2])	{
				if (new_compat[c2]==0)	{
					character_compats[c2] <- character_compats[c2]-1
					} else if (new_compat[c2]==1)	{
					character_compats[c2] <- character_compats[c2]+1
					}
				}	# case of mismatch
			}	# modify compatibility matrix of cha}racters affected by change
		}
	}
#output <- list(d,simcompmat)
rownames(steps_v_compat) <- rep("",dim(steps_v_compat)[1])
colnames(steps_v_compat) <- c("Steps","Compat")
return (steps_v_compat)
}

# simulate character evolution up to a certain compatibility
# returns # steps required.  This will be a fraction if compatibility is passed
#evolve_to_particular_compatibilty <- function(ttl_compat,venn_tree,branchings,nchars,states,types,maxsteps,UNKNOWN,INAP,repl=-1)	{
evolve_to_particular_compatibilty_old <- function(init_chmatrix,ttl_compat,venn_tree,branchings,states,types,UNKNOWN,INAP,repl=-1)	{
# init_chmatrix: initial character matrix to be emulated
# ttl_compat: compatibility to which to evolve
# venn_tree: a tree giving all of the observed descendants (direct & indirect)
# branchings: length of each branch
# nchars: number of characters
# states: number of states for each character
# types: types for each character (0: unordered; 1: ordered)
# maxsteps: maximum number of steps to do
# UNKNOWN: numeric code for "?"
# INAP: numeric code for "-"
# repl: replication number
notu <- dim(venn_tree)[2]	# number of observed taxa in tree
nodes <- dim(venn_tree)[1]
nchars <- ncol(init_chmatrix)
simchmatrix <- init_chmatrix
simchmatrix[simchmatrix>0] <- 0	# added 2017-08-10

unique_ab <- (1:length(branchings))[branchings>0]
#ab <- vector(length=branches)		# available branches, with branches with 2+ species listed 2+ times
branch_changes <- vector(length=notu+nodes)
char_changes <- vector(length=nchars)

ab <- c()
for (b in 1:length(unique_ab))	{
	br <- unique_ab[b]
	ab <- c(ab,rep(br,branchings[br]))
	}
tb <- length(ab)		# total branches; = branches above, so why both?
fb <- length(unique_ab)		# free branches, with each branch that can change listed only once
bcc <- matrix(0,nchars,fb)	# branch changes per character
st_rich <- matrix(0,nchars,max(states))
mx_ch <- vector(length=nchars)
for (ch in 1:nchars)	{
	for (s in 1:notu)	{
		if (chmatrix[s,ch]!=UNKNOWN && chmatrix[s,ch]!=INAP)	st_rich[ch,1] <- st_rich[ch,1]+1
		}
	if (states[ch]>=2)	mx_ch[ch] <- st_rich[ch,1]
	}

rich <- vector(length=nodes)
for (n in 1:nodes)  rich[n] <- sum(venn_tree[n,]>0)
#for (n in 1:nodes)  rich[n] <- length(subset(venn_tree[n,],venn_tree[n,]>0))
pabm <- get_possible_branches_for_all_characters(simchmatrix,branchings,ab,venn_tree) # possible available branches matrix
pcc <- vector(length=nchars)  # possible character changes

for (ch in 1:nchars)  {
	x <- unique(permute(subset(pabm[ch,],pabm[ch,]>0)))
	pcc[ch] <- length(x)
	pabm[ch,1:pcc[ch]] <- x
	# this causes problems sometimes and should not be necessary
#	if (pcc[ch] < tb)	pabm[ch,(pcc[ch]+1):tb] <- 0
	}

# first, make sure that all states appear
for (ch in 1:nchars)	{
	if (states[ch]>=2)	{
    use_br <- sort(pabm[ch,1:(states[ch]-1)],decreasing=TRUE)   # do high nodes first, then low nodes, then otus
		for (c in 1:(states[ch]-1))	{
			br <- use_br[c]
			if (br<=notu)	{
				simchmatrix[br,ch] <- c
				} else {
				ndd <- br-notu
				for (r in 1:rich[ndd])	{
					s <- venn_tree[ndd,r]
					if (simchmatrix[s,ch]==0)	simchmatrix[s,ch] <- c
					}
				}
			char_changes[ch] <- 1+char_changes[ch]
			branch_changes[br] <- branch_changes[br]+1
			}	# add state c to matrix
		}	# only do variant characters
	}

simstates <- vector(length=nchars)
for (ch in 1:nchars)	simstates[ch] <- sum(unique(simchmatrix[,ch])>=0)

# makes sure that all branches have 1+ change
#  PUT THIS ROUTINE IN OTHER MODULES!
needy <- unique(ab)[branch_changes[unique(ab)]==0]
nb <- length(needy)
b <- 1
while (b < nb)	{
	# get first branch in need of a change
	br <- needy[b]
	# target characters where it would change earliest
	candidates <- which(pabm==br,arr.ind=TRUE)
	cn <- 1
	ch <- candidates[cn,1]
	while (states[ch]<2
		   || char_changes[ch]>=mx_ch[ch]
		   || char_changes[ch]>=candidates[cn,2]
		   || is.na(match(br,pabm[ch,])))	{
		cn <- cn+1
		ch <- candidates[cn,1]
		}	# make sure that this is an appropriate character
	# determine shifts
	if (states[ch]==2)	{
		dstates <- c(1,0)
		} else if (types[ch]==0)	{
		dstates <- scramble_multistates(states[ch])	
		} 
	# routine for tips
	if (br<=notu)	{
		if (simchmatrix[br,ch]!=UNKNOWN && simchmatrix[br,ch]!=INAP)	{
			cc <- simchmatrix[br,ch]
			simchmatrix[br,ch] <- dstates[cc+1]
			} else {
			ndg <- br-notu
			for (d in 1:length(subset(venn_tree[ndg,],venn_tree[ndg,]>0)))	{
				sp <- venn_tree[ndg,d]
				if (simchmatrix[sp,ch]!=UNKNOWN && simchmatrix[sp,ch]!=INAP)	simchmatrix[sp,ch] <- dstates[1+simchmatrix[sp,ch]]
				}	
			}
		#change_character_on_branch(ch,branchings[b])
		}
	branch_changes[br] <- 1
	char_changes[ch] <- char_changes[ch]+1
	# now, flip over branches
	pabm[ch,candidates[cn,2]] <- pabm[ch,char_changes[ch]]
	pabm[ch,char_changes[ch]] <- br
	b <- b+1
	}	# end case where branch needed changes

#third, tally compatibility at this point
delta <- sum(char_changes)
simcompmat <- compatibility_matrix(simchmatrix,states,types,UNKNOWN,INAP)
simcharacter_compats <- rowSums(simcompmat)-1
simcompat <- vector(length=delta)
simcompat[delta] <- sum(simcharacter_compats)/2
#simcompat[delta] <- total_compatibility(simchmatrix,nchars,states,types,notu,UNKNOWN,INAP)
steps_v_compat <- c(delta,simcompat[delta])
if (repl>0)	print(c(repl,delta,simcompat[delta]))

d <- delta
counter <- 0
simcompat_test <- simcompat	# for debugging.
while (simcompat[d]>ttl_compat)	{
	counter <- counter+1
	d <- d+1
	br <- 0
	br_counter <- 0
	while (br==0 && br_counter<100)	{
		ch <- ceiling(nchars*runif(1))
		while (ch>nchars || char_changes[ch]>=mx_ch[ch])	ch <- ceiling(nchars*runif(1))
		c <- char_changes[ch]+1
		br <- pabm[ch,c]
		br_counter <- br_counter+1
		}
	if (br_counter <= 99)	{

		prior_char <- simchmatrix[,ch]
		if (states[ch]==2)	{
			dstates <- c(1,0)       # this is a transition matrix: 0->1, 1->0
			} else if (types[ch]==0)	{
			dstates <- scramble_multistates(states[ch])	# this transition matrix will have either 0->1+1->2+2->0 or 0->2+1->0+2->1 for a 3 state character
			} 
#	prior_compat <- compatibility_of_a_character <- (ch1,simchmatrix,nchars,states,types,notu,UNKNOWN,INAP)
#		prior_compat <- simcompmat[ch,]
		prior_compat <- simcompmat
		if (br<=notu)	{
			if (simchmatrix[br,ch]!=UNKNOWN && simchmatrix[br,ch]!=INAP)	{
				simchmatrix[br,ch] <- dstates[1+simchmatrix[br,ch]]
				} else {
				ndh <- br-notu
				for (d in 1:length(subset(venn_tree[ndh,],venn_tree[ndh,]>0)))	{
					sp <- venn_tree[ndh,d]
					if (simchmatrix[sp,ch]!=UNKNOWN && simchmatrix[sp,ch]!=INAP)	simchmatrix[sp,ch] <- dstates[1+simchmatrix[sp,ch]]
					}	
				}
		#change_character_on_branch(ch,branchings[b])
			}
	# new vector of compatibilities
		new_compat <- compatibility_of_a_character(ch,simchmatrix,states,types,UNKNOWN,INAP)
	#	simcompmat <- compatibility_matrix(simchmatrix,states,types,UNKNOWN,INAP)
		simcompmat[,ch] <- simcompmat[ch,] <- new_compat
		simcharacter_compats <- rowSums(simcompmat)-1
		simcompat <- c(simcompat,sum(simcharacter_compats)/2)
		simcompat_test[d] <- total_compatibility(chmatrix=simchmatrix,states,types,UNKNOWN,INAP)
		# update compatibility matrix
		# tally new matrix compatibility
		#simcompat[d] <- (sum(simcompmat)-nchars)/2
		newmatcomp <- simcompat[d-1]-(sum(prior_compat)-sum(new_compat))
	#	simcompat <- c(simcompat,simcompat[d-1]-(sum(prior_compat)-sum(new_compat)))
		if (repl>0 && counter%%10==repl)	print(c(d,simcompat[d]))
		steps_v_compat <- rbind(steps_v_compat,c(d,simcompat[d]))
		if (simcompat[d]>ttl_compat)	{
		# update character compatibilities that might have been altered
			branch_changes[br] <- branch_changes[br]+1
			char_changes[ch] <- char_changes[ch]+1
			simcompmat[ch,] <- simcompmat[,ch] <- new_compat
	#		for (c2 in 1:nchars)	{
	#			if (new_compat[c2]!=prior_compat[c2])	{
	#				if (new_compat[c2]==0)	{
	#					character_compats[c2] <- character_compats[c2]-1
	#					} else if (new_compat[c2]==1)	{
	#					character_compats[c2] <- character_compats[c2]+1
	#					}
	#				}	# case of mismatch
	#			}
			} else	if (simcompat[d]<ttl_compat)	{
				
			}
		}
	xxx <- cbind(simcompat,simcompat_test)
	print(c(xxx[d,],ttl_compat))
	}

	#simcompat[d-1]
	#simcompat[d]
	#d+(ttl_compat-simcompat[d])/(simcompat[d-1]-ttl_compat)
	### START HERE!!!!
return (d+(ttl_compat-simcompat[d])/(simcompat[d-1]-ttl_compat))
}

# simulate character evolution up to a certain compatibility
evolve_up_to_compatibility <- function(init_chmatrix,ttl_compat,venn_tree,branchings,nchars,states,types,UNKNOWN,INAP,repl=-1)	{
# ttl_compat: compatibility to which to evolve
# venn_tree: a tree giving all of the observed descendants (direct & indirect)
# branchings: length of each branch
# nchars: number of characters
# states: number of states for each character
# types: types for each character (0: unordered; 1: ordered)
# maxsteps: maximum number of steps to do
# UNKNOWN: numeric code for "?"
# INAP: numeric code for "-"
# repl: replication number
notu <- dim(venn_tree)[2]	# number of observed taxa in tree
nodes <- dim(venn_tree)[1]
simchmatrix <- init_chmatrix
for (ch in 1:nchars)	simchmatrix[simchmatrix[,ch]>=0,ch] <- 0
unique_ab <- (1:length(branchings))[branchings>0]
#ab <- vector(length=branches)		# available branches, with branches with 2+ species listed 2+ times
branch_changes <- vector(length=notu+nodes)
char_changes <- vector(length=nchars)
#if (maxsteps==-1)	{
#	taxa_scored <- count_scored_otu_per_character(init_chmatrix)
#	maxsteps_per_char <- accersi_maximum_parsimony_steps_per_character(init_chmatrix,states)
#	maxsteps <- ceiling((sum(maxsteps_per_char)+sum(taxa_scored))/2)
#	}
#simcompat <- vector(length=maxsteps)
ab <- c()
for (b in 1:length(unique_ab))	{
	br <- unique_ab[b]
	ab <- c(ab,rep(br,branchings[br]))
#	if (b==1)	{
#		ab <- rep(br,branchings[br])
#		}	else	{
#		ab <- c(ab,rep(br,branchings[br]))
#		}
	}
tb <- length(ab)		# total branches; = branches above, so why both?
fb <- length(unique_ab)		# free branches, with each branch that can change listed only once
bcc <- matrix(0,nchars,fb)	# branch changes per character
st_rich <- matrix(0,nchars,max(states))
mx_ch <- vector(length=nchars)
for (ch in 1:nchars)	{
	for (s in 1:notu)	{
		if (chmatrix[s,ch]!=UNKNOWN && chmatrix[s,ch]!=INAP)	st_rich[ch,1] <- st_rich[ch,1]+1
		}
	if (states[ch]>=2)	mx_ch[ch] <- st_rich[ch,1]
	}

rich <- vector(length=nodes)
for (n in 1:nodes)  rich[n] <- sum(venn_tree[n,]>0)
#for (n in 1:nodes)  rich[n] <- length(subset(venn_tree[n,],venn_tree[n,]>0))
pabm <- get_possible_branches_for_all_characters(simchmatrix,branchings,ab,venn_tree) # possible available branches matrix
pcc <- vector(length=nchars)  # possible character changes

for (ch in 1:nchars)  {
	x <- unique(permute(subset(pabm[ch,],pabm[ch,]>0)))
	pcc[ch] <- length(x)
	for (c in 1:pcc[ch])  		pabm[ch,c] <- x[c]
	if (pcc[ch]<length(ab))   for (c in (pcc[ch]+1):length(ab))   pabm[ch,c] <- 0
	}
#rm(ab,list=subset(bcc[1,],bcc[1,]>0))

#delta <- sum(char_changes)
#simcompatmat <- compatibility_matrix(simchmatrix,states,types,UNKNOWN,INAP)
simcompatmat <- matrix(1,nchars,nchars)
prior_compat <- character_compats <- nchars-1
#simcompat[delta] <- sum(character_compats)/2
d <- 0
simcompat <- c()
# first, make sure that all states appear
for (ch in 1:nchars)	{
	if (states[ch]>=2)	{
		prior_compat <- simcompatmat[ch,]
    	use_br <- sort(pabm[ch,1:(states[ch]-1)],decreasing=TRUE)   # do high nodes first, then low nodes, then otus
		for (c in 1:(states[ch]-1))	{
			br <- use_br[c]
			if (br<=notu)	{
				simchmatrix[br,ch] <- c
				} else {
				ndd <- br-notu
				for (r in 1:rich[ndd])	{
					s <- venn_tree[ndd,r]
					if (simchmatrix[s,ch]==0)	simchmatrix[s,ch] <- c
					}
				}
			char_changes[ch] <- 1+char_changes[ch]
			branch_changes[br] <- branch_changes[br]+1
			d <- d+1
			new_compat <- compatibility_of_a_character(ch,simchmatrix,states,types,UNKNOWN,INAP)
			simcompatmat[ch,] <- simcompatmat[,ch] <- new_compat
			simcompat <- c(simcompat,sum(simcompatmat[lower.tri(simcompatmat)]))
#			if (d>1)	{
#				simcompat <- c(simcompat,simcompat[d-1]-(sum(prior_compat)-sum(new_compat)))
#				}	else	simcompat <- ((nchars^2)-nchars)/2
			}	# add state c to matrix
		}	# only do variant characters
	}

simstates <- vector(length=nchars)
for (ch in 1:nchars)	simstates[ch] <- sum(unique(simchmatrix[,ch])>=0)

# second, make sure that all branches have change
counter <- 0
if (min(branch_changes[unique_ab])==0)	{
	needy_brs <- unique_ab[branch_changes[unique_ab]==0]
	for (b in 1:length(needy_brs))	{
		br <- needy_brs[b]	### make sure that is in all such routines!!!!
#		if (branch_changes[br]==0)	{
		counter <- counter+1
		d <- d+1
		ch <- ceiling(nchars*runif(1))
		while (states[ch]<2 || char_changes[ch]>=mx_ch[ch])	ch <- ceiling(nchars*runif(1))
		prior_compat <- simcompatmat[ch,]
		if (states[ch]==2)	{
			dstates <- c(1,0)
			} else if (types[ch]==0)	{
			dstates <- scramble_multistates(states[ch])	
			} 
		if (br<=notu)	{
			if (simchmatrix[br,ch]!=UNKNOWN && simchmatrix[br,ch]!=INAP)	{
				cc <- simchmatrix[br,ch]
				simchmatrix[br,ch] <- dstates[cc+1]
				} else {
				ndg <- br-notu
				for (d in 1:length(subset(venn_tree[ndg,],venn_tree[ndg,]>0)))	{
					sp <- venn_tree[ndg,d]
					if (simchmatrix[sp,ch]!=UNKNOWN && simchmatrix[sp,ch]!=INAP)	simchmatrix[sp,ch] <- dstates[1+simchmatrix[sp,ch]]
					}	
				}
			#change_character_on_branch(ch,branchings[b])
			}
		branch_changes[br] <- branch_changes[br]+1
		char_changes[ch] <- char_changes[ch]+1
		new_compat <- compatibility_of_a_character(ch,simchmatrix,states,types,UNKNOWN,INAP)
		simcompat[d] <- simcompat[d-1]-(sum(prior_compat)-sum(new_compat))
	#		bcc[ch,char_changes[ch]] <- br
		}	# end case where branch needed changes
	}

#third, tally compatibility at this point
#simcompat[delta] <- total_compatibility(simchmatrix,nchars,states,types,notu,UNKNOWN,INAP)
steps_v_compat <- c(delta,simcompat[delta])
if (repl>0)	print(c(repl,delta,simcompat[delta]))
if (simcompat[d]<ttl_compat)	{
	### find where we overshot
	
	}	else	{
	counter <- 1
	while (simcompat[d]>ttl_compat)	{
	### while loop added 2017-06-03
		br <- 0
		while (br==0)	{
			ch <- ceiling(nchars*runif(1))
			while (ch>nchars || char_changes[ch]>=mx_ch[ch])	ch <- ceiling(nchars*runif(1))
			c <- char_changes[ch]+1
			br <- pabm[ch,c]
			}
		prior_char <- simchmatrix[,ch]
		prior_compat <- simcompatmat[ch,]
	
		if (states[ch]==2)	{
			dstates <- c(1,0)       # this is a transition matrix: 0->1, 1->0
			} else if (types[ch]==0)	{
			dstates <- scramble_multistates(states[ch])	# this transition matrix will have either 0->1+1->2+2->0 or 0->2+1->0+2->1 for a 3 state character
			} 
	#	prior_compat <- compatibility_of_a_character <- (ch1,simchmatrix,nchars,states,types,notu,UNKNOWN,INAP)
		if (br<=notu)	{
			if (simchmatrix[br,ch]!=UNKNOWN && simchmatrix[br,ch]!=INAP)	{
				simchmatrix[br,ch] <- dstates[1+simchmatrix[br,ch]]
				} else {
				ndh <- br-notu
				for (d in 1:length(subset(venn_tree[ndh,],venn_tree[ndh,]>0)))	{
					sp <- venn_tree[ndh,d]
					if (simchmatrix[sp,ch]!=UNKNOWN && simchmatrix[sp,ch]!=INAP)	simchmatrix[sp,ch] <- dstates[1+simchmatrix[sp,ch]]
					}	
				}
			#change_character_on_branch(ch,branchings[b])
			}
		# new vector of compatibilities
#		ccc <- compatibility_matrix(simchmatrix,states,types,UNKNOWN,INAP)
		new_compat <- compatibility_of_a_character(ch,simchmatrix,states,types,UNKNOWN,INAP)
		# update compatibility matrix
		# tally new matrix compatibility
		#simcompat[d] <- (sum(simcompatmat)-nchars)/2
		d <- d+1
#		simcompat <- c(simcompat,simcompat[d-1]-(sum(prior_compat)-sum(new_compat)))
		simcompatmat[ch,] <- new_compat
		simcompatmat[,ch] <- new_compat
		simcompat <- c(simcompat,sum(simcompatmat[lower.tri(simcompatmat)]))
		if (repl>0 && counter%%10==repl)	print(c(d,simcompat[d]))
		counter <- counter+1
		steps_v_compat <- rbind(steps_v_compat,c(d,simcompat[d]))
		# update character compatibilities that might have been altered
		branch_changes[br] <- branch_changes[br]+1
		char_changes[ch] <- char_changes[ch]+1
		for (c2 in 1:nchars)	{
			if (new_compat[c2]!=prior_compat[c2])	{
				if (new_compat[c2]==0)	{
					character_compats[c2] <- character_compats[c2]-1
					} else if (new_compat[c2]==1)	{
					character_compats[c2] <- character_compats[c2]+1
					}
				}	# case of mismatch
			}
		}
	#output <- list(d,simcompatmat)
	}

return ((d-1)+abs(simcompat[d-1]-ttl_compat)/abs(simcompat[d-1]-simcompat[d]))
}
	
# simulate character evolution up to a certain compatibility
emulate_observed_compatibility <- function(ttl_compat,init_chmatrix,venn_tree,branchings,nchars,states,types,maxsteps,UNKNOWN,INAP,repl=-1)	{
# ttl_compat: compatibility to which to evolve
# venn_tree: a tree giving all of the observed descendants (direct & indirect)
# branchings: length of each branch
# nchars: number of characters
# states: number of states for each character
# types: types for each character (0: unordered; 1: ordered)
# maxsteps: maximum number of steps to do
# UNKNOWN: numeric code for "?"
# INAP: numeric code for "-"
# repl: replication number
notu <- dim(venn_tree)[2]	# number of observed taxa in tree
nodes <- dim(venn_tree)[1]
simchmatrix <- init_chmatrix
simchmatrix[simchmatrix>0] <- 0	# added 2017-08-10
unique_ab <- (1:length(branchings))[branchings>0]
#ab <- vector(length=branches)		# available branches, with branches with 2+ species listed 2+ times
branch_changes <- vector(length=notu+nodes)
char_changes <- vector(length=nchars)
simcompat <- vector(length=maxsteps)

for (b in 1:length(unique_ab))	{
	br <- unique_ab[b]
	if (b==1)	{
		ab <- rep(br,branchings[br])
		}	else	{
		ab <- c(ab,rep(br,branchings[br]))
		}
	}
tb <- length(ab)		# total branches; = branches above, so why both?
fb <- length(unique_ab)		# free branches, with each branch that can change listed only once
bcc <- matrix(0,nchars,fb)	# branch changes per character
st_rich <- matrix(0,nchars,max(states))
mx_ch <- vector(length=nchars)
for (ch in 1:nchars)	{
	for (s in 1:notu)	{
		if (chmatrix[s,ch]!=UNKNOWN && chmatrix[s,ch]!=INAP)	st_rich[ch,1] <- st_rich[ch,1]+1
		}
	if (states[ch]>=2)	mx_ch[ch] <- st_rich[ch,1]
	}

rich <- vector(length=nodes)
for (n in 1:nodes)  rich[n] <- sum(venn_tree[n,]>0)
#for (n in 1:nodes)  rich[n] <- length(subset(venn_tree[n,],venn_tree[n,]>0))
pabm <- get_possible_branches_for_all_characters(simchmatrix,branchings,ab,venn_tree) # possible available branches matrix
pcc <- vector(length=nchars)  # possible character changes

for (ch in 1:nchars)  {
	x <- unique(permute(subset(pabm[ch,],pabm[ch,]>0)))
	pcc[ch] <- length(x)
	for (c in 1:pcc[ch])  		pabm[ch,c] <- x[c]
	if (pcc[ch]<length(ab))   for (c in (pcc[ch]+1):length(ab))   pabm[ch,c] <- 0
	}
#rm(ab,list=subset(bcc[1,],bcc[1,]>0))
# first, make sure that all states appear
for (ch in 1:nchars)	{
	if (states[ch]>=2)	{
    use_br <- sort(pabm[ch,1:(states[ch]-1)],decreasing=TRUE)   # do high nodes first, then low nodes, then otus
		for (c in 1:(states[ch]-1))	{
			br <- use_br[c]
			if (br<=notu)	{
				simchmatrix[br,ch] <- c
				} else {
				ndd <- br-notu
				for (r in 1:rich[ndd])	{
					s <- venn_tree[ndd,r]
					if (simchmatrix[s,ch]==0)	simchmatrix[s,ch] <- c
					}
				}
			char_changes[ch] <- 1+char_changes[ch]
			branch_changes[br] <- branch_changes[br]+1
			}	# add state c to matrix
		}	# only do variant characters
	}

simstates <- vector(length=nchars)
for (ch in 1:nchars)	simstates[ch] <- sum(unique(simchmatrix[,ch])>=0)

# second, make sure that all states appear
for (b in 1:fb)	{
	br <- ab[b]
	if (branch_changes[br]==0)	{
		ch <- ceiling(nchars*runif(1))
		while (states[ch]<2 || char_changes[ch]>=mx_ch[ch])	ch <- ceiling(nchars*runif(1))
		if (states[ch]==2)	{
			dstates <- c(1,0)
			} else if (types[ch]==0)	{
			dstates <- scramble_multistates(states[ch])	
			} 
		if (br<=notu)	{
			if (simchmatrix[br,ch]!=UNKNOWN && simchmatrix[br,ch]!=INAP)	{
				cc <- simchmatrix[br,ch]
				simchmatrix[br,ch] <- dstates[cc+1]
				} else {
				ndg <- br-notu
				for (d in 1:length(subset(venn_tree[ndg,],venn_tree[ndg,]>0)))	{
					sp <- venn_tree[ndg,d]
					if (simchmatrix[sp,ch]!=UNKNOWN && simchmatrix[sp,ch]!=INAP)	simchmatrix[sp,ch] <- dstates[1+simchmatrix[sp,ch]]
					}	
				}
			#change_character_on_branch(ch,branchings[b])
			}
		branch_changes[br] <- branch_changes[br]+1
		char_changes[ch] <- char_changes[ch]+1
#		bcc[ch,char_changes[ch]] <- br
		}	# end case where branch needed changes
	}

#third, tally compatibility at this point
delta <- sum(char_changes)
simcompmat <- compatibility_matrix(simchmatrix,states,types,UNKNOWN,INAP)
character_compats <- rowSums(simcompmat)-1
simcompat[delta] <- sum(character_compats)/2
#simcompat[delta] <- total_compatibility(simchmatrix,nchars,states,types,notu,UNKNOWN,INAP)
steps_v_compat <- c(delta,simcompat[delta])
if (repl>0)	print(c(repl,delta,simcompat[delta]))

d <- delta
counter <- 0
while (simcompat[d]>ttl_compat)	{
	counter <- counter+1
	d <- d+1
	ch <- ceiling(nchars*runif(1))
	while (ch>nchars || char_changes[ch]>=mx_ch[ch])	ch <- ceiling(nchars*runif(1))
	c <- char_changes[ch]+1
	br <- pabm[ch,c]
	
	prior_char <- simchmatrix[,ch]
	
	if (states[ch]==2)	{
		dstates <- c(1,0)       # this is a transition matrix: 0->1, 1->0
		} else if (types[ch]==0)	{
		dstates <- scramble_multistates(states[ch])	# this transition matrix will have either 0->1+1->2+2->0 or 0->2+1->0+2->1 for a 3 state character
		} 
#	prior_compat <- compatibility_of_a_character <- (ch1,simchmatrix,nchars,states,types,notu,UNKNOWN,INAP)
	prior_compat <- simcompmat[ch,]
	if (br<=notu)	{
		if (simchmatrix[br,ch]!=UNKNOWN && simchmatrix[br,ch]!=INAP)	{
			simchmatrix[br,ch] <- dstates[1+simchmatrix[br,ch]]
			} else {
			ndh <- br-notu
			for (d in 1:length(subset(venn_tree[ndh,],venn_tree[ndh,]>0)))	{
				sp <- venn_tree[ndh,d]
				if (simchmatrix[sp,ch]!=UNKNOWN && simchmatrix[sp,ch]!=INAP)	simchmatrix[sp,ch] <- dstates[1+simchmatrix[sp,ch]]
				}	
			}
		#change_character_on_branch(ch,branchings[b])
		}
	# new vector of compatibilities
	new_compat <- compatibility_of_a_character(ch,simchmatrix,states,types,UNKNOWN,INAP)
	# update compatibility matrix
	# tally new matrix compatibility
	#simcompat[d] <- (sum(simcompmat)-nchars)/2
	simcompat[d] <- simcompat[d-1]-(sum(prior_compat)-sum(new_compat))
	if (repl>0 && counter%%10==0)	print(c(d,simcompat[d]))
	steps_v_compat <- rbind(steps_v_compat,c(d,simcompat[d]))
	if (simcompat[d]>ttl_compat)	{
	# update character compatibilities that might have been altered
		branch_changes[br] <- branch_changes[br]+1
		char_changes[ch] <- char_changes[ch]+1
		simcompmat[ch,] <- new_compat
		simcompmat[,ch] <- new_compat
		for (c2 in 1:nchars)	{
			if (new_compat[c2]!=prior_compat[c2])	{
				if (new_compat[c2]==0)	{
					character_compats[c2] <- character_compats[c2]-1
					} else if (new_compat[c2]==1)	{
					character_compats[c2] <- character_compats[c2]+1
					}
				}	# case of mismatch
			}
		}	else if (simcompat[d]<ttl_compat) {
		simchmatrix[,ch] <- prior_char
		d <- d-1
		}
	}

#output <- list(d,simcompmat)
return (steps_v_compat)
}

# simulate character evolution replicating (or nearly so!) observed compatibility
### modified from evolve_compatibility_over_N_changes on 2017-08-29
## two steps: 1) race up to compatibility;
##			  2) if overshat, backup and try 100 times till it's found
emuli_observed_compatibility <- function(init_chmatrix,ttl_compat,venn_tree,branchings,states,types,UNKNOWN,INAP,repl=-1)	{
# initial_character matrix
# ttl_compat: compatibility to which to evolve
# venn_tree: a tree giving all of the observed descendants (direct & indirect)
# branchings: length of each branch
# nchars: number of characters
# states: number of states for each character
# types: types for each character (0: unordered; 1: ordered)
# maxsteps: maximum number of steps to do
# UNKNOWN: numeric code for "?"
# INAP: numeric code for "-"
# repl: replication number
nchars <- ncol(init_chmatrix)
notu <- nrow(init_chmatrix)	# number of observed taxa in tree
nodes <- nrow(venn_tree)
simchmatrix <- init_chmatrix
for (ch in 1:nchars)	simchmatrix[simchmatrix[,ch]>=0,ch] <- 0
unique_ab <- (1:length(branchings))[branchings>0]

branch_changes <- vector(length=notu+nodes)
char_changes <- vector(length=nchars)

ab <- c()
for (b in 1:length(unique_ab))	{
	br <- unique_ab[b]
	ab <- c(ab,rep(br,branchings[br]))
	}
tb <- length(ab)		# total branches; = branches above, so why both?
fb <- length(unique_ab)		# free branches, with each branch that can change listed only once
bcc <- matrix(0,nchars,fb)	# branch changes per character

rich <- vector(length=nodes)
for (n in 1:nodes)  rich[n] <- sum(venn_tree[n,]>0)
pabm <- get_possible_branches_for_all_characters(simchmatrix,branchings,ab,venn_tree) # possible available branches matrix
mx_ch <- vector(length=nchars)
for (ch in 1:nchars)	mx_ch[ch] <- sum(pabm[ch,]>0)
pcc <- vector(length=nchars)  # possible character changes

for (ch in 1:nchars)  {
	x <- unique(permute(subset(pabm[ch,],pabm[ch,]>0)))
	pcc[ch] <- length(x)
	for (c in 1:pcc[ch])  		pabm[ch,c] <- x[c]
	if (pcc[ch]<length(ab))   for (c in (pcc[ch]+1):length(ab))   pabm[ch,c] <- 0
	}

simcompatmat <- matrix(1,nchars,nchars)
prior_compat <- character_compats <- nchars-1

d <- 0
simcompat <- c()

# first, make sure that all states appear
for (ch in 1:nchars)	{
	if (states[ch]>=2)	{
		prior_compat <- simcompatmat[ch,]
    	use_br <- sort(pabm[ch,1:(states[ch]-1)],decreasing=TRUE)   # do high nodes first, then low nodes, then otus
		for (c in 1:(states[ch]-1))	{
			br <- use_br[c]
			if (br<=notu)	{
				simchmatrix[br,ch] <- c
				} else {
				ndd <- br-notu
				for (r in 1:rich[ndd])	{
					s <- venn_tree[ndd,r]
					if (simchmatrix[s,ch]==0)	simchmatrix[s,ch] <- c
					}
				}
			char_changes[ch] <- 1+char_changes[ch]
			branch_changes[br] <- branch_changes[br]+1
			d <- d+1
			new_compat <- compatibility_of_a_character(ch,simchmatrix,states,types,UNKNOWN,INAP)
			simcompatmat[ch,] <- simcompatmat[,ch] <- new_compat
			simcompat <- c(simcompat,sum(simcompatmat[lower.tri(simcompatmat)]))
			}	# add state c to matrix
		}	# only do variant characters
	}

simstates <- vector(length=nchars)
for (ch in 1:nchars)	simstates[ch] <- sum(unique(simchmatrix[,ch])>=0)

# second, make sure that all branches have change
counter <- 0
if (min(branch_changes[unique_ab])==0)	{
	needy_brs <- unique_ab[branch_changes[unique_ab]==0]
	for (b in 1:length(needy_brs))	{
		br <- needy_brs[b]	### make sure that is in all such routines!!!!
		counter <- counter+1
		d <- d+1
		ch <- ceiling(nchars*runif(1))
		while (states[ch]<2 || char_changes[ch]>=mx_ch[ch])	ch <- ceiling(nchars*runif(1))
		prior_compat <- simcompatmat[ch,]
		if (states[ch]==2)	{
			dstates <- c(1,0)
			} else if (types[ch]==0)	{
			dstates <- scramble_multistates(states[ch])	
			} 
		if (br<=notu)	{
			if (simchmatrix[br,ch]!=UNKNOWN && simchmatrix[br,ch]!=INAP)	{
				cc <- simchmatrix[br,ch]
				simchmatrix[br,ch] <- dstates[cc+1]
				} else {
				ndg <- br-notu
				for (d in 1:length(subset(venn_tree[ndg,],venn_tree[ndg,]>0)))	{
					sp <- venn_tree[ndg,d]
					if (simchmatrix[sp,ch]!=UNKNOWN && simchmatrix[sp,ch]!=INAP)	simchmatrix[sp,ch] <- dstates[1+simchmatrix[sp,ch]]
					}
				}
			}
		branch_changes[br] <- branch_changes[br]+1
		char_changes[ch] <- char_changes[ch]+1
		new_compat <- compatibility_of_a_character(ch,simchmatrix,states,types,UNKNOWN,INAP)
		simcompat[d] <- simcompat[d-1]-(sum(prior_compat)-sum(new_compat))
		}	# end case where branch needed changes
	}

#third, tally compatibility at this point
#steps_v_compat <- c(d,simcompat[d])
d <- length(simcompat)
if (repl>0)	print(c(repl,d,simcompat[d]))
counter <- 1
while (simcompat[d]>ttl_compat)	{
### while loop added 2017-06-03
#		simchmatrix_last <- simchmatrix
	br <- 0
	while (br==0)	{
		ch <- ceiling(nchars*runif(1))
		while (ch>nchars || char_changes[ch]>=mx_ch[ch])	ch <- ceiling(nchars*runif(1))
		c <- char_changes[ch]+1
		br <- pabm[ch,c]
		}
	prior_char <- simchmatrix[,ch]
	prior_compat <- simcompatmat[ch,]

	if (states[ch]==2)	{
		dstates <- c(1,0)       # this is a transition matrix: 0->1, 1->0
		} else if (types[ch]==0)	{
		dstates <- scramble_multistates(states[ch])	# this transition matrix will have either 0->1+1->2+2->0 or 0->2+1->0+2->1 for a 3 state character
		} 
#	prior_compat <- compatibility_of_a_character <- (ch1,simchmatrix,nchars,states,types,notu,UNKNOWN,INAP)
	if (br<=notu)	{
		if (simchmatrix[br,ch]!=UNKNOWN && simchmatrix[br,ch]!=INAP)	{
			simchmatrix[br,ch] <- dstates[1+simchmatrix[br,ch]]
			} else {
			ndh <- br-notu
			for (d in 1:length(subset(venn_tree[ndh,],venn_tree[ndh,]>0)))	{
				sp <- venn_tree[ndh,d]
				if (simchmatrix[sp,ch]!=UNKNOWN && simchmatrix[sp,ch]!=INAP)	simchmatrix[sp,ch] <- dstates[1+simchmatrix[sp,ch]]
				}	
			}
		#change_character_on_branch(ch,branchings[b])
		}
	# new vector of compatibilities
#		ccc <- compatibility_matrix(simchmatrix,states,types,UNKNOWN,INAP)
	new_compat <- compatibility_of_a_character(ch,simchmatrix,states,types,UNKNOWN,INAP)
	# update compatibility matrix
	# tally new matrix compatibility
	#simcompat[d] <- (sum(simcompatmat)-nchars)/2
	d <- d+1
#		simcompat <- c(simcompat,simcompat[d-1]-(sum(prior_compat)-sum(new_compat)))
	simcompatmat[ch,] <- new_compat
	simcompatmat[,ch] <- new_compat
	simcompat <- c(simcompat,sum(simcompatmat[lower.tri(simcompatmat)]))
	if (repl>0 && counter%%10==repl)	print(c(d,simcompat[d]))
	counter <- counter+1
#	steps_v_compat <- rbind(steps_v_compat,c(d,simcompat[d]))
	# update character compatibilities that might have been altered
	branch_changes[br] <- branch_changes[br]+1
	char_changes[ch] <- char_changes[ch]+1
	for (c2 in 1:nchars)	{
		if (new_compat[c2]!=prior_compat[c2])	{
			if (new_compat[c2]==0)	{
				character_compats[c2] <- character_compats[c2]-1
				} else if (new_compat[c2]==1)	{
				character_compats[c2] <- character_compats[c2]+1
				}
			}	# case of mismatch
		}
	}

if (simcompat[d]<ttl_compat)	{
	simchmatrix_closest <- simchmatrix
	simchmatrix_closest[,ch] <- prior_char
	closest <- simcompat[d-1]-ttl_compat
	max_br <- ncol(pabm)
	}	# set up closest compatibility with too much: we'll use that in a pinch

counter <- 0
while (simcompat[d]!=ttl_compat && counter<100)	{
	# Undo prior run (if we improved, then we are keeping it the same here)
	simcompatmat[,ch] <- simcompatmat[ch,] <- prior_compat
	simchmatrix[,ch] <- prior_char
	branch_changes[br] <- branch_changes[br]-1
	char_changes[ch] <- char_changes[ch]-1
	br <- 0
	while (br==0)	{
		ch <- ceiling(nchars*runif(1))
		while (ch>nchars || char_changes[ch]>=mx_ch[ch])	ch <- ceiling(nchars*runif(1))
		c <- char_changes[ch]+1
		brbr <- pabm[ch,c:mx_ch[ch]]
		rem_poss_br <- pabm[ch,c:mx_ch[ch]][pabm[ch,c:mx_ch[ch]]>0]
		if (length(rem_poss_br)>0)	{
			cc <- 0
			while (cc==0)	cc <- ceiling(runif(1)*length(rem_poss_br))
			br <- rem_poss_br[cc]
			}
		}
	prior_char <- simchmatrix[,ch]
	prior_compat <- simcompatmat[ch,]

	if (states[ch]==2)	{
		dstates <- c(1,0)       # this is a transition matrix: 0->1, 1->0
		} else if (types[ch]==0)	{
		dstates <- scramble_multistates(states[ch])	# this transition matrix will have either 0->1+1->2+2->0 or 0->2+1->0+2->1 for a 3 state character
		} 
#	prior_compat <- compatibility_of_a_character <- (ch1,simchmatrix,nchars,states,types,notu,UNKNOWN,INAP)
	if (br<=notu  && br>0)	{
		if (simchmatrix[br,ch]!=UNKNOWN && simchmatrix[br,ch]!=INAP)	{
			simchmatrix[br,ch] <- dstates[1+simchmatrix[br,ch]]
			} else if (br>notu) {
			ndh <- br-notu
			for (d in 1:length(subset(venn_tree[ndh,],venn_tree[ndh,]>0)))	{
				sp <- venn_tree[ndh,d]
				if (simchmatrix[sp,ch]!=UNKNOWN && simchmatrix[sp,ch]!=INAP)	simchmatrix[sp,ch] <- dstates[1+simchmatrix[sp,ch]]
				}	
			}
		#change_character_on_branch(ch,branchings[b])
		}
	# new vector of compatibilities
#		ccc <- compatibility_matrix(simchmatrix,states,types,UNKNOWN,INAP)
	new_compat <- compatibility_of_a_character(ch,simchmatrix,states,types,UNKNOWN,INAP)
	# update compatibility matrix
	# tally new matrix compatibility
	#simcompat[d] <- (sum(simcompatmat)-nchars)/2
#		d <- d+1
#		simcompat <- c(simcompat,simcompat[d-1]-(sum(prior_compat)-sum(new_compat)))
	simcompatmat[ch,] <- new_compat
	simcompatmat[,ch] <- new_compat
	simcompat[d] <- sum(simcompatmat[lower.tri(simcompatmat)])
	if (repl>0 && counter%%10==repl)	print(c(d,simcompat[d]))
	counter <- counter+1
#	steps_v_compat[dd,2] <- simcompat[d]
	# update character compatibilities that might have been altered
	branch_changes[br] <- branch_changes[br]+1
	char_changes[ch] <- char_changes[ch]+1
	for (c2 in 1:nchars)	{
		if (new_compat[c2]!=prior_compat[c2])	{
			if (new_compat[c2]==0)	{
				character_compats[c2] <- character_compats[c2]-1
				} else if (new_compat[c2]==1)	{
				character_compats[c2] <- character_compats[c2]+1
				}
			}	# case of mismatch
		}
#	print(c(simcompat[d],ttl_compat))
	
	# if improvement without getting it exactly right
	if (closest > (simcompat[d]-ttl_compat) && (simcompat[d]-ttl_compat)>0)	{
		closest <- simcompat[d]-ttl_compat
		simchmatrix_closest <- simchmatrix
		prior_char <- simchmatrix[,ch]
		prior_compat <- simcompatmat[ch,]
		d <- d+1
		simcompat <- c(simcompat,ttl_compat-1)
#		dd <- dd+1
		branch_changes[br] <- branch_changes[br]+1	# these will both be decremented above
		char_changes[ch] <- char_changes[ch]+1		# these will both be decremented above
		}
	}

output <- list(simchmatrix,simcompatmat,d)
names(output) <- c("Simulated_Char_Matrix","Simulated_Compatibility_Matrix","Steps")
return (output)
}

# routine to list all branches on which character change can be reconstructed for each character
get_possible_branches_for_all_characters <- function(simchmatrix,branchings,ab,venn_tree)	{
# find branches on which each character can change.  (Accommodates unknowns & inapplicables)
# simchmatrix: character matrix
# branchings: 
# rich: number of nodes
# ab: number of available branches (e.g., some scored states above it)
nchars <- ncol(simchmatrix)
notu <- nrow(simchmatrix)
poss_brs <- length(ab)
pabm <- matrix(0,nchars,poss_brs)
for (ch in 1:nchars)	{
	a <- 1
	for (b in 1:poss_brs)	{
		if (ab[b]<=notu)	{
			s <- ab[b]
			if (simchmatrix[s,ch]!=UNKNOWN && simchmatrix[s,ch]!=INAP)	{
				pabm[ch,a]<-s
				a <- a+1
				}	# case where species is scored
			} else {
			nda <- ab[b]-notu  ### PROBLEM APPEARS HERE!!!!
			accept <- 0
##			  for (r in 1:rich[nda])	{
			for (r in 1:sum(venn_tree[nda,]>0))	{
				s <- venn_tree[nda,r]
				if (simchmatrix[s,ch]!=UNKNOWN && simchmatrix[s,ch]!=INAP)	{
					accept <- 1
				  	r <- sum(venn_tree[nda,])
					}	# case where a clade member is scored
				}
			if (accept==1)	{
				pabm[ch,a]<-nda+notu
				a <- a+1
				}
			}	# end case of node
		}
	}
char_names <- vector(length=nchars)
for (c in 1:nchars)	{
	if (nchars<100)	{
		if (c<10)	{
			char_names[c] <- paste("ch_0",c,sep="")
			}	else	{
			char_names[c] <- paste("ch_",c,sep="")	
			}
		}	else	{
		if (c<10)	{
			char_names[c] <- paste("ch_00",c,sep="")
			} else if (c<100)	{
			char_names[c] <- paste("ch_0",c,sep="")
			}	else	{
			char_names[c] <- paste("ch_",c,sep="")
			}
		}
	}

brnch_names <- vector(length=poss_brs)
for (b in 1:poss_brs)	{
	if (poss_brs<100)	{
		if (b<10)	{
			brnch_names[b] <- paste("br_0",b,sep="")
			}	else	{
			brnch_names[b] <- paste("br_",b,sep="")	
			}
		}	else	{
		if (b<10)	{
			brnch_names[b] <- paste("br_00",b,sep="")
			} else if (c<100)	{
			brnch_names[b] <- paste("br_0",b,sep="")
			}	else	{
			brnch_names[b] <- paste("br_",b,sep="")
			}
		}
	}
rownames(pabm) <- char_names
colnames(pabm) <- brnch_names
return(pabm)
}

# routine to create unique a <-> b <-> c transitions
scramble_multistates <- function(nstates)	{
dstates <- vector(length=nstates)
for (i in 1:nstates)	dstates[i]<-i-1
for (i in 1:nstates)	{
	p <- i+ceiling((nstates-i)*runif(1))
	q <- dstates[i]
	dstates[i]<-dstates[p]
	dstates[p]<-q
	}
return(dstates)
}

# routine to create P[compatibility | steps]
# Ugh: I messed with this by mistake.  Restore it!
accersi_expected_compatibility_given_steps <- function(init_chmatrix,runs=100,lambda=0.6,mu=0.5,freqrat=0.5,bifurcation=1,contemporaneous=FALSE,hidden_reversals=TRUE,maxsteps=-1,temp_prec=0.1,UNKNOWN,INAP)	{
notu <- nrow(init_chmatrix)
nchars <- ncol(init_chmatrix)
states <- c()
for (c in 1:nchars)	{
	dummy <- init_chmatrix[(1:notu)[init_chmatrix[,c]>=0],c]
	states <- c(states,1+max(dummy)-min(dummy))
	}
types <- rep(0,nchars)
if (maxsteps==-1)	{
	taxa_scored <- count_scored_otu_per_character(init_chmatrix)
	maxsteps_per_char <- accersi_maximum_parsimony_steps_per_character(init_chmatrix,states)
	maxsteps <- ceiling((sum(maxsteps_per_char)+sum(taxa_scored))/2)
	}
sim_results <- c()
for (r in 1:runs)	{
	apprise <- paste("doing run",r,sep=" ")
	print(apprise)
	if (contemporaneous)	{
		simulation <- evolve_to_standing_richness_S(S=notu,lambda,mu,bifurcation,temp_prec)
		} else	{
		simulation <- evolve_to_sampled_richness_S(S=notu,lambda,mu,freqrat,bifurcation=0,temp_prec=0.1)
		}
#	pt_done <- 1
	venn_tree <- simulation$Venn_Tree
	nNodes <- dim(venn_tree)[1]
	branchings <- simulation$Branchings[1:(notu+nNodes)]
#	simulation$Branchings[(notu+1):(notu+nNodes)]
	stc <- evolve_compatibility_over_N_changes(N=maxsteps,init_chmatrix,venn_tree,branchings,nchars,states,types,hidden_reversals,UNKNOWN,INAP,repl=0)
#	pt_done <- 2
#	new_results <- c(rep(0,(min(stc[,1])-1)),stc[,2])
	sim_results <- rbind(sim_results,new_results)
	}
colnames(sim_results) <- 1:maxsteps
#sim_results <- sim_results/runs
return(sim_results)
}

# routine to create P[compatibility | steps] by finding X-1 & X steps crossing observed compatibility
accersi_ave_steps_to_reach_compatibility <- function(init_chmatrix,runs=100,lambda=0.6,mu=0.5,freqrat=0.5,bifurcation=1,contemporaneous=FALSE,hidden_reversals=TRUE,maxsteps=-1,temp_prec=0.1,UNKNOWN,INAP)	{
notu <- nrow(init_chmatrix)
nchars <- ncol(init_chmatrix)
states <- c()
for (c in 1:nchars)	{
	dummy <- init_chmatrix[(1:notu)[init_chmatrix[,c]>=0],c]
	states <- c(states,1+max(dummy)-min(dummy))
	}
types <- rep(0,nchars)
obs_comp_matrix <- compatibility_matrix(init_chmatrix,states,types,UNKNOWN,INAP)
obs_compat <- sum(obs_comp_matrix[lower.tri(obs_comp_matrix)])
if (maxsteps==-1)	{
	taxa_scored <- count_scored_otu_per_character(init_chmatrix)
	maxsteps_per_char <- accersi_maximum_parsimony_steps_per_character(init_chmatrix,states)
	maxsteps <- ceiling((sum(maxsteps_per_char)+sum(taxa_scored))/2)
	}
sim_results <- c()
for (r in 1:runs)	{
	apprise <- paste("doing run",r,sep=" ")
	print(apprise)
	if (contemporaneous)	{
		simulation <- evolve_to_standing_richness_S(S=notu,lambda,mu,bifurcation,temp_prec)
		} else	{
		simulation <- evolve_to_sampled_richness_S(S=notu,lambda,mu,freqrat,bifurcation=0,temp_prec=0.1)
		}
#	pt_done <- 1
	venn_tree <- simulation$Venn_Tree
	nNodes <- dim(venn_tree)[1]
	branchings <- simulation$Branchings[1:(notu+nNodes)]
	stc <- evolve_to_particular_compatibilty(init_chmatrix,ttl_compat=obs_compat,venn_tree,branchings,states,types,UNKNOWN,INAP,repl=0)
	sim_results <- c(sim_results,stc)
#	simulation$Branchings[(notu+1):(notu+nNodes)]
#	stc <- evolve_compatibility_over_N_changes(N=maxsteps,init_chmatrix,venn_tree,branchings,nchars,states,types,hidden_reversals,UNKNOWN,INAP,repl=0)
#	pt_done <- 2
#	new_results <- c(rep(0,(min(stc[,1])-1)),stc[,2])
#	stc <- evolve_to_particular_compatibilty(ttl_compat,venn_tree,branchings,nchars,states,types,maxsteps,UNKNOWN,INAP,repl=1)
#	sim_results <- rbind(sim_results,new_results)
	}
#colnames(sim_results) <- 1:maxsteps
#sim_results <- sim_results/runs
return(sim_results)
}

# routine to create P[compatibility | steps]
accersi_most_prob_compatibility_given_steps <- function(init_chmatrix,runs=100,lambda=0.6,mu=0.5,freqrat=0.5,bifurcation=1,contemporaneous=FALSE,hidden_reversals=TRUE,maxsteps=-1,temp_prec=0.1,UNKNOWN,INAP,repl=-1)	{
notu <- nrow(init_chmatrix)
nchars <- ncol(init_chmatrix)
states <- c()
for (c in 1:nchars)	{
	dummy <- init_chmatrix[(1:notu)[init_chmatrix[,c]>=0],c]
	states <- c(states,1+max(dummy)-min(dummy))
	}
types <- rep(0,nchars)
if (maxsteps==-1)	{
	taxa_scored <- count_scored_otu_per_character(init_chmatrix)
	maxsteps_per_char <- accersi_maximum_parsimony_steps_per_character(init_chmatrix,states)
	maxsteps <- ceiling((sum(maxsteps_per_char)+sum(taxa_scored))/2)
	}
sim_results <- c()
for (r in 1:runs)	{
	apprise <- paste("doing run",r,sep=" ")
	print(apprise)
	if (contemporaneous)	{
		simulation <- evolve_to_standing_richness_S(S=notu,lambda,mu,bifurcation,temp_prec)
		} else	{
		simulation <- evolve_to_sampled_richness_S(S=notu,lambda,mu,freqrat,bifurcation=0,temp_prec=0.1)
		}
#	pt_done <- 1
	venn_tree <- simulation$Venn_Tree
	nNodes <- dim(venn_tree)[1]
	branchings <- simulation$Branchings[1:(notu+nNodes)]
#	simulation$Branchings[(notu+1):(notu+nNodes)]
	stc <- evolve_up_to_compatibility(ttl_compat,init_chmatrix,venn_tree,branchings,nchars,states,types,UNKNOWN,INAP,repl)
#	pt_done <- 2
	new_results <- c(rep(0,(min(stc[,1])-1)),stc[,2])
#	stc <- evolve_to_particular_compatibilty(ttl_compat,venn_tree,branchings,nchars,states,types,maxsteps,UNKNOWN,INAP,repl=1)
	sim_results <- rbind(sim_results,new_results)
	}
colnames(sim_results) <- 1:maxsteps
sim_results <- sim_results/runs
return(sim_results)
}

# routine to find standing richness
tally_richness_from_continuous_ranges <- function(ranges,temp_prec=0.1)	{
stg_ranges <- ceiling(ranges/temp_prec)
stg_richness <- vector(length=max(stg_ranges))
etus <- dim(stg_ranges)[1]
#test <- c()
for (i in 1:etus)	{
	stg_richness[(stg_ranges[i,1]:stg_ranges[i,2])] <- stg_richness[(stg_ranges[i,1]:stg_ranges[i,2])]+1
#	if (stg_ranges[i,1]<=present && stg_ranges[i,2]>=present)
#		test <- c(test,i)
	}
return(stg_richness)
}
