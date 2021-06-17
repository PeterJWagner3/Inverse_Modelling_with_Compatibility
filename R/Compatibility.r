# written by Peter Wagner & Peter D Smits
#' Find the swing between three adjacent state pairs
#'
#' When determining compatibility between unordered multistate characters it is 
#' necessary to determine the position of the swing pair when there are three 
#' observed. This function is used internally.
#'
#' @param pa 3x2 matrix of state pairs
#' @param u1 state identities (values) of the first column
#' @param u2 state identities (values) of the second column

count_scored_otu_per_character <- function(chmatrix,UNKNOWN=-11,INAP=-22)	{
nch <- dim(chmatrix)[2]
notu <- dim(chmatrix)[1]
scored <- vector(length=nch)
for (c in 1:nch)
	scored[c] <- notu - (sum(chmatrix[,c]==UNKNOWN)+sum(chmatrix[,c]==INAP))
return(scored)
}

reduce_to_unique_pairs <- function (ch1, ch2, st1, st2, UNKNOWN, INAP)	{
charpair <- cbind(ch1,ch2)
# Remove pairs with unknowns or inapplicables: they do not count
if (any(charpair == UNKNOWN) || any(charpair == INAP) || any(charpair < 0)) {
	rmz <- which(charpair == UNKNOWN | charpair == INAP | charpair < 0, arr.ind = T)[, 1]
	charpair <- charpair[-rmz, ]
	}
# reduce pairs to unique ones.  Then sort on first then second column  ADDED 2017-05-31
charpair <- unique(charpair)
pair_order <- c()
for (i in 1:nrow(charpair))	{
	pair_order <- c(pair_order,1+(charpair[i,1]*st2)+charpair[i,2])
	}
charpair <- charpair[order(pair_order),]
# get all of the pairs that we have
return (charpair)
}

#shifts <- c(1,2,3,4)
#pshifts_norm <- dnorm(shifts) 
#pshifts_exp <- pshifts_norm[1]*exp((min(log(pshifts_norm))-max(log(pshifts_norm)))/(max(shifts)-min(shifts)))^(shifts-1)
#plot(shifts,log(pshifts_norm),pch=21,bg="red",cex=2)
#points(shifts,log(pshifts_exp),pch=21,bg="green",cex=2)
# for unordered multistates.  This makes sure that pairs with only 3 combos have at least one
#	set that is not intermediate between two others.
#
find_swing_pairs <- function (pa, u1, u2) {
  
if (any(pa[, 1] == u1[1] & pa[, 2] == u2[1]) & any(pa[, 1] == u1[2] & pa[, 2] == u2[2])) {
	# heterogeneous swing
	# one of the maxes and one of the mins
	# can't be both
	p1 <- ifelse(sum(pa[, 1] == u1[1]) > sum(pa[, 1] == u1[2]), 1, 2)
	p2 <- ifelse(sum(pa[, 2] == u2[1]) > sum(pa[, 2] == u2[2]), 1, 2)
    swing <- pa[pa[, 1] == u1[p1] & pa[, 2] == u2[p2], ]
    }	else {
	# homogeneous swing
	# either both of the maxes or both of the mins
	# can't be both
	o <- list()
	o[[1]] <- pa[, 1] == max(u1) & pa[, 2] == max(u2)
	o[[2]] <- pa[, 1] == min(u1) & pa[, 2] == min(u2)
	y <- ifelse(any(o[[1]]), 1, 2)
	swing <- pa[o[[y]], ]
	}
return(swing)
}

#' pair compatibility
#'
#' Calculate the compatibility between any two unordered characters. 
#' Compatibility is a tree-free method of determining the plausibility 
#' of of two characters having evolved with homoplasy. This function does this 
#' as a hybrid of various methods devised in the 80s to account for any 
#' number of character states for unordered characters.
#'
#' @param ch1 either a vector of scored states for a character or a two column 
#' matrix representing two characters
#' @param ch2 a vector of scored states for a character. optional if ch1 is a 
#' two column matrix.
#' @param st1 number of states in the first character (optional)
#' @param st2 number of states in the second character (optional)
#' @param t1 type of the first character (0 = unordered/binary, 1 = ordered) 
#' DOES NOTHING. DO NOT CHANGE
#' @param t2 type of the second character
#' @param UNKNOWN value of unknown character state scores
#' @param INAP value of inapplicable character state scores
#' @keywords
#' @export
#' @examples
pair_compatibility <- function(ch1,ch2,st1,st2,t1,t2,UNKNOWN=-11,INAP=-22) {
  #  Returns:
  #    value of 0 or 1
  #    0 is incompatible
  #    1 is compatible
  
  # make sure the input is in the right format
  # weirdness?
comp <- 1;
if (!missing(ch2) && identical(length(dim(ch1)), 1)) {
	stop("ch2 exists but ch1 is not a 1-d vector! one or the other please");
	}
  # ch1 is a pair?
if (missing(ch2) && identical(dim(ch1)[2], 2)) {
	charpair <- ch1;
	}	else if (missing(ch2) && !identical(dim(ch1)[2], 2))	{
	stop("char pair is not of the right shape. 2 column matrix please.");
	}
# ch1 and ch2 both exist? are they the same length?
if (!missing(ch2) && identical(length(ch1), length(ch2))) {
	charpair <- cbind(ch1, ch2);
	} else if (!missing(ch2) && !identical(length(ch1), length(ch2))) {
	stop("ch1 and ch2 are not of equal length");
	}
  
# handle states
if (missing(st1)) st1 <- max(charpair[, 1]) + 1;
if (missing(st2)) st2 <- max(charpair[, 2]) + 1;
  
# Remove pairs with unknowns or inapplicables: they do not count
if (any(charpair == UNKNOWN) || any(charpair == INAP)) {
	rmz <- which(charpair == UNKNOWN | charpair == INAP, arr.ind = T)[, 1];
	charpair <- charpair[-rmz, ];
	}
  
# get all of the pairs that we have
all_pairs <- unique(charpair);

# this is relevant only if there are 3+ pairs!
if (is.matrix(all_pairs))	{
	# binary is really easy.
	if (st1 <= 2 && st2 <= 2) {
		comp <- ifelse(nrow(all_pairs) > 3, 0, 1)
		}	else if (st1 == 1 || st2 == 1)	{
    	comp <- 1
  		} else {  # multistate
  	# all state pair combinations
  		sp1 <- combn(st1, 2) - 1;
  		sp2 <- combn(st2, 2) - 1;
    
    # if ordered multistate
    	if (t1 == 1)
    		sp1 <- sp1[, sp1[1, ] == sp1[2, ] - 1];
    	if (t2 == 1)
    		sp2 <- sp2[, sp2[1, ] == sp2[2, ] - 1];

       	pai <- c() # one of the few times you need to intialize
       	swing <- matrix(, ncol = 2);
       	for (ii in seq(ncol(sp1))) {
       		for (jj in seq(ncol(sp2))) {
       			c1 <- which(all_pairs[, 1] == sp1[1, ii] | 
                    all_pairs[, 1] == sp1[2, ii], 
                    arr.ind = T);
    			c2 <- which(all_pairs[, 2] == sp2[1, jj] | 
                    all_pairs[, 2] == sp2[2, jj], 
                    arr.ind = T);
    			pa <- all_pairs[intersect(c1, c2), ];
    			p <- nrow(as.matrix(pa))
    			n <- ifelse(p > 3, 0, 1)
    			pai <- c(pai, n)
    			if (any(c(t1, t2) == 0) & p == 3) {
    				u1 <- sp1[, ii];
    				u2 <- sp2[, jj];
    				swing <- na.omit(rbind(swing, find_swing_pairs(pa, u1, u2)));
    				}
    			}
    		}
    	circuit <- FALSE;
    	if (any(c(t1, t2) == 0) & nrow(swing) > 3) {
    		swing <- unique(swing)
    		state_diagram <- matrix(0, nrow = st1, ncol = st2)
    		state_diagram[swing + 1] <- 1
    		ci <- state_diagram
    		if (any(rowSums(state_diagram) < 2) | any(colSums(state_diagram) < 2)) {
    			ci <- state_diagram[rowSums(state_diagram) >= 2, 
                            colSums(state_diagram) >= 2]
      			}
    		circuit <- (sum(ci) == sum(state_diagram)) | 
    			(sum(ci) %% 2 == 0 && sum(ci) >= 4)
    		}
    	comp <- ifelse(any(pai == 0), 0, 1);
    	comp <- ifelse(circuit, 0, comp);
  		}
	}
return(comp)
}

pair_compatibility_fuzzy <- function(ch1,ch2,st1,st2,t1,t2,UNKNOWN=-11,INAP=-22) {
  #  Returns:
  #    value of 0 or 1
  #    0 is incompatible
  #    1 is compatible
  
  # make sure the input is in the right format
  # weirdness?
if (st2 > st1)	{
	chd <- ch1;
	td <- t1;
	std <- st1;
	ch1 <- ch2;
	t1 <- t2;
	st1 <- st2;
	ch2 <- chd;
	t2 <- td;
	st2 <- std;
	}
comp <- 1;
if (!missing(ch2) && identical(length(dim(ch1)), 1)) {
	stop("ch2 exists but ch1 is not a 1-d vector! one or the other please");
	}
  # ch1 is a pair?
if (missing(ch2) && identical(dim(ch1)[2], 2)) {
	charpair <- ch1;
	}	else if (missing(ch2) && !identical(dim(ch1)[2], 2))	{
	stop("char pair is not of the right shape. 2 column matrix please.");
	}

# ch1 and ch2 both exist? are they the same length?
if (!missing(ch2) && identical(length(ch1), length(ch2))) {
	charpair <- cbind(ch1, ch2);
	} else if (!missing(ch2) && !identical(length(ch1), length(ch2))) {
	stop("ch1 and ch2 are not of equal length");
	}
  
# handle states
if (missing(st1)) st1 <- max(charpair[, 1]) + 1;
if (missing(st2)) st2 <- max(charpair[, 2]) + 1;
  
# Remove pairs with unknowns or inapplicables: they do not count
if (any(charpair == UNKNOWN) || any(charpair == INAP)) {
	rmz <- which(charpair == UNKNOWN | charpair == INAP, arr.ind = T)[, 1];
	charpair <- charpair[-rmz, ];
	}
  
# get all of the pairs that we have
all_pairs <- unique(charpair);
all_pairs[,1] <- match(all_pairs[,1],sort(unique(all_pairs[,1])))-1;
all_pairs[,2] <- match(all_pairs[,2],sort(unique(all_pairs[,2])))-1;

pair_tests_1 <- c();
for (cs1 in 1:(st1-1))	for (cs2 in (cs1+1):st1)	pair_tests_1 <- cbind(pair_tests_1,c(cs1-1,cs2-1));
pair_tests_2 <- c();
for (cs1 in 1:(st2-1))	for (cs2 in (cs1+1):st2)	pair_tests_2 <- cbind(pair_tests_2,c(cs1-1,cs2-1));

compatible_subset <- comparable_subset <- c();
for (cp1 in 1:ncol(pair_tests_1))	{
	relv_pairs_1 <- (1:nrow(all_pairs))[all_pairs[,1] %in% pair_tests_1[,cp1]]
	for (cp2 in 1:ncol(pair_tests_2))	{
		relv_pairs_2 <- (1:nrow(all_pairs))[all_pairs[,2] %in% pair_tests_2[,cp2]]
		relv_pairs <- relv_pairs_1[relv_pairs_1 %in% relv_pairs_2];
		comparable_subset <- c(comparable_subset,ifelse (length(relv_pairs)>1,1,0));
		compatible_subset <- c(compatible_subset,ifelse ((length(relv_pairs)>1 && length(relv_pairs)<4),1,0));
		}
	}

return(c(sum(compatible_subset),sum(comparable_subset)));
}

compatibility_matrix <- function(chmatrix,states,types,UNKNOWN=-11,INAP=-22,fuzzy=F)	{
#	chmatrix: matrix of characters
#	types: 0 for unordered, 1 for ordered
#	fuzzy: if "true", then multistate characters will get compatibilities calculated for multiple breakdowns
#	e.g.: 00 10 11  00 02 12
nchars <- ncol(chmatrix);
comp_matrix <- matrix(0,nchars,nchars);
for (c1 in 1:(nchars-1))  {
	comp_matrix[c1,c1] <- 1;
	ch1 <- chmatrix[,c1];
	st1 <- states[c1];
	t1 <- types[c1];
	for (c2 in (c1+1):nchars)	{
		ch2 <- chmatrix[,c2];
		st2 <- states[c2];
		t2 <- types[c2];
		if (fuzzy)	{
			pair_comp <- pair_compatibility_fuzzy(ch1,ch2,st1,st2,t1,t2,UNKNOWN,INAP);
			if (pair_comp[2]>0)	{
				comp_matrix[c2,c1] <- comp_matrix[c1,c2] <- pair_comp[1]/pair_comp[2];
				} else	{
				comp_matrix[c2,c1] <- comp_matrix[c1,c2] <- 0.5;
				}
			} else	{
			comp_matrix[c2,c1] <- comp_matrix[c1,c2] <- pair_compatibility(ch1,ch2,st1,st2,t1,t2,UNKNOWN,INAP);
			}
		}
#	c1 <- c1+1
	}
comp_matrix[nchars,nchars] <- 1;
return (comp_matrix)
}

compatibility_per_character <- function(chmatrix,nchars,states,types,UNKNOWN=-11,INAP=-22)	{
character_compatibility <- vector(length=nchars)
for (c1 in 1:(nchars-1))  {
	for (c2 in (c1+1):nchars)  {
		c <- pair_compatibility(ch1=chmatrix[,c1],ch2=chmatrix[,c2],st1=states[c1],st2=states[c2],t1=types[c1], t2=types[c2],UNKNOWN,INAP)
		character_compatibility[c1] <- character_compatibility[c1]+c
		character_compatibility[c2] <- character_compatibility[c2]+c
		}
	}
return(character_compatibility)
}

# get the compatibility of a single character with all others
compatibility_of_a_character <- function(c1,chmatrix,states,types,UNKNOWN=-11,INAP=-22)	{
nch <- length(states)
character_compatibility <- vector(length=nch)
for (c2 in 1:nch)  {
	if (c2==c1)  {
		character_compatibility[c1] <- 1
		} else {
		character_compatibility[c2] <- pair_compatibility(chmatrix[,c1],chmatrix[,c2],st1=states[c1],st2=states[c2],t1=types[c1],t2=types[c2],UNKNOWN,INAP)
#		 <- character_compatibility[c2]+c
		}
	}
return(character_compatibility)
}

total_compatibility <- function(chmatrix,states,types,UNKNOWN=-11,INAP=-22)	{
nchars <- length(states)
ttl_compat <- 0
for (c1 in 1:(nchars-1))  {
	for (c2 in (c1+1):nchars)  {
		c <- pair_compatibility(ch1=chmatrix[,c1],ch2=chmatrix[,c2],st1=states[c1],st2=states[c2],t1=types[c1], t2=types[c2],UNKNOWN,INAP)
		ttl_compat <- ttl_compat+c
		}
	}
return(ttl_compat)
}

#ch1 <- c(0,0,0,0,0,1,1,1,0,2)
#ch2 <- c(0,0,0,1,1,1,1,1,2,1)
#fad <- c(1,1,2,2,2,3,3,4,4,5)
#lad <- c(1,2,3,2,4,4,3,4,5,5)
#ranges <- cbind(fad,lad)
#st1 <- 3
#st2 <- 3
#notu <- 10

date_character_pairs_raw <- function(notu, ranges, ch1, ch2, st1, st2, UNKNOWN=-11,INAP=-22)	{
pairs <- matrix(0,st1,st2)
pair_ranges <- matrix(-1,st1*st2,2)
for (s in 1:notu)	{
	if ((ch1[s]>=0 && ch2[s]>=0) && ((ch1[s]!=UNKNOWN && ch1[s]!=INAP) && (ch2[s]!=UNKNOWN && ch2[s]!=INAP)))	{
		pr <- 1+(ch1[s]*st2)+ch2[s]		# pr=1 for 00, pr=2 for 01; if st2=2, then pr=3 for 10 and 4 for 11
		if (pair_ranges[pr,1]==-1)	pair_ranges[pr,1] <- max(ranges[,2])
		pairs[ch1[s]+1,ch2[s]+1] <- pr
		if (ranges[s,1]<pair_ranges[pr,1])	pair_ranges[pr,1] <- ranges[s,1]
		if (ranges[s,2]>pair_ranges[pr,2])	pair_ranges[pr,2] <- ranges[s,2]
		#		if (ranges[1])
		}
	}
return (pair_ranges)
}

date_character_pairs <- function(ranges, pairs, ch1, ch2, st1, st2, UNKNOWN=-11,INAP=-22)	{
## added 2017-05-31
notu <- length(ch1)
charpair <- cbind(ch1,ch2)
pair_order <- c()
# number the pairs in their Xidecimal. 2 binary has 1 for 00, 2 for 01, 3 for 10, 4 for 11
for (i in 1:nrow(charpair))	{
	if (ch1[i]>=0 && ch2[i]>=0)	{
		pair_order <- c(pair_order,1+(charpair[i,1]*st2)+charpair[i,2])
		}	else	pair_order <- c(pair_order,0)
	}
# list numbers tied to unique pairs
unique_pairs <- sort(unique(pair_order[pair_order>0]))

# now, get the ranges of each pair by id'ing notus with each pair & getting their ranges
up <- 1
pair_ranges <- c()
while (up <= length(unique_pairs))	{
	stotu <- (1:notu)[pair_order==unique_pairs[up]]
	pair_ranges <- rbind(pair_ranges,c(min(ranges[stotu,1]),max(ranges[stotu,2])))
	up <- up+1
#	ranges[stotu,]
	}
return (pair_ranges)
}

reduce_to_unique_pairs_old <- function (ch1, ch2, st1, st2, notu, UNKNOWN=-11,INAP=-22)	{
found_pairs <- vector(length=(st1*st2))
#pairs <- matrix(-1,st1*st2,2)
fnd_prs <- 0
for (s in 1:notu)	{
	if ((ch1[s]>=0 && ch2[s]>=0) && ((ch1[s]!=UNKNOWN && ch1[s]!=INAP) && (ch2[s]!=UNKNOWN && ch2[s]!=INAP)))	{
		pr <- 1+(ch1[s]*st2)+ch2[s]		# pr=1 for 00, pr=2 for 01; if st2=2, then pr=3 for 10 and 4 for 11
		if (found_pairs[pr]==0)	{
			fnd_prs <- fnd_prs+1
			if (fnd_prs==1) {	
				pairs <- c(ch1[s],ch2[s])
			} else {
				pairs <- rbind(pairs,c(ch1[s],ch2[s]))
			}
			#			pairs[fnd_prs,1] <- ch1[s]
			#			pairs[fnd_prs,2] <- ch2[s]
			found_pairs[pr] <- 1
			}
		}
	}
return (pairs)
}

stratigraphic_compatibility <- function(ch1, ch2, st1, st2, ranges, UNKNOWN=-11,INAP=-22)	{
	
compats <- stratcompat <- divergstratcompat <- hierstratcompat <- 0
# first, find all pairs of characters
found_pairs <- vector(length=(st1*st2))
#pairs <- matrix(-1,st1*st2,2)
fnd_prs <- 0
pairs <- reduce_to_unique_pairs(ch1, ch2, st1, st2, UNKNOWN, INAP)
#pairs <- unique(cbind(ch1,ch2))
fnd_prs <- length(pairs[,1])
# now, find stratigraphic ranges of taxa with those pairs
pair_ranges <- date_character_pairs(ranges, pairs, ch1, ch2, st1, st2, UNKNOWN, INAP)

testprv <- vector(length=3)
testpr <- matrix(0,3,2)
tested <- 0
for (a in 1:(fnd_prs-2))	{
	testpr[1,] <- pairs[a,]
	testprv[1] <- 1+(testpr[1,1]*st2)+testpr[1,2]
	for (b in (a+1):(fnd_prs-1))	{
		testpr[2,] <- pairs[b,]
		testprv[2] <- 1+(testpr[2,1]*st2)+testpr[2,2]
		for (c in (b+1):fnd_prs)	{
			testpr[3,] <- pairs[c,]	
			testprv[3] <- 1+(testpr[3,1]*st2)+testpr[3,2]
			if (length(unique(testpr[,1]))==2 && length(unique(testpr[,2]))==2)	{
				compats <- compats+1
				swp <- find_swing_pairs(testpr,unique(testpr[,1]),unique(testpr[,2]))	# find the pair linking the other two pairs
				swppr <- 1+(swp[1]*st2)+swp[2]	# this gives the number to match with pair-ranges from all possible pairs
#				swingpr <- match(swppr,testprv)	# this gives which of the three considered pairs is the swing pair
				others <- subset(testprv,testprv!=swppr)	# vector of the two pairs (numbered from all pairs) other than the swing
				if (ranges[swppr,1]<=ranges[others[1],1] || ranges[swppr,1]<=ranges[others[2],1])	{
					stratcompat <- stratcompat+1
					if (ranges[swppr,1]<ranges[others[1],1] && ranges[swppr,1]<ranges[others[2],1])	{
						divergstratcompat <- divergstratcompat+1
						}  else if (ranges[swppr,1]>ranges[others[1],1] || ranges[swppr,1]>ranges[others[2],1])	{
						hierstratcompat <- hierstratcompat+1
						}   else if (ranges[swppr,1]==ranges[others[1],1] || ranges[swppr,1]==ranges[others[2],1])	{
						divergstratcompat <- divergstratcompat+0.5
						hierstratcompat <- hierstratcompat+0.5
						}	# end possible types of stratigraphic compatibility
					}	# end case of stratigraphic compatibility
#				tested <- tested+1
#				if (tested==1)	{
#					tested_pairs <- c(a,b,c)
#					} else tested_pairs <- rbind(tested_pairs,c(a,b,c))
				}	# end case where we have a compatible pair with 3 combinations
			} # go through third pairs
		} # go through second pairs
	} # go through first pairs

return (c(compats,stratcompat,divergstratcompat,hierstratcompat))
}

stratigraphic_compatibility_simple <- function(ch1, ch2, st1, st2, ranges, UNKNOWN, INAP)	{
# first, find all pairs of characters
#found_pairs <- vector(length=(st1*st2))
#pairs <- matrix(-1,st1*st2,2)
fnd_prs <- 0
pairs <- reduce_to_unique_pairs(ch1, ch2, st1, st2, UNKNOWN, INAP)

stratcompat <- poss_strat_compat <- 0
#compats <- stratcompat <- divergstratcompat <- hierstratcompat <- 0
if (is.matrix(pairs) && nrow(pairs)>2)	{
	fnd_prs <- length(pairs[,1])
#	} else	{
#	fnd_prs <- 1
#	}
#if (nrow(pairs)>2)	{
	# now, find stratigraphic ranges of taxa with those pairs
	pair_ranges <- date_character_pairs(ranges, pairs, ch1, ch2, st1, st2, UNKNOWN, INAP)
	testprv <- vector(length=3)
	testpr <- matrix(0,3,2)
	for (a in 1:(fnd_prs-2))	{
		testpr[1,] <- pairs[a,]
#		testprv[1] <- 1+(testpr[1,1]*st2)+testpr[1,2]	# get unique number of 1st pair
		for (b in (a+1):(fnd_prs-1))	{
			testpr[2,] <- pairs[b,]
#			testprv[2] <- 1+(testpr[2,1]*st2)+testpr[2,2]	# get unique number of 2nd pair
			for (c in (b+1):fnd_prs)	{
				testpr[3,] <- pairs[c,]
#			testprv[3] <- 1+(testpr[3,1]*st2)+testpr[3,2]	# get unique number of 3rd pair
			# compare if the three have two states each.
				if (length(unique(testpr[,1]))==2 && length(unique(testpr[,2]))==2)	{
				#print(paste("doing ",paste(a,b,c,sep=" "),sep=""))
				#print(testpr)
					poss_strat_compat <- poss_strat_compat+1
					swp <- find_swing_pairs(testpr,unique(testpr[,1]),unique(testpr[,2]))	# find the pair linking the other two pairs
					swppr <- c(a,b,c)[c(a,b,c)*(testpr[,1] %in% swp[1])*(testpr[,2] %in% swp[2])>0]
					others <- c(a,b,c)[c(a,b,c)*(testpr[,1] %in% swp[1])*(testpr[,2] %in% swp[2])==0]
					if (pair_ranges[swppr,1]<=pair_ranges[others[1],1] || pair_ranges[swppr,1]<=pair_ranges[others[2],1])	{
						stratcompat <- stratcompat+1	# end case of stratigraphic compatibility
					### move this stuff to full routine above.
#					if (pair_ranges[swppr,1]<pair_ranges[others[1],1] && pair_ranges[swppr,1]<pair_ranges[others[2],1])	{
#						divergstratcompat <- divergstratcompat+1
#						}  else if (pair_ranges[swppr,1]>pair_ranges[others[1],1] || pair_ranges[swppr,1]>pair_ranges[others[2],1])	{
#						hierstratcompat <- hierstratcompat+1
#						}   else if (pair_ranges[swppr,1]==pair_ranges[others[1],1] || pair_ranges[swppr,1]==pair_ranges[others[2],1])	{
#						divergstratcompat <- divergstratcompat+0.5
#						hierstratcompat <- hierstratcompat+0.5
#						}	# end possible types of stratigraphic compatibility
						}
					}	# end case where we have a compatible pair with 3 combinations
				} # go through third pairs
			} # go through second pairs
		} # go through first pairs
	}
return(c(stratcompat,poss_strat_compat))
}

stratigraphic_compatibility_simple_effed <- function(ch1, ch2, st1, st2, ranges, UNKNOWN=-11,INAP=-22)	{
# first, find all pairs of characters
#found_pairs <- vector(length=(st1*st2))
#pairs <- matrix(-1,st1*st2,2)
fnd_prs <- 0
pairs <- reduce_to_unique_pairs(ch1, ch2, st1, st2, UNKNOWN, INAP)

stratcompat <- poss_strat_compat <- 0
#compats <- stratcompat <- divergstratcompat <- hierstratcompat <- 0
if (is.matrix(pairs) && nrow(pairs)>2)	{
	fnd_prs <- length(pairs[,1])
	# now, find stratigraphic ranges of taxa with those pairs
	pair_ranges <- date_character_pairs(ranges, pairs, ch1, ch2, st1, st2, UNKNOWN, INAP)
	testprv <- vector(length=3)
	testpr <- matrix(0,3,2);
#	if (nrow(pair_ranges)<3)	{
		
#		}
	for (a in 1:(fnd_prs-2))	{
		testpr[1,] <- pairs[a,]
#		testprv[1] <- 1+(testpr[1,1]*st2)+testpr[1,2]	# get unique number of 1st pair
		for (b in (a+1):(fnd_prs-1))	{
			testpr[2,] <- pairs[b,]
#			testprv[2] <- 1+(testpr[2,1]*st2)+testpr[2,2]	# get unique number of 2nd pair
			for (c in (b+1):fnd_prs)	{
				testpr[3,] <- pairs[c,]
#			testprv[3] <- 1+(testpr[3,1]*st2)+testpr[3,2]	# get unique number of 3rd pair
			# compare if the three have two states each.
				if (length(unique(testpr[,1]))==2 && length(unique(testpr[,2]))==2)	{
					poss_strat_compat <- poss_strat_compat+1
					swp <- find_swing_pairs(testpr,unique(testpr[,1]),unique(testpr[,2]))	# find the pair linking the other two pairs
					swppr <- c(a,b,c)[c(a,b,c)*(testpr[,1] %in% swp[1])*(testpr[,2] %in% swp[2])>0]
					others <- c(a,b,c)[c(a,b,c)*(testpr[,1] %in% swp[1])*(testpr[,2] %in% swp[2])==0]
					if (pair_ranges[swppr,1]<=pair_ranges[others[1],1] || pair_ranges[swppr,1]<=pair_ranges[others[2],1])	{
						stratcompat <- stratcompat+1	# end case of stratigraphic compatibility
						}
					}	# end case where we have a compatible pair with 3 combinations
				} # go through third pairs
			} # go through second pairs
		} # go through first pairs
	} else	{
	stratcompat <- 1;	# end case of stratigraphic compatibility
	poss_strat_compat <- 1;
	}	# case where we have only 1 or 2 of 4 pairs
return(c(stratcompat,poss_strat_compat))
}

stratigraphic_incompatibility_simple <- function(ch1, ch2, st1, st2, ranges, UNKNOWN=-11,INAP=-22)	{
# first, find all pairs of characters
#found_pairs <- vector(length=(st1*st2))
#pairs <- matrix(-1,st1*st2,2)
fnd_prs <- 0
pairs <- reduce_to_unique_pairs(ch1, ch2, st1, st2, UNKNOWN, INAP)

stratcompat <- poss_strat_compat <- 0
#compats <- stratcompat <- divergstratcompat <- hierstratcompat <- 0
if (is.matrix(pairs) && nrow(pairs)>2)	{
	fnd_prs <- length(pairs[,1])
	# now, find stratigraphic ranges of taxa with those pairs
	pair_ranges <- date_character_pairs(ranges, pairs, ch1, ch2, st1, st2, UNKNOWN, INAP)
	testprv <- vector(length=3)
	testpr <- matrix(0,3,2);
#	if (nrow(pair_ranges)<3)	{
		
#		}
	for (a in 1:(fnd_prs-2))	{
		testpr[1,] <- pairs[a,]
#		testprv[1] <- 1+(testpr[1,1]*st2)+testpr[1,2]	# get unique number of 1st pair
		for (b in (a+1):(fnd_prs-1))	{
			testpr[2,] <- pairs[b,]
#			testprv[2] <- 1+(testpr[2,1]*st2)+testpr[2,2]	# get unique number of 2nd pair
			for (c in (b+1):fnd_prs)	{
				testpr[3,] <- pairs[c,]
#			testprv[3] <- 1+(testpr[3,1]*st2)+testpr[3,2]	# get unique number of 3rd pair
			# compare if the three have two states each.
				if (length(unique(testpr[,1]))==2 && length(unique(testpr[,2]))==2)	{
					poss_strat_compat <- poss_strat_compat+1
					swp <- find_swing_pairs(testpr,unique(testpr[,1]),unique(testpr[,2]))	# find the pair linking the other two pairs
					swppr <- c(a,b,c)[c(a,b,c)*(testpr[,1] %in% swp[1])*(testpr[,2] %in% swp[2])>0]
					others <- c(a,b,c)[c(a,b,c)*(testpr[,1] %in% swp[1])*(testpr[,2] %in% swp[2])==0]
					if (pair_ranges[swppr,1]<=pair_ranges[others[1],1] || pair_ranges[swppr,1]<=pair_ranges[others[2],1])	{
						stratcompat <- stratcompat+1	# end case of stratigraphic compatibility
						}
					}	# end case where we have a compatible pair with 3 combinations
				} # go through third pairs
			} # go through second pairs
		} # go through first pairs
	} else	{
	stratcompat <- 1;	# end case of stratigraphic compatibility
	poss_strat_compat <- 1;
	}	# case where we have only 1 or 2 of 4 pairs
return(c(stratcompat,poss_strat_compat))
}

stratigraphic_compatibility_matrix <- function(chmatrix,states,types,ranges,UNKNOWN=-11,INAP=-22)	{
comp_matrix <- compatibility_matrix(chmatrix,states,types,UNKNOWN,INAP)
nchars <- ncol(chmatrix)
notu <- nrow(chmatrix)
strat_compat_matrix <- data.frame(matrix(0,nchars,nchars))
for (c1 in 1:(nchars-1))	{
	ch1 <- chmatrix[,c1];
	st1 <- states[c1];
	# get remaining characters
	cc <- ((c1+1):nchars)[comp_matrix[c1,(c1+1):nchars]>0];
	cb <- 1;
	coded1 <- (1:nchars)[ch1!=UNKNOWN][(1:nchars)[ch1!=UNKNOWN] %in% (1:nchars)[ch1!=INAP]];
	while (cb <= length(cc))	{
		c2 <- cc[cb]
		ch2 <- chmatrix[,c2]
		st2 <- states[c2]
		coded2 <- (1:nchars)[ch2!=UNKNOWN][(1:nchars)[ch2!=UNKNOWN] %in% (1:nchars)[ch2!=INAP]];
		
		if (sum(coded1 %in% coded2)>0)	{
			sc <- stratigraphic_compatibility_simple(ch1, ch2, st1, st2, ranges, UNKNOWN, INAP);
			}	else	{
			sc <- c(0,0);
			}
		
		if (sc[2]>0)	{
			strat_compat_matrix[c1,c2] <- strat_compat_matrix[c2,c1] <- sc[1]/sc[2];
			}	else	{
			strat_compat_matrix[c1,c2] <- strat_compat_matrix[c2,c1] <- 0.5;
			}
		cb <- cb+1;
		}
	}
#	for (c3 in (c1+1):nchars)	{
#		if(comp_matrix[c1,c2]==1)	{
#			ch2 <- chmatrix[,c2]
#			st2 <- states[c2]
#			}
#		}
return(strat_compat_matrix)
}

mutual_compatibility <- function(compat_matrix,denom="min")	{
# routine to calculate mutual compatibility, as described in O'Keefe & Wagner (2002: Syst. Biol.).
# if denom == all, then divide mutuals by total characters
# if denom == max, then divide mutuals by maximum compatibilities of two characters
# if denom == min, then divide mutuals by minimum compatibilities of two characters
nchars <- nrow(compat_matrix)
mutsim <- matrix(1,nchars,nchars)
compsim <- matrix(0,nchars,nchars)
char_comps <- rowSums(compat_matrix)-1
mutual <- c()
for (c1 in 1:nchars)	{
	if (denom=="all") {
		mutsim[c1,c1] <- char_comps[c1]/(nchars-1)
		} else	{
		mutsim[c1,c1] <- 1.0
		}
	compsim[c1,c1] <- char_comps[c1]
	}
for (c1 in 1:(nchars-1))	{
	for (c2 in (c1+1):nchars)	{
		m <- 0
		for (c3 in 1:nchars)	{
			if (c3!=c1 && c3!=c2)	{
				m <- m+(compat_matrix[c1,c3]*compat_matrix[c2,c3])
				}
			}
		if (char_comps[c1]>=char_comps[c2])	{
			pair <- c(c1,c2)
			} else	{
			pair <- c(c2,c1)
			}
		mutual <- rbind(mutual,c(pair,compat_matrix[c1,c2],max(char_comps[c1],char_comps[c2]),min(char_comps[c1],char_comps[c2]),m))
		if (c1==1 && c2==3)	colnames(mutual) <- c("ch1","ch2","compatible","compatibility_ch1","compatibility_ch2","mutual_compatibilities")
		compsim[c1,c2] <- compsim[c2,c1] <- m
		if (denom=="all") {
			mutsim[c1,c2] <- mutsim[c2,c1] <- m/(nchars-1)
			} else if (denom=="max") {
			mutsim[c1,c2] <- mutsim[c2,c1] <- m/max(char_comps[c1],char_comps[c2])
			} else if (denom=="min")	{
			mutsim[c1,c2] <- mutsim[c2,c1] <- m/min(char_comps[c1],char_comps[c2])
			}
		}
	}
robin <- list(mutual,compsim,mutsim)
names(robin) <- c("Mutual_Compatibility","Compatibility_Similarity","Mutual_Similarity")
return(robin)
}