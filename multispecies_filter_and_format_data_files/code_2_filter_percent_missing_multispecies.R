# NOTE: Set the maximum percent missing data you will allow. Give this in percentage, not proportion.
max_percent_missing <- 20
min_w_minor_allele <- 3

# Read in simplified VCF file.
RAD <- as.matrix( read.table("output_from_step_one_multispecies.txt", 
	header = F, stringsAsFactor = F) )

# Clip off the first 10 columns of the matrix. These are the position and quality info, etc.
RAD2 <- RAD[ , 10:dim(RAD)[2] ]

# Save REF and ALT columns as vectors.
REF <- RAD[, 4]
ALT <- RAD[, 5]

# Save CHROM and POS columns as vector.
CHROM <- RAD[, 1]
POS <- RAD[, 2]

# Pull the first 3 characters out of the genotype strings (just the genotypes).

GT_mat <- RAD2[-1, ]
GT_nrow <- nrow(GT_mat)
GT_ncol <- ncol(GT_mat)
GT_mat_as_vector <- as.vector(GT_mat)
GT_mat_as_vector_missing <- which(GT_mat_as_vector == ".")

library(stringr)
GT_mat_as_vector[GT_mat_as_vector_missing] <- "."
GT_mat_as_vector_ONLY_GT <- str_sub(GT_mat_as_vector, 1, 3)
GT_mat_only_GT <- matrix(nrow = GT_nrow, ncol = GT_ncol, GT_mat_as_vector_ONLY_GT)

# Make names vector for value matching.
sample_names <- RAD2[1, ]

# Remove individuals with high missing data (10 > 90%, 1 ~= 50 %).
high_missing_data_indivs_and_dups <- read.table("high_missing_data_samples_and_duplicates.txt", stringsAsFactors = F)[, 1]
remove_for_high_missing_data_and_dup <- which( !(sample_names %in% high_missing_data_indivs_and_dups) )
GT_mat_only_GT <- GT_mat_only_GT[, remove_for_high_missing_data_and_dup]
sample_names <- sample_names[remove_for_high_missing_data_and_dup]

# Count the cells with missing data for each locus, and record the ..
# .. number of unique values to later test whether the locus is invariant.
# Also calculate the frequency with which the major allele is the only allele.
n_missing <- NA
n_unique <- NA
n_with_minor_allele_present <- NA
prop_with_minor_allele_present_non_missing <- NA
for(i in 1:dim(GT_mat_only_GT)[1])
{
	row_i <- GT_mat_only_GT[i, ]
	n_unique[i] <- length( unique(row_i) )
	n_missing[i] <- ifelse("." %in% row_i, length( which( row_i == "." ) ), 0)
	n_not_missing <- length(row_i) - n_missing[i]
	n_with_minor_allele_present[i] <- length( which( (row_i != "0/0") & (row_i != ".") ) )
	if( is.na(n_with_minor_allele_present[i]) )
	{
		n_with_minor_allele_present[i] <- 0 
	}
	prop_with_minor_allele_present_non_missing[i] <- n_with_minor_allele_present[i] / n_not_missing
}


# Calculate the proportion of samples missing data for each locus.
prop_missing <- n_missing / dim(GT_mat_only_GT)[2]

# Test whether the locus is invariant.
is_invariant <- (n_unique == 1) | ((n_unique == 2) & (n_missing > 0))

# Filter data down to just SNPs that have less than the specified percent missing data, that are not invariant ..
# .. and have at least the minimum minor allele frequency.
prop_max <- max_percent_missing / 100
keep <- which( (prop_missing <= prop_max) & (!is_invariant) & (n_with_minor_allele_present >= min_w_minor_allele) )
df <- GT_mat_only_GT[ keep,  ]
REF <- REF[ keep ]
ALT <- ALT[ keep ]
CHROM <- CHROM[ keep ]
POS <- POS[ keep ]
df <- cbind(CHROM, POS, REF, ALT, df)
colnames(df) <- c("CHROM","POS", "REF", "ALT", sample_names)

# Save the filtered data as a table.
write.table(file = "output_from_step_two_multispecies.txt",
	x = df, row.names = F, col.names = T, quote = F)

