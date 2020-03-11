library(stringr)

# Read in simplified VCF file.
RAD <- read.table("output_from_step_two_4N.txt", 
	header = F, stringsAsFactor = F)
RAD <- as.matrix(RAD)

# Clip off the first 2 columns of the matrix. These are the CHROM, REF and ALT columns.
RAD2 <- as.matrix( RAD[ , 4:dim(RAD)[2] ] )

# Make names vector for value matching.
sample_names <- RAD2[1, ]

# Eliminate the first row (the sample names).
RAD3 <- RAD2[-1, ]

# Read in table of population/taxa data.
ref_table <- read.table("A_gland_samples.csv", sep = ",", header = T, stringsAsFactors = F)

# Read in coordinate data.
coords <- read.csv("positions.csv")

# Make CHROM, REF and ALT vectors.
CHROM <- RAD[-1, 1]
REF <- RAD[-1, 2]
ALT <- RAD[-1, 3]

# Calculate the lengths of the REF strings. Use this to eliminate the ones with more ..
# .. than one base given for each. Also eliminate those with more than three ALT alleles.
REF_length <- nchar(REF)
longer_than_one_ref <- REF_length > 1
more_than_three_alt <- nchar(ALT) > 5
lose_these <- longer_than_one_ref | more_than_three_alt
GT_mat <- RAD3[-which(lose_these), ]
CHROM <- CHROM[ -which(lose_these) ]
REF <- REF[ -which(lose_these) ]
ALT <- ALT[ -which(lose_these) ]

# Toss out any SNPs that have an N listed for the REF or ALT alleles.
N_listed <- (REF == "N") | (ALT == "N")
keep_these <- which( !(N_listed) )
CHROM <- CHROM[keep_these]
REF <- REF[keep_these]
ALT <- ALT[keep_these]
GT_mat <- GT_mat[keep_these, ]

# Remove SNPs that have more numbers in genotypes than entries in REF and ALT together.
library(stringr)
has_a_0 <- apply(GT_mat, 1, FUN = function(x){ any(str_detect(x, pattern = "0")) } )
has_a_1 <- apply(GT_mat, 1, FUN = function(x){ any(str_detect(x, pattern = "1")) } )
has_a_2 <- apply(GT_mat, 1, FUN = function(x){ any(str_detect(x, pattern = "2")) } )
has_a_3 <- apply(GT_mat, 1, FUN = function(x){ any(str_detect(x, pattern = "3")) } )
okay_for_2 <- (has_a_0 & has_a_1 & has_a_2 & !has_a_3) & (nchar(ALT) == 3)
okay_for_3 <- (has_a_0 & has_a_1 & has_a_2 & has_a_3) & (nchar(ALT) ==  5)
okay_for_1 <- has_a_0 & has_a_1 & !has_a_2 & !has_a_3 & (nchar(ALT) ==  1)
no_other_characters <- apply( GT_mat, 1, FUN = function(x){ all( unlist( strsplit(x, "/") ) %in% c(".", "0", "1", "2", "3") ) } )
okay <- (okay_for_1 | okay_for_2 | okay_for_3) & no_other_characters

PASS <- ifelse(okay, "YES", "NO") 
REF <- REF[which(okay)]
ALT <- ALT[which(okay)]
GT_mat <- GT_mat[which(okay), ]

# Test if each RAD fragment entry is the same as the previous entry.
n_SNPs <- dim(GT_mat)[1]
last_one <- CHROM[1:(n_SNPs - 1)]
this_one <- CHROM[2:n_SNPs]
is_same_frag <- this_one == last_one

# Keep first entries per RAD fragment, and lose the first 10 columns.
keep_this <- c(TRUE, !(is_same_frag))
GT_mat <- GT_mat[keep_this, ]

# Reformat into STRUCTURE format.
n_indivs <- dim(GT_mat)[2]
for(i in 1:n_indivs)
{
        if(i == 1)
        {
                first <- substr(GT_mat[, i], 1, 1)
                which_missing <- which( first == "." )
                first[which_missing] <- "-9"
                second <- substr(GT_mat[, i], 3, 3)
                second[which_missing] <- "-9"
		third <- substr(GT_mat[, i], 5, 5)
		third[which_missing] <- "-9"
		fourth <- substr(GT_mat[, i], 7, 7)
		fourth[which_missing] <- "-9"
                formatted <- rbind(first, second, third, fourth)
        }else
        {
                first <- substr(GT_mat[, i], 1, 1)
                which_missing <- which( first == "." )
                first[which_missing] <- "-9"
                second <- substr(GT_mat[, i], 3, 3)
                second[which_missing] <- "-9"
		third <- substr(GT_mat[, i], 5, 5)
                third[which_missing] <- "-9"
                fourth <- substr(GT_mat[, i], 7, 7)
                fourth[which_missing] <- "-9"
                formatted <- rbind(formatted, first, second, third, fourth)
        }
}


# Make numeric IDs for samples.
samples_fac <- as.factor(sample_names)
samples_num <- as.numeric(samples_fac)

# Make subspecies numeric IDs.
subsp <- ref_table[match(sample_names, ref_table$sample), 3]
subsp_fac <- as.factor(subsp)
subsp_num <- as.numeric(subsp_fac)

# Make lat and long vectors.
latitude <- coords[match(sample_names, coords$sample_ID), 2]
longitude <- coords[match(sample_names, coords$sample_ID), 3]

# Attach numeric sample IDs and subspecies to the STRUCTURE format data.
individual <- rep(samples_num, each = 4)
population <- rep(subsp_num, each = 4) 
final_formatted <- cbind(individual, population, formatted)

# Write STRUCTURE format file as a table.
write.table(file = "A_glandulosa_4N.str", x = final_formatted,
        quote = F, row.names = F, col.names = F)

# Write a table relating the sample names and subspecies determinations to the numeric IDs.
df <- data.frame(samples_num, sample_names, subsp_num, subsp, latitude, longitude)
write.table(file = "A_glandulosa_structure_key_table_4N.txt", x = df,
        quote = F, row.names = F, col.names = T)
