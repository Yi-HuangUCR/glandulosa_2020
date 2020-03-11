library(stringr)

# Read in simplified VCF file.
RAD <- read.table("output_from_step_two_2N.txt", 
	header = F, stringsAsFactor = F)
RAD <- as.matrix(RAD)

# Clip off the first 3 columns of the matrix. These are the CHROM, REF and ALT columns.
RAD2 <- as.matrix( RAD[ , 4:dim(RAD)[2] ] )

# Make names vector for value matching.
sample_names <- RAD2[1, ]

# Eliminate the first row (the sample names).
RAD3 <- RAD2[-1, ]

# Read in table of population/taxa data.
ref_table <- read.table("A_gland_samples.csv", sep = ",", header = T, stringsAsFactors = F)

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
REF <- REF[ -which(lose_these) ]
ALT <- ALT[ -which(lose_these) ]

# Toss out any SNPs that have an N listed for the REF or ALT alleles.
N_listed <- (REF == "N") | (ALT == "N")
keep_these <- which( !(N_listed) )
REF <- REF[keep_these]
ALT <- ALT[keep_these]
GT_mat <- GT_mat[keep_these, ]

# Remove SNPs that have more numbers in genotypes than entries in REF and ALT together.
has_a_0 <- apply(GT_mat, 1, FUN = function(x){ any(str_detect(x, pattern = "0")) } )
has_a_1 <- apply(GT_mat, 1, FUN = function(x){ any(str_detect(x, pattern = "1")) } )
has_a_2 <- apply(GT_mat, 1, FUN = function(x){ any(str_detect(x, pattern = "2")) } )
has_a_3 <- apply(GT_mat, 1, FUN = function(x){ any(str_detect(x, pattern = "3")) } )
okay_for_2 <- (has_a_0 & has_a_1 & has_a_2 & !has_a_3) & (nchar(ALT) == 3)
okay_for_3 <- (has_a_0 & has_a_1 & has_a_2 & has_a_3) & (nchar(ALT) ==  5)
okay_for_1 <- has_a_0 & has_a_1 & !has_a_2 & !has_a_3 & (nchar(ALT) ==  1)
no_other_characters <- apply( GT_mat, 1, FUN = function(x){ all( unlist( strsplit(x, "/") ) %in% c(".", "0", "1", "2", "3") ) } )
okay <- (okay_for_1 | okay_for_2 | okay_for_3) & no_other_characters
REF <- REF[which(okay)]
ALT <- ALT[which(okay)]
GT_mat <- GT_mat[which(okay), ]

# Read in the IUPAC base code table.
IUPAC_table <- read.table(file = "IUPAC_base_code_table.txt", header = T, stringsAsFactor = F)

# Replace any missing data with ?.
GT_nrow <- nrow(GT_mat)
GT_ncol <- ncol(GT_mat)
GT_mat_as_vector <- as.character(GT_mat)
if("." %in% GT_mat_as_vector)
{
	missing_at <- which(GT_mat_as_vector == ".")
	GT_mat_as_vector[ missing_at ] <- "?"
}

# Shove back into a matrix.
GT_mat <- matrix(nrow = GT_nrow, ncol = GT_ncol, GT_mat_as_vector)

# Loop through the rows (loci).
A <- function(x){ strsplit(x, "&") }
B <- function(x){ lapply(x, FUN = function(x){sort(unique(x))} ) }
C <- function(x){ lapply(x, FUN = function(x){paste0(x, collapse = "&")} ) }
reform <- function(x){ unlist( C( B( A(x) ) ) ) }
alt_nchar <- nchar(ALT)

for(i in 1:GT_nrow)
{
	locus_i <- as.character( GT_mat[i, ] )
	ref <- REF[i] 
	alt1 <- substr(ALT[i], 1, 1)
	locus_i <- str_replace_all(locus_i, "/", "&")
        locus_i <- str_replace_all(locus_i, "0", ref)
        locus_i <- str_replace_all(locus_i, "1", alt1)
	if(alt_nchar[i] == 3)
	{
		alt2 <- substr(ALT[i], 3, 3)		
        	locus_i <- str_replace_all(locus_i, "2", alt2)
	}
        if(alt_nchar[i] == 5)
        {
 		alt2 <- substr(ALT[i], 3, 3)
                locus_i <- str_replace_all(locus_i, "2", alt2)
                alt3 <- substr(ALT[i], 5, 5)
                locus_i <- str_replace_all(locus_i, "3", alt3)
		locus_i <- reform(locus_i)
       }
	locus_i <- reform(locus_i)
	GT_mat[i, ] <- IUPAC_table[ match(locus_i, IUPAC_table[, 2]), 1 ]
	if(any(is.na(GT_mat[i, ])))
	{
		print(paste("NA produced at:", i))
		print(locus_i)
		print(match(locus_i, IUPAC_table[, 2]))
	}
}

# Looped through the IUPACked matrix and process columns into nexus alignments.
for(i in 1:GT_ncol)
{
	if(i == 1)
	{
		nexus <- paste(GT_mat[, i], collapse = "")
	}else
	{
		nexus <- rbind( nexus, paste(GT_mat[, i], collapse = "") )
	}
}

# Create taxa labels for nexus format file.
subsp <- ref_table[match(sample_names, ref_table$sample), 3]

taxa_labels <- NA
seen_it <- "whoopdiddyscoop"
for(i in 1:length(subsp))
{
	if(subsp[i] %in% seen_it)
	{
		n_times_seen <- length( which(seen_it == subsp[i]) )
		taxa_labels[i] <- paste0(subsp[i], "_", n_times_seen + 1)
	}else
	{
		taxa_labels[i] <- paste0(subsp[i], "_1")
	}
	seen_it <- c(seen_it, subsp[i])
}	

# Calculate the proportion missing genotypes for each sample/taxon.
prop_missing_per_sample <- apply(nexus, 1, FUN = function(x){ return( str_count(x, "\\?") / nchar(x) ) })
sort_order <- rev( order(prop_missing_per_sample) )
Sample <- taxa_labels[sort_order]
Prop_missing <- prop_missing_per_sample[sort_order]
missing_report <- data.frame(Sample, Prop_missing)
write.table(file = "proportion_missing_table_2N.txt",
        x = missing_report,
        row.name = F, col.name = F, quote = F)

# Save table with nexus taxon label and sample code number.y
relate_table <- data.frame(taxa_labels, sample_names)
colnames(relate_table) <- c("Name_in_nex_phy", "Sample_name")
write.csv(file = "nexus_phy_taxa_label_ref_table_2N.csv",
        x = relate_table)

# Paste the individual tags ("taxa") onto the nexus alignment.
nexus <- matrix(nrow = length(nexus), ncol = 1,
	paste(taxa_labels, as.character(nexus) )
)

n_tax <- paste0("ntax=", dim(nexus)[1])
n_char <- paste0("nchar=", dim(GT_mat)[1], ";")
Dimensions <- paste("Dimensions", n_tax, n_char)
header <- rbind(
	"#NEXUS", 
	"Begin data;", 
	Dimensions, 
	"Format datatype=dna missing=? gap=-;",
	"Matrix")
footer <- rbind(
	";",
	"End;")
formatted_nexus <- rbind(
	header,
	nexus,
	footer
)

formatted_phylip <- rbind(
	paste(n_tax, n_char),
	nexus	
)

# Save the nexus format file.
write.table(file = "A_glandulosa_2N.nex",
	x = formatted_nexus,
	row.name = F, col.name = F, quote = F)

# Save the phylip format file.
write.table(file = "A_glandulosa_2N.phy",
        x = formatted_phylip,
        row.name = F, col.name = F, quote = F)
