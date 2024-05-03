#######Gene Function Analysis

# Load required libraries
# Load required libraries
library(dplyr)
library(ggplot2)
library(knitr)
library(igraph)

######## Load csv files. These are the annotated data from RAST
RE22 <- read.csv("~/BIO539/Final Project/190893.139 - Vibrio coralliilyticus RE22-.csv")
BAA_450 <- read.csv("~/BIO539/Final Project/190893.140 - Vibrio coralliilyticus ATCC BAA-450.csv")

#Calculate the size of each feature
data_RE22 <- mutate(RE22, size = abs(stop - start))
data_BAA_450 <- mutate(BAA_450, size = abs(stop - start))

# Summary statistics
summary_stats_RE22 <- data_RE22 %>%
  summarise(mean_size = mean(size), min_size = min(size), max_size = max(size))

summary_stats_BAA_450 <- data_BAA_450 %>%
  summarise(mean_size = mean(size), min_size = min(size), max_size = max(size))

# Summary Stats into table
df <- data.frame(
  Name = c("Mean Size", "Minimum Size", "Maximum Size"),
  RE22 = c(973.8127, 73, 20552),
  BAA450 = c(933.7458, 73, 17090)
)

# Display the table using kable
kable(
  df,
  knitr.table.format = 'html',
  align = c("l","c", "c"),  # Align columns left, center, and right
  col.names = c("Name", "RE22", "ATCC BAA-450"),  # Rename columns
  row.names = FALSE  # Hide row names
)
######## Analyze the most common gene functions
common_functions_RE22<- RE22 %>%
  filter(type == "peg") %>%  # Filter for protein-encoding genes
  count(function., sort = TRUE) %>%
  top_n(10, n)  # Select the top 10 most common functions

common_functions_BAA_450 <- BAA_450 %>%
  filter(type == "peg") %>%  # Filter for protein-encoding genes
  count(function., sort = TRUE) %>%
  top_n(10, n)  # Select the top 10 most common functions

# Combine the data frames
combined_data <- rbind(mutate(common_functions_RE22, strain = "V. coralliilyticus RE22"),
                       mutate(common_functions_BAA_450, strain = "V. coralliilyticus ATCC BAA-450"))

# Create the ggplot
ggplot(combined_data, aes(x = reorder(function., n), y = n)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Top 10 most common gene functions",
       x = "Function", y = "Frequency") +
  facet_wrap(~ strain, scales = "free", ncol = 1) +
  theme(plot.title = element_text(hjust = 0.5),  # Center the title
        plot.margin = margin(t = 1, r = 4, b = 1, l = 1, unit = "cm"),  # Adjust the plot margins
        axis.text.y = element_text(size = 4))  # Set the size of the y-axis labels

#######Orientation and Location Analysis

# Orientation analysis
orientation_count <- RE22 %>%
  count(strand)

orientation_count2 <- BAA_450 %>%
  count(strand)

orientation_count_combined <- rbind(mutate(orientation_count, strain = "Vibrio coralliilyticus RE22"),
                                    mutate(orientation_count2, strain = "Vibrio coralliilyticus ATCC BAA-450"))

# Create ggplot
ggplot(orientation_count_combined, aes(x = strand, y = n)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Gene orientation distribution", x = "Strand", y = "Count") +
  facet_wrap(~ strain, scales = "free", ncol = 1)

# Clustering of Genes along Contigs (simple visualization)
# Clustering of Genes along Contigs (simple visualization)
RE22 %>%
  ggplot(aes(x = start, y = as.factor(contig_id), color = strand)) +
  geom_point() +
  labs(title = "Gene distribution along contigs - Vibrio coralliilyticus RE22", x = "Start position", y = "Contig")

BAA_450 %>%
  ggplot(aes(x = start, y = as.factor(contig_id), color = strand)) +
  geom_point() +
  labs(title = "Gene distribution along contigs - Vibrio coralliilyticus ATCC BAA-450", x = "Start position", y = "Contig")

###### Antimicrobial Resistance Genes (ARG)
# Load csv files containing annotations from RAST
RE22_annotations <- read.csv("~/BIO539/Final Project/190893.139 - Vibrio coralliilyticus RE22-.csv")
BAA_450_annotations <- read.csv("~/BIO539/Final Project/190893.140 - Vibrio coralliilyticus ATCC BAA-450.csv")

# Filter data based on a simple keyword search for "resistance" (resistance genes)
RE22_resistance_genes <- subset(RE22_annotations, grepl("resistance", function. , ignore.case = TRUE))
BAA_450_resistance_genes <- subset(BAA_450_annotations, grepl("resistance", function. , ignore.case = TRUE))

# Identify the common and unique resistance genes.
common_genes <- intersect(RE22_resistance_genes$function., BAA_450_resistance_genes$function.)
unique_RE22 <- setdiff(RE22_resistance_genes$function., BAA_450_resistance_genes$function.)
unique_BAA_450 <- setdiff(BAA_450_resistance_genes$function., RE22_resistance_genes$function.)

# Combine the data frames
combined__resistance_data <- rbind(mutate(RE22_resistance_genes, strain = "Vibrio coralliilyticus RE22"),
                                   mutate(BAA_450_resistance_genes, strain = "Vibrio coralliilyticus ATCC BAA-450")) %>% 
  select(contig_id, feature_id, location, start, stop, strand, function., nucleotide_sequence, aa_sequence, strain)

#Save output as csv File
write.csv(combined__resistance_data, "~/BIO539/Final Project/combined__resistance_data.csv")

# Output results
list_of_genes <- list(common = common_genes, uniqueToRE22 = unique_RE22, uniqueToBAA_450 = unique_BAA_450)

all_genes <- unique(c(RE22_annotations$function., BAA_450_annotations$function.))

# Create a matrix showing presence (1) or absence (0)
presence_absence_matrix <- data.frame(
  Gene = all_genes,
  RE22 = as.integer(all_genes %in% RE22_annotations$function.),
  BAA_450 = as.integer(all_genes %in% BAA_450_annotations$function.)
)

# Genes present in both strains
common_genes <- presence_absence_matrix$Gene[presence_absence_matrix$RE22 == 1 & presence_absence_matrix$BAA_450 == 1]

# Genes unique to each strain
unique_RE22 <- presence_absence_matrix$Gene[presence_absence_matrix$RE22 == 1 & presence_absence_matrix$BAA_450 == 0]
unique_BAA_450 <- presence_absence_matrix$Gene[presence_absence_matrix$RE22 == 0 & presence_absence_matrix$BAA_450 == 1]

# save output as csv file
write.csv(presence_absence_matrix, "~/BIO539/Final Project/presence_absence_matrix.csv")

# Number of unique and common genes
data <- data.frame(
  Category = c("Unique to RE22", "Unique to BAA_450"),
  Count = c(length(unique_RE22), length(unique_BAA_450))
)

# Bar plot
ggplot(data, aes(x = Category, y = Count)) +
  geom_bar(stat = "identity") 
theme_minimal() +
  labs(title = "Comparison of ARGs between Two Bacterial Strains", x = "Strain", y = "Number of Genes")
######Presence and absent matrix specific to ARGs
presence_absence_matrix_ARGs <- subset(presence_absence_matrix, grepl("resistance", Gene, ignore.case = TRUE))

# Genes present in both strains
common_genes_arg <- presence_absence_matrix_ARGs$Gene[presence_absence_matrix_ARGs$RE22 == 1 & presence_absence_matrix_ARGs$BAA_450 == 1]

# Genes unique to each strain
unique_RE22_arg <- presence_absence_matrix_ARGs$Gene[presence_absence_matrix_ARGs$RE22 == 1 & presence_absence_matrix_ARGs$BAA_450 == 0]
unique_BAA_450_arg <- presence_absence_matrix_ARGs$Gene[presence_absence_matrix_ARGs$RE22 == 0 & presence_absence_matrix_ARGs$BAA_450 == 1]


# Number of unique and common ARGs
data_arg <- data.frame(
  Category = c("Common Genes", "Unique to RE22", "Unique to BAA_450"),
  Count = c(length(common_genes_arg), length(unique_RE22_arg), length(unique_BAA_450_arg))
)

# Bar plot
#ggplot(data_arg, aes(x = Category, y = Count)) +
#  geom_bar(stat = "identity") 
#  theme_minimal() +
#  labs(title = "Comparison of ARGs between Two Bacterial Strains", x = "Strain", y = "Number of Genes")

# Create edges and nodes for the graph based on gene presence
edges <- data.frame(
  from = c(rep("RE22", length(unique_RE22_arg)), rep("BAA_450", length(unique_BAA_450_arg))),
  to = c(unique_RE22_arg, unique_BAA_450_arg)
)


# Add common genes connections
edges <- rbind(edges, data.frame(from = "RE22", to = common_genes_arg))
edges <- rbind(edges, data.frame(from = "BAA_450", to = common_genes_arg))

# Create graph
graph <- graph_from_data_frame(edges, directed = FALSE)

# Plot
plot(graph, vertex.size = 10, vertex.label.cex = 0.5, edge.arrow.size = 5, main = "Network of ARGs Shared and Unique to Strains")
