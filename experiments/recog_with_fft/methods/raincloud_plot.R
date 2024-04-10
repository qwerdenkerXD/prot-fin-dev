args <- commandArgs(trailingOnly = TRUE)

if ((length(args) - 1) %% 2 != 0 | (length(args) - 1) < 2) {
    message("Usage: Rscript raincloud_plot.R (<value-title> <one-line-values.csv>)... <out-file.png>\n")
    message("Set title with environment variable TITLE and label of x axis with X_LABEL\n")
    message("Example:")
    message("X_LABEL=Lengths TITLE=\"Length Distribution\" Rscript raincloud_plot.R \"Case 1\" <(echo 1,2,3) \"Case 2\" <(echo 2,1,2,3) out-file.png")
    q()
}

# Load necessary libraries
library(ggplot2)
library(ggdist)
library(tibble)

# Read the FASTA file
value_titles <- args[seq(1, length(args) - 1, 2)]
values <- args[seq(2, length(args) - 1, 2)]
values <- lapply(values, function(v) {
  sapply(read.csv(v, header=FALSE, nrows=1, sep=","), function(v) {as.numeric(v) + 1}) # pseudo count
})
groups <- lapply(seq_along(values), function(i) {tibble(value = values[[i]], group = value_titles[i])})

data <- dplyr::bind_rows(groups)

plt <- ggplot(data, aes(x = group, y = value, fill = group)) +
  geom_jitter(colour = "darkgray", position = position_jitter(width = 0.1, seed = 1), size = 0.1) +
  stat_halfeye(width = 2, position = position_nudge(x = 0.11, y=0), alpha = 0.8) +
  coord_flip() +
  theme(plot.title = element_text(hjust = 0.5), axis.ticks.y = element_blank(), legend.position = "none") +
  scale_y_continuous(n.breaks = 6, trans = "log10") +
  labs(x = "Frequencies", y = Sys.getenv("X_LABEL", "Values")) +
  ggtitle(Sys.getenv("TITLE", "Distribution of values"))

ggsave(args[length(args)], plot = plt, width = 8, height = 3, dpi = 300)
