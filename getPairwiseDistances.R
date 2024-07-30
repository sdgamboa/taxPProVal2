library(taxPPro)
library(bugphyzz)
library(castor)
library(dplyr)
library(purrr)

ltp <- ltp()
tree <- ltp$tree

bp <- importBugphyzz()

taxids <- map(bp, ~ {
    .x |>
        filter(Evidence != "asr") |>
        pull(NCBI_ID)
}) |>
    unlist() |>
    unique()

tip_data <- ltp$tip_data |>
    group_by(taxid) |>
    mutate(n = n()) |>
    ungroup() |>
    filter(n == 1) |>
    select(-n) |>
    filter(taxid %in% taxids)
tips <- tip_data$tip_label

tim <- system.time({
    distances <- castor::get_all_pairwise_distances(tree, tips)
})

rownames(distances) <- tips
colnames(distances) <- tips


distPos <- which(upper.tri(distances, diag = FALSE), arr.ind = TRUE)

node1 <- rownames(distances)[distPos[, 1]]
node2 <- colnames(distances)[distPos[, 2]]
distance <- distances[distPos]

# Create the data frame
df <- data.frame(node1, node2, distance)

vroom::vroom_write(
    x = df, file = "pairwise_distances.tsv", delim = "\t",
    num_threads = 10
)

