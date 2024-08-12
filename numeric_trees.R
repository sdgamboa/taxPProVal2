library(bugphyzz)
library(taxPPro)
library(purrr)
library(dplyr)
library(tidyr)

bp <- importBugphyzz(v = 0)
bpModified <- bp |>
    map(~{
        .x |>
            mutate(NCBI_ID = as.character(NCBI_ID)) |>
            select(
                taxid = NCBI_ID, Attribute, Attribute_value, Evidence,
                Score, Validation
            ) |>
            as_tibble()
    }) |>
    keep(~ sum(.x$Evidence == "asr") > 0)


ltp <- ltp()
tr <- ltp$tree
tipData <- ltp$tip_data |>
    select(-accession) |>
    as_tibble()
nodeData <- ltp$node_data |>
    as_tibble()

dat <- bpModified$aerophilicity

tipDataAnnotated <- bpModified |>
    map(~{
        left_join(
            tipData, .x, by = "taxid", relationship = "many-to-many"
        )
    })

taxRanks <- c("class", "order", "family")
taxRanks <- paste0(taxRanks, "_taxid")

countsPerRank <- function(x) {
    # taxRanks <- c("class", "order", "family")
    # taxRanks <- paste0(taxRanks, "_taxid")
    output <- map(taxRanks, ~ {
        cols <- c("taxid", .x, "Evidence")
        df <- x[, cols, drop = FALSE]
        df |>
            drop_na() |>
            # filter(Evidence != "asr") |>
            count(!!sym(.x),  Evidence) |>
            group_by(!!sym(.x)) |>
            mutate(totalN = sum(n)) |>
            ungroup() |>
            mutate(tempCol = ifelse(Evidence == "asr", 0, n)) |>
            group_by(!!sym(.x)) |>
            mutate(notASR = sum(tempCol)) |>
            ungroup() |>
            select(-tempCol) |>
            filter(totalN >= 50 & totalN <= 100)
            # filter(notASR >= 40)
    })
    names(output) <- sub("_taxid$", "", taxRanks)
    output |>
        discard(~ nrow(.x) == 0)
}

tipDataCountPerRank <- map(tipDataAnnotated, countsPerRank) |>
    discard(~ !length(.x))


map(y, names)





