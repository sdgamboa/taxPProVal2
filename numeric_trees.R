library(bugphyzz)
library(taxPPro)
library(purrr)
library(dplyr)
library(tidyr)
library(ape)
library(ggplot2)
library(ggtree)
library(ggtreeExtra)

bp <- importBugphyzz(v = 0, excludeRarely = FALSE)


df <- bp$aerophilicity
x <- df |>
    filter(Evidence == "asr")

df |>
    filter(Attribute_value == "aerotolerant")

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
            filter(totalN >= 50 & totalN <= 100) |>
            filter(notASR >= 10)
    })
    names(output) <- sub("_taxid$", "", taxRanks)
    output |>
        discard(~ nrow(.x) == 0)
}

tipDataCountPerRank <- map(tipDataAnnotated, countsPerRank) |>
    discard(~ !length(.x))

orders <- list_flatten(map(tipDataCountPerRank, ~ .x$order))
orderIdsL <- map(orders, ~ unique(.x$order_taxid))

getSubTr <- function(id) {
    idRx <- paste0("\\b", id, "\\b")
    targetClade <- grep(idRx, tr$node.label, value = TRUE)
    subTr <- extract.clade(tr, targetClade)
    return(subTr)
}

numberTips <- map(orderIdsL, \(x) {
    ntips <- map_int(x, \(y) Ntip(getSubTr(y)))
    names(ntips) <- x
    return(ntips)
})

lgList <- map(numberTips, ~ .x <= 100)

selectIDs <- map2(orderIdsL, lgList, ~ {
    .x[.y]
}) |>
    discard(~ !length(.x))

subTrees <- map(selectIDs, \(x) {
    subTree <- getSubTr(x[[1]])
    output <- list(subTree)
    names(output) <- x[[1]]
    output
})



## growth temperature -----------------------------------------------------
attr_name <- "coding genes"
t <- subTrees[[attr_name]][[1]]
taxid <- names(subTrees[[attr_name]])
tbl <- tipDataAnnotated[[attr_name]]
data <- filter(tbl, .data$order_taxid == .env$taxid) |>
    rename(id = tip_label)
myData <- data
subData <- select(myData, id, attr = Attribute_value)
p <- ggtree(t) %<+% myData +
    geom_tippoint(mapping = aes(color = Evidence)) +
    geom_tiplab(mapping = aes(label = Taxon_name), align = TRUE)
p2 <- ggtree::facet_plot(
    p = p,
    mapping = aes(x = attr, fill=Evidence),
    data = subData, geom = geom_bar, panel = attr_name,
    stat = "identity",
    orientation = 'y', width = 0.9
) +
    xlim_tree(0.2) +
    ggtitle("Order X") +
    theme_tree2()
p3 <- facet_widths(p2, widths = c(4, 1))

tail(t$tip.label)


tbl |>
    filter(tip_label == "g__60892")
taxid

