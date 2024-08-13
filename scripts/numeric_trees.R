library(bugphyzz)
library(taxPPro)
library(purrr)
library(dplyr)
library(tidyr)
library(ape)
library(ggplot2)
library(ggtree)
library(ggtreeExtra)

## Bugphyzz data
bp <- importBugphyzz(v = 0, excludeRarely = FALSE)
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

## Tree data
ltp <- ltp()
tr <- ltp$tree
tipData <- ltp$tip_data |>
    select(-accession) |>
    as_tibble()
nodeData <- ltp$node_data |>
    as_tibble()

## Annotate tips
tipDataAnnotated <- bpModified |>
    map(~{
        left_join(
            tipData, .x, by = "taxid", relationship = "many-to-many"
        )
    })

taxRanks <- c("class", "order", "family")
taxRanks <- paste0(taxRanks, "_taxid")

countsPerRank <- function(x) {
    ## Function for generating counts per clade.
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

## Order
order <- list_flatten(map(tipDataCountPerRank, ~ .x$order)) |>
    discard(~ !length(.x))
orderIDL <- map(order, ~ unique(.x$order_taxid))

getSubTr <- function(id) {
    ## Function to get a subtree based on an internal node taxid.
    idRx <- paste0("\\b", id, "\\b")
    targetClade <- grep(idRx, tr$node.label, value = TRUE)
    subTr <- extract.clade(tr, targetClade)
    return(subTr)
}

numberTips <- map(orderIDL, \(x) {
    ntips <- map_int(x, \(y) Ntip(getSubTr(y)))
    names(ntips) <- x
    return(ntips)
})

## A 100 tips should be easy to view.
## Work on number of tips of the subtree, not the number of rows
lgL <- map(numberTips, ~ .x <= 100)

selectIDs <- map2(orderIDL, lgL, ~ {
    .x[.y]
}) |>
    discard(~ !length(.x))

## Get a subtree for each taxid
subTrees <- map(selectIDs, \(x) {
    l <- map(x,\(y) getSubTr(y))
    names(l) <- x
    l
})


## There are five numeric attributes that we should check in detail
## 1. Coding genes.
## 2. Genome size.
## 3. Growth temperature.
## 4. Width.
## 5. Length.

## Let's generate a few subtrees (clades) to explore the results
## The subtrees shoould correspond to clades based on the taxonomic
## classification; however, since phylogeny and taxonomy not always correspond
## to each other, some discrepancies are expected. Still, the visualization
## of these trees should be helpful.
plotCladeNum <- function(attrName, treeIndex = 1, cladeName = NULL) {
    clade <- subTrees[[attrName]][[treeIndex]]
    taxid <- names(subTrees[[attrName]])[[treeIndex]]
    ann <- tipDataAnnotated[[attrName]] |>
        filter(.data[["order_taxid"]] == .env[["taxid"]]) |>
        rename(id = tip_label)
    annSubset <- select(ann, id, attr = Attribute_value)
    annSubset # for the bar plots, a different dataset is needed
    p <- ggtree(clade) %<+% ann +
        geom_tippoint(mapping = aes(color = Evidence)) +
        geom_tiplab(mapping = aes(label = taxid), align = TRUE)
    p2 <- ggtree::facet_plot(
        p = p,
        mapping = aes(x = attr, fill=Evidence),
        data = annSubset, geom = geom_bar, panel = attrName,
        stat = "identity",
        orientation = 'y', width = 0.9
    ) +
        xlim_tree(0.2) +
        theme_tree2()
    if (is.null(cladeName)) {
        p2 <- p2 +
            ggtitle(paste0(attrName, " - ", taxid))
    } else {
        p2 <- p2 +
            ggtitle(paste0(attrName, " - ", cladeName))

    }

}

## Growth temperature
gtPlot <- plotCladeNum("growth temperature", cladeName = "Deinococcales")
gtPlot <- facet_widths(gtPlot, widths = c(4, 1))
ggsave(
    filename = "numeric_clades_plots/growth_temperature_Deinococcales.png",
    plot = gtPlot, width = 20, height = 15
)

## optimal ph
opPlot <- plotCladeNum("optimal ph", cladeName = "Halobacteriales")
opPlot <- facet_widths(opPlot, widths = c(4, 1))
ggsave(
    filename = "numeric_clades_plots/optimal_ph_Halobacteriales.png",
    plot = opPlot, width = 20, height = 15
)

## coding genes
cgPlot <- plotCladeNum("coding genes", cladeName = "Entomoplasmatales")
cgPlot <- facet_widths(cgPlot, widths = c(4, 1))
ggsave(
    filename = "numeric_clades_plots/codinge_genes_Entoplasmatales.png",
    plot = cgPlot, width = 20, height = 15
)

## genome size
gzPlot <- plotCladeNum("genome size", cladeName = "Rickettsiales")
gzPlot <- facet_widths(gzPlot, widths = c(4, 1))
ggsave(
    filename = "numeric_clades_plots/genome_size_Rickettsiales.png",
    plot = gzPlot, width = 20, height = 15
)

## width
wPlot <- plotCladeNum("width", cladeName = "Desulfobacterales")
wPlot <- facet_widths(wPlot, widths = c(4, 1))
ggsave(
    filename = "numeric_clades_plots/width_Desulfobacterales.png",
    plot = wPlot, width = 20, height = 15
)

## length
lPlot <- plotCladeNum("length", cladeName = "Desulfobacterales")
lPlot <- facet_widths(lPlot, widths = c(4, 1))
ggsave(
    filename = "numeric_clades_plots/length_Desulfobacterales.png",
    plot = lPlot, width = 20, height = 15
)
