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


# Numeric trees -----------------------------------------------------------

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
        geom_tippoint(mapping = aes(color = Evidence), size = 4) +
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
gtPlot <- plotCladeNum("growth temperature")
gtPlot <- facet_widths(gtPlot, widths = c(4, 1))
ggsave(
    filename = "numeric_clade_plots/growth_temperature.png",
    plot = gtPlot, width = 20, height = 15
)

## optimal ph
opPlot <- plotCladeNum("optimal ph")
opPlot <- facet_widths(opPlot, widths = c(4, 1))
ggsave(
    filename = "numeric_clade_plots/optimal_ph.png",
    plot = opPlot, width = 20, height = 15
)

## coding genes
cgPlot <- plotCladeNum("coding genes")
cgPlot <- facet_widths(cgPlot, widths = c(4, 1))
ggsave(
    filename = "numeric_clade_plots/codinge_genes.png",
    plot = cgPlot, width = 20, height = 15
)

## genome size
gzPlot <- plotCladeNum("genome size")
gzPlot <- facet_widths(gzPlot, widths = c(4, 1))
ggsave(
    filename = "numeric_clade_plots/genome_size.png",
    plot = gzPlot, width = 20, height = 15
)

## width
wPlot <- plotCladeNum("width")
wPlot <- facet_widths(wPlot, widths = c(4, 1))
ggsave(
    filename = "numeric_clade_plots/width.png",
    plot = wPlot, width = 20, height = 15
)

## length
lPlot <- plotCladeNum("length")
lPlot <- facet_widths(lPlot, widths = c(4, 1))
ggsave(
    filename = "numeric_clade_plots/length.png",
    plot = lPlot, width = 20, height = 15
)

# Discrete trees ----------------------------------------------------------

plotCladeDis <- function(attrName, treeIndex = 1, cladeName = NULL) {
    clade <- subTrees[[attrName]][[treeIndex]]
    taxid <- names(subTrees[[attrName]])[[treeIndex]]
    ann <- tipDataAnnotated[[attrName]] |>
        filter(.data[["order_taxid"]] == .env[["taxid"]]) |>
        rename(id = tip_label)
    treeAnn <- ann |>
        select(
            id, taxid, Taxon_name, Rank, NCBI_ID, Evidence
        )
    barPlotAnn <- ann |>
        select(
            id, Score, Attribute_value
        )
    p <- ggtree(clade) %<+% treeAnn +
        geom_tippoint(mapping = aes(color = Evidence), size = 4) +
        geom_tiplab(mapping = aes(label = taxid), align = TRUE)
    p2 <- ggtree::facet_plot(
        p = p,
        mapping = aes(x = Score, fill = Attribute_value),
        data = barPlotAnn, geom = geom_col, panel = attrName,
        orientation = 'y', width = 0.9
    ) +
        xlim_tree(0.2) +
        theme_tree2()
    p2
}

## Animal pathogen
# [1] "animal pathogen"     "host-associated"     "motility"            "plant pathogenicity"
# [5] "spore formation"     "aerophilicity"       "arrangement"         "biosafety level"
# [9] "gram stain"          "shape"

## Animal pathogen
apPlot <- plotCladeDis("animal pathogen")
apPlot <- facet_widths(apPlot, widths = c(4, 1))
ggsave(
    filename = "discrete_clade_plots/animal_pathogen.png",
    plot = apPlot, width = 20, height = 15
)

## Host-associated
haPlot <- plotCladeDis("host-associated")
haPlot <- facet_widths(haPlot, widths = c(4, 1))
ggsave(
    filename = "discrete_clade_plots/host_associated.png",
    plot = haPlot, width = 20, height = 15
)

## motility
moPlot <- plotCladeDis("motility")
moPlot <- facet_widths(moPlot, widths = c(4, 1))
ggsave(
    filename = "discrete_clade_plots/motility.png",
    plot = moPlot, width = 20, height = 15
)

## plant plant pathogenicity
ppPlot <- plotCladeDis("plant pathogenicity")
ppPlot <- facet_widths(ppPlot, widths = c(4, 1))
ggsave(
    filename = "discrete_clade_plots/plant_pathogenicity.png",
    plot = ppPlot, width = 20, height = 15
)

## spore formation
sfPlot <- plotCladeDis("spore formation")
sfPlot <- facet_widths(sfPlot, widths = c(4, 1))
ggsave(
    filename = "discrete_clade_plots/spore_formation.png",
    plot = sfPlot, width = 20, height = 15
)

## aerophilicity
aePlot <- plotCladeDis("aerophilicity")
aePlot <- facet_widths(aePlot, widths = c(4, 1))
ggsave(
    filename = "discrete_clade_plots/aerophilicity.png",
    plot = aePlot, width = 20, height = 15
)

## arrangement
arPlot <- plotCladeDis("arrangement")
arPlot <- facet_widths(arPlot, widths = c(4, 1))
ggsave(
    filename = "discrete_clade_plots/arrangement.png",
    plot = arPlot, width = 20, height = 15
)

## biosafety level
blPlot <- plotCladeDis("biosafety level")
blPlot <- facet_widths(blPlot, widths = c(4, 1))
ggsave(
    filename = "discrete_clade_plots/biosafety_level.png",
    plot = blPlot, width = 20, height = 15
)

## gram stain
gsPlot <- plotCladeDis("gram stain")
gsPlot <- facet_widths(gsPlot, widths = c(4, 1))
ggsave(
    filename = "discrete_clade_plots/gram_stain.png",
    plot = gsPlot, width = 20, height = 15
)

## shape
sPlot <- plotCladeDis("shape")
sPlot <- facet_widths(sPlot, widths = c(4, 1))
ggsave(
    filename = "discrete_clade_plots/shape.png",
    plot = sPlot, width = 20, height = 15
)


# Roseburia - Animal Pathogen ---------------------------------------------
## Lachnospirales        order 3085636

## These are annotations for the genus Rosenbura
## But these annotations could not be mapped to the LTP tree
"301301" %in% tipData$taxid
"166486" %in% tipData$taxid
"360807" %in% tipData$taxid

nodeName <- "186803"
sub_tree <- getSubTr(nodeName)

tipData |>
    filter(genus_taxid == "841")

tree_ann <- tipDataAnnotated[["animal pathogen"]] |>
    filter(.data[["family_taxid"]] == nodeName) |>
    rename(id = tip_label)
tree_ann1 <- tree_ann |>
    select(id, Evidence, taxid)

tree_ann2 <- tree_ann |>
    mutate(Attribute_value = as.character(Attribute_value)) |>
    select(id, Attribute_value, Score)

ggpt <- ggtree(sub_tree) %<+% tree_ann1 +
    geom_tippoint(mapping = aes(color = Evidence), size = 4) +
    geom_tiplab(mapping = aes(label = taxid), align = TRUE)
ggpt2 <- ggtree::facet_plot(
    p = ggpt,
    mapping = aes(x = Score, fill = Attribute_value),
    data = tree_ann2, geom = geom_col, panel = "Animal pathogen",
    orientation = 'y', width = 0.9
) +
    xlim_tree(0.23) +
    theme_tree2()

ggpt3 <- facet_widths(ggpt2, widths = c(6, 1))
ggsave(
    filename = "discrete_clade_plots/animal_pathogen_roseburia.pdf",
    plot = ggpt3, width = 20, height = 45
)
