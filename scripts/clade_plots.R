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
        xlim_tree(0.23) +
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
gtTaxNames <- names(subTrees[["growth temperature"]])
for (i in seq_along(gtTaxNames)) {
    gtPlot <- plotCladeNum("growth temperature", treeIndex = i)
    gtPlot <- facet_widths(gtPlot, widths = c(4, 1))
    gtPlotName <- file.path(
        "numeric_clade_plots",
        paste0("growth_temperature_", gtTaxNames[i], ".png")
    )
    ggsave(filename = gtPlotName, plot = gtPlot, width = 20, height = 15)
}
# gtPlot <- plotCladeNum("growth temperature")
# gtPlot <- facet_widths(gtPlot, widths = c(4, 1))
# ggsave(
#     filename = "numeric_clade_plots/growth_temperature.png",
#     plot = gtPlot, width = 20, height = 15
# )

## optimal ph
opTaxNames <- names(subTrees[["optimal ph"]])
for (i in seq_along(opTaxNames)) {
    opPlot <- plotCladeNum("optimal ph", treeIndex = i)
    opPlot <- facet_widths(opPlot, widths = c(4, 1))
    opPlotName <- file.path(
        "numeric_clade_plots",
        paste0("optimal_ph_", opTaxNames[i], ".png")
    )
    ggsave(filename = opPlotName, plot = opPlot, width = 20, height = 15)
}
# opPlot <- plotCladeNum("optimal ph")
# opPlot <- facet_widths(opPlot, widths = c(4, 1))
# ggsave(
#     filename = "numeric_clade_plots/optimal_ph.png",
#     plot = opPlot, width = 20, height = 15
# )

## coding genes
cgTaxNames <- names(subTrees[["coding genes"]])
for (i in seq_along(cgTaxNames)) {
    cgPlot <- plotCladeNum("coding genes", treeIndex = i)
    cgPlot <- facet_widths(cgPlot, widths = c(4, 1))
    cgPlotName <- file.path(
        "numeric_clade_plots",
        paste0("coding_genes_", cgTaxNames[i], ".png")
    )
    ggsave(filename = cgPlotName, plot = cgPlot, width = 20, height = 15)
}
# cgPlot <- plotCladeNum("coding genes")
# cgPlot <- facet_widths(cgPlot, widths = c(4, 1))
# ggsave(
#     filename = "numeric_clade_plots/codinge_genes.png",
#     plot = cgPlot, width = 20, height = 15
# )

## genome size
gzTaxNames <- names(subTrees[["genome size"]])
for (i in seq_along(gzTaxNames)) {
    gzPlot <- plotCladeNum("genome size", treeIndex = i)
    gzPlot <- facet_widths(gzPlot, widths = c(4, 1))
    gzPlotName <- file.path(
        "numeric_clade_plots",
        paste0("genome_size_", gzTaxNames[i], ".png")
    )
    ggsave(filename = gzPlotName, plot = gzPlot, width = 20, height = 15)
}
# gzPlot <- plotCladeNum("genome size")
# gzPlot <- facet_widths(gzPlot, widths = c(4, 1))
# ggsave(
#     filename = "numeric_clade_plots/genome_size.png",
#     plot = gzPlot, width = 20, height = 15
# )

## width
wTaxNames <- names(subTrees[["width"]])
for (i in seq_along(wTaxNames)) {
    wPlot <- plotCladeNum("width", treeIndex = i)
    wPlot <- facet_widths(wPlot, widths = c(4, 1))
    wPlotName <- file.path(
        "numeric_clade_plots",
        paste0("width_", wTaxNames[i], ".png")
    )
    ggsave(filename = wPlotName, plot = wPlot, width = 20, height = 15)
}
# wPlot <- plotCladeNum("width")
# wPlot <- facet_widths(wPlot, widths = c(4, 1))
# ggsave(
#     filename = "numeric_clade_plots/width.png",
#     plot = wPlot, width = 20, height = 15
# )

## length
lTaxNames <- names(subTrees[["length"]])
for (i in seq_along(lTaxNames)) {
    lPlot <- plotCladeNum("length", treeIndex = i)
    lPlot <- facet_widths(lPlot, widths = c(4, 1))
    lPlotName <- file.path(
        "numeric_clade_plots",
        paste0("length_", lTaxNames[i], ".png")
    )
    ggsave(filename = lPlotName, plot = lPlot, width = 20, height = 15)
}
# lPlot <- plotCladeNum("length")
# lPlot <- facet_widths(lPlot, widths = c(4, 1))
# ggsave(
#     filename = "numeric_clade_plots/length.png",
#     plot = lPlot, width = 20, height = 15
# )

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
        xlim_tree(0.23) +
        theme_tree2()
    p2
}

## Animal pathogen
# [1] "animal pathogen"     "host-associated"     "motility"            "plant pathogenicity"
# [5] "spore formation"     "aerophilicity"       "arrangement"         "biosafety level"
# [9] "gram stain"          "shape"

## Animal pathogen
apTaxNames <- names(subTrees[["animal pathogen"]])
for (i in seq_along(apTaxNames)) {
    apPlot <- plotCladeDis("animal pathogen", treeIndex = i)
    apPlot <- facet_widths(apPlot, widths = c(4, 1))
    apPlotName <- file.path(
        "discrete_clade_plots",
        paste0("animal_pathogen_", apTaxNames[i], ".png")
    )
    ggsave(filename = apPlotName, plot = apPlot, width = 20, height = 15)
}
# apPlot <- plotCladeDis("animal pathogen")
# apPlot <- facet_widths(apPlot, widths = c(4, 1))
# ggsave(
#     filename = "discrete_clade_plots/animal_pathogen.png",
#     plot = apPlot, width = 20, height = 15
# )

## Host-associated
haTaxNames <- names(subTrees[["host-associated"]])
for (i in seq_along(haTaxNames)) {
    haPlot <- plotCladeDis("host-associated", treeIndex = i)
    haPlot <- facet_widths(haPlot, widths = c(4, 1))
    haPlotName <- file.path(
        "discrete_clade_plots",
        paste0("host_associated_", haTaxNames[i], ".png")
    )
    ggsave(filename = haPlotName, plot = haPlot, width = 20, height = 15)
}
# haPlot <- plotCladeDis("host-associated")
# haPlot <- facet_widths(haPlot, widths = c(4, 1))
# ggsave(
#     filename = "discrete_clade_plots/host_associated.png",
#     plot = haPlot, width = 20, height = 15
# )

## motility
moTaxNames <- names(subTrees[["motility"]])
for (i in seq_along(moTaxNames)) {
    moPlot <- plotCladeDis("motility", treeIndex = i)
    moPlot <- facet_widths(moPlot, widths = c(4, 1))
    moPlotName <- file.path(
        "discrete_clade_plots",
        paste0("motility_", moTaxNames[i], ".png")
    )
    ggsave(filename = moPlotName, plot = moPlot, width = 20, height = 15)
}
# moPlot <- plotCladeDis("motility")
# moPlot <- facet_widths(moPlot, widths = c(4, 1))
# ggsave(
#     filename = "discrete_clade_plots/motility.png",
#     plot = moPlot, width = 20, height = 15
# )

## plant plant pathogenicity
ppTaxNames <- names(subTrees[["plant pathogenicity"]])
for (i in seq_along(ppTaxNames)) {
    ppPlot <- plotCladeDis("plant pathogenicity", treeIndex = i)
    ppPlot <- facet_widths(ppPlot, widths = c(4, 1))
    ppPlotName <- file.path(
        "discrete_clade_plots",
        paste0("plant_pathogenicity_", ppTaxNames[i], ".png")
    )
    ggsave(filename = ppPlotName, plot = ppPlot, width = 20, height = 15)
}
# ppPlot <- plotCladeDis("plant pathogenicity")
# ppPlot <- facet_widths(ppPlot, widths = c(4, 1))
# ggsave(
#     filename = "discrete_clade_plots/plant_pathogenicity.png",
#     plot = ppPlot, width = 20, height = 15
# )

## spore formation
sfTaxNames <- names(subTrees[["spore formation"]])
for (i in seq_along(sfTaxNames)) {
    sfPlot <- plotCladeDis("spore formation", treeIndex = i)
    sfPlot <- facet_widths(sfPlot, widths = c(4, 1))
    sfPlotName <- file.path(
        "discrete_clade_plots",
        paste0("spore_formation_", sfTaxNames[i], ".png")
    )
    ggsave(filename = sfPlotName, plot = sfPlot, width = 20, height = 15)
}
# sfPlot <- plotCladeDis("spore formation")
# sfPlot <- facet_widths(sfPlot, widths = c(4, 1))
# ggsave(
#     filename = "discrete_clade_plots/spore_formation.png",
#     plot = sfPlot, width = 20, height = 15
# )

## aerophilicity
aeTaxNames <- names(subTrees[["aerophilicity"]])
for (i in seq_along(aeTaxNames)) {
    aePlot <- plotCladeDis("aerophilicity", treeIndex = i)
    aePlot <- facet_widths(aePlot, widths = c(4, 1))
    aePlotName <- file.path(
        "discrete_clade_plots",
        paste0("aerophilicity_", aeTaxNames[i], ".png")
    )
    ggsave(filename = aePlotName, plot = aePlot, width = 20, height = 15)
}
# aePlot <- plotCladeDis("aerophilicity")
# aePlot <- facet_widths(aePlot, widths = c(4, 1))
# ggsave(
#     filename = "discrete_clade_plots/aerophilicity.png",
#     plot = aePlot, width = 20, height = 15
# )

## arrangement
arTaxNames <- names(subTrees[["arrangement"]])
for (i in seq_along(arTaxNames)) {
    arPlot <- plotCladeDis("arrangement", treeIndex = i)
    arPlot <- facet_widths(arPlot, widths = c(4, 1))
    arPlotName <- file.path(
        "discrete_clade_plots",
        paste0("arrangement_", arTaxNames[i], ".png")
    )
    ggsave(filename = arPlotName, plot = arPlot, width = 20, height = 15)
}
# arPlot <- plotCladeDis("arrangement")
# arPlot <- facet_widths(arPlot, widths = c(4, 1))
# ggsave(
#     filename = "discrete_clade_plots/arrangement.png",
#     plot = arPlot, width = 20, height = 15
# )

## biosafety level
blTaxNames <- names(subTrees[["biosafety level"]])
for (i in seq_along(blTaxNames)) {
    blPlot <- plotCladeDis("biosafety level", treeIndex = i)
    blPlot <- facet_widths(blPlot, widths = c(4, 1))
    blPlotName <- file.path(
        "discrete_clade_plots",
        paste0("biosafety_level_", blTaxNames[i], ".png")
    )
    ggsave(filename = blPlotName, plot = blPlot, width = 20, height = 15)
}
# blPlot <- plotCladeDis("biosafety level")
# blPlot <- facet_widths(blPlot, widths = c(4, 1))
# ggsave(
#     filename = "discrete_clade_plots/biosafety_level.png",
#     plot = blPlot, width = 20, height = 15
# )

## gram stain
gsTaxNames <- names(subTrees[["gram stain"]])
for (i in seq_along(gsTaxNames)) {
    gsPlot <- plotCladeDis("gram stain", treeIndex = i)
    gsPlot <- facet_widths(gsPlot, widths = c(4, 1))
    gsPlotName <- file.path(
        "discrete_clade_plots",
        paste0("gram_stain_", gsTaxNames[i], ".png")
    )
    ggsave(filename = gsPlotName, plot = gsPlot, width = 20, height = 15)
}
# gsPlot <- plotCladeDis("gram stain")
# gsPlot <- facet_widths(gsPlot, widths = c(4, 1))
# ggsave(
#     filename = "discrete_clade_plots/gram_stain.png",
#     plot = gsPlot, width = 20, height = 15
# )

## shape
sTaxNames <- names(subTrees[["shape"]])
for (i in seq_along(sTaxNames)) {
    sPlot <- plotCladeDis("shape", treeIndex = i)
    sPlot <- facet_widths(sPlot, widths = c(4, 1))
    sPlotName <- file.path(
        "discrete_clade_plots",
        paste0("shape_", sTaxNames[i], ".png")
    )
    ggsave(filename = sPlotName, plot = sPlot, width = 20, height = 15)
}
# sPlot <- plotCladeDis("shape")
# sPlot <- facet_widths(sPlot, widths = c(4, 1))
# ggsave(
#     filename = "discrete_clade_plots/shape.png",
#     plot = sPlot, width = 20, height = 15
# )

# Roseburia - Animal Pathogen ---------------------------------------------
## These are a few annotations for species in the genus Rosenburia.
## However, these annotations could not be mapped to the LTP tree.
## So these annotations were not used for the ASR.
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
