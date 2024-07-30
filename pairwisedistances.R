library(bugphyzz)
library(taxPPro)
library(dplyr)
library(castor)
library(tibble)
library(readr)
library(ggplot2)
library(purrr)
library(tidyr)


# Functions ---------------------------------------------------------------
removeASR <- function(x) {
    x |>
        dplyr::filter(Evidence != "asr") |>
        mutate(NCBI_ID = as.character(NCBI_ID))
}
addAttributes <- function(attrDat, tipDat) {
    dplyr::left_join(
        tipDat, attrDat,
        by = c("taxid" = "NCBI_ID"), relationship = "many-to-many"
    ) |>
        filter(!is.na(Attribute_value))
}

removeDups <- function(x) {
    x |>
        dplyr::rowwise() |>
        dplyr::mutate(
            sorted = paste(sort(c(node1, node2)), collapse = "_")
        ) |>
        dplyr::ungroup() |>
        dplyr::filter(!duplicated(sorted)) |>
        dplyr::select(-sorted)
}
getClosestTips <- function(tr) {
    purrr::map(seq_along(tr$tip.label), function(i) {
        if ((i %% 5000) == 0)
            message(i)
        res <- castor::find_nearest_tips(tr, target_tips = tr$tip.label[i])
        tip_distances <- res[["nearest_distance_per_tip"]]
        names(tip_distances) <- tr$tip.label
        tip_distances <- tip_distances[-i]
        minValue <- min(tip_distances)
        tip_distances <- tip_distances[which(tip_distances == minValue)]
        data.frame(
            node2 = names(tip_distances),
            distance = unname(tip_distances)
        )
    }) |>
        purrr::set_names(tr$tip.label) |>
        dplyr::bind_rows(.id = "node1") |>
        removeDups()
}
# -------------------------------------------------------------------------

ltp <- ltp()
tr <- ltp$tree
tip_data <- ltp$tip_data |>
    group_by(taxid) |>
    mutate(n = n()) |>
    ungroup() |>
    filter(n == 1) |>
    select(-n)
bp <- importBugphyzz()

dat_name <- "spore shape"
dat <- bp[[dat_name]] |>
    removeASR() |>
    select(NCBI_ID, Attribute_value)

tipDat <- left_join(
    select(tip_data, tip_label, taxid),
    dat,
    by = c("taxid" = "NCBI_ID")) |>
    filter(!is.na(Attribute_value))

annotationsL <- split(tipDat, tipDat$Attribute_value) |>
    map(~ pull(.x, tip_label))

closeTips <- getClosestTips(tr)

for (i in seq_along(annotationsL)) {
    annName <- names(annotationsL)[i]

    colName1 <- paste0(annName, "_node1")
    colName2 <- paste0(annName, "_node2")

    closeTips[[colName1]] <- closeTips$node1 %in% annotationsL[[i]]
    closeTips[[colName2]] <- closeTips$node2 %in% annotationsL[[i]]

    closeTips[[annName]] <- paste0(closeTips[[colName1]], "|", closeTips[[colName2]])
}


tbl <- closeTips |>
    select(-matches("_node\\d$")) |>
    pivot_longer(
        names_to = "Attribute", values_to = "logical",
        cols = aerobic:last_col()
    )



attrs <- c(
    "aerobic", "anaerobic", "facultatively anaerobic",
    "elliptic", "endospore", "coccus"
)

tbl |>
    filter(logical != "FALSE|FALSE") |>
    filter(Attribute %in% attrs) |>
    mutate(Attribute = forcats::fct_inorder(Attribute)) |>
    mutate(
        grp = case_when(
            Attribute %in% c("aerobic", "anaerobic", "facultatively anaerobic") ~ "aer",
            TRUE ~ "shape"
        )
    ) |>
    mutate(
        logical = ifelse(logical != "TRUE|TRUE", FALSE, logical)
    ) |>
    # group_by(grp, logical) |>
    # summarise(
    #     mean = mean(distance),
    #     sd = sd(distance)
    # )
    ggplot(aes(x = grp, y = distance)) +
    geom_boxplot(aes(fill = logical))








