## Functions

# Discrete ----------------------------------------------------------------
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
    ## This function removes symmetric duplicates.
    ## Maybe is better not to do this since the pairwise distances are not
    ## really symmetrical.
    x |>
        dplyr::rowwise() |>
        dplyr::mutate(
            sorted = paste(sort(c(tip1, tip2)), collapse = "_")
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
            tip2 = names(tip_distances),
            distance = unname(tip_distances)
        )
    }) |>
        purrr::set_names(tr$tip.label)
        # dplyr::bind_rows(.id = "tip1")
        # removeDups()
}
