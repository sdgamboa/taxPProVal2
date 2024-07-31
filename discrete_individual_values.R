## Check "phylogenetic signal" for discrete values

library(bugphyzz)
library(taxPPro)
library(dplyr)
library(castor)
library(tibble)
library(readr)
library(ggplot2)
library(purrr)
library(tidyr)

source("fun.R")

valFname <- system.file(
    "extdata", "validation_summary.tsv", package = "bugphyzz", mustWork = TRUE
)
valData <- read_tsv(valFname, show_col_types = FALSE)

ltp <- ltp()
tr <- ltp$tree

## Drop tips matching more than one taxid
tipData <- ltp$tip_data |>
    group_by(taxid) |>
    mutate(n = n()) |>
    ungroup() |>
    filter(n == 1) |>
    select(-n)

## Get pairwise distances between closest tips only
fname <- "close_tips_distance.tsv"
if (file.exists(fname)) {
    closeTips <- readr::read_tsv(fname, show_col_types = FALSE)
} else {
    closeTips <- getClosestTips(tr)
    readr::write_tsv(x = closeTips, file = fname, num_threads = 10)
}

bp <- importBugphyzz()

## Discrete attributes with ASR annotations

attrMulti <- c(
    ## Multistate
    "spore shape",
    "hemolysis",
    "arrangement",
    "biosafety level",
    "cogem pathogenicity rating",
    "aerophilicity",
    "gram stain",
    "shape"
)

attrBin <- c(
    "extreme environment",
    "host-associated",
    "animal pathogen",
    "antimicrobial sensitivity",
    "motility",
    "plant pathogenity",
    "spore formation"
)

attrNames <- c(attrMulti, attrBin)

physDatL <- bp[attrNames] |>
    map(~ {
        .x |>
            removeASR() |>
            # mutate(
            #     Attribute_value = paste0(Attribute, "|", Attribute_value)
            # ) |>
            select(NCBI_ID, Attribute_value) |>
            as_tibble() |>
            group_by(NCBI_ID) |>
            mutate(
                n = n()
            ) |>
            ungroup() |>
            filter(n == 1) |>
            select(-n)
    })


tipDatL <- physDatL |>
    map(~ {
        left_join(
            x = select(tipData, tip_label, taxid),
            y = .x,
            by = c("taxid" = "NCBI_ID")
        ) |>
            filter(!is.na(Attribute_value))
    })

annotationsL <- tipDatL |>
    map(~ {
        split(.x, .x$Attribute_value)
    }) |>
    map(.f = \(y) map(y, ~ pull(.x, tip_label)))

logicalTbl <- imap(annotationsL, ~ {
    uniqLabs <- unique(unlist(.x, use.names = FALSE))
    ctSub <- closeTips |>
        filter(
            node1 %in% uniqLabs & node2 %in% uniqLabs
        )
    names(.x) <- paste0(.y, "|", names(.x))

    for (i in seq_along(.x)) {
        annName <- names(.x)[i]
        colName1 <- paste0(annName, "_node1")
        colName2 <- paste0(annName, "_node2")

        ctSub[[colName1]] <- ctSub$node1 %in% .x[[i]]
        ctSub[[colName2]] <- ctSub$node2 %in% .x[[i]]

        ctSub[[annName]] <- paste0(ctSub[[colName1]], "|", ctSub[[colName2]])
    }
    return(ctSub)
})


x <- map(logicalTbl, ~ {
    .x |>
        select(matches("_node\\d")) |>
        as.data.frame() |>
        rowSums() |>
        unique()
})




tbl <- map(logicalTbl,  ~ {
    .x |>
        select(-matches("_node\\d$")) |>
        pivot_longer(
            names_to = "Attribute", values_to = "logical",
            cols = 4:last_col()
        ) |>
        filter(logical != "FALSE|FALSE") |>
        arrange(Attribute, logical, distance) |>
        mutate(
            Attribute = forcats::fct_inorder(Attribute)
        ) |>
        mutate(
            grp = sub("\\|.*$", "", Attribute)
        ) |>
        # separate(
        #     col = "Attribute", into = c("grp", "Attribute"), sep = "\\|"
        # ) |>
        mutate(
            logical = ifelse(grepl("FALSE", logical), "FALSE", "TRUE")
        )
}) |>
    bind_rows()

## By Attribute value and by Attribute as a whole
## Number of pairs
## Number of tips # done
## Number of FALSE # done
## Number of TRUE # done
## Mean, median, sd distances # done






uniqNodCounts <- split(tbl, tbl$Attribute) |>
    imap(~ {
        node1 <- pull(.x, node1)
        node2 <- pull(.x, node2)
        val <- length(unique(c(node1, node2)))
        grp <- unique(.x$grp)
        attr <- .y
        tibble(
            grp = grp,
            Attribute = attr,
            nNodes = val
        )
    }) |>
    bind_rows()

ratios <- tbl |>
    count(grp, Attribute, logical) |>
    pivot_wider(
        names_from = "logical", values_from = "n"
    ) |>
    drop_na() |>
    relocate(Attribute, `TRUE`, `FALSE`) |>
    mutate(
        TFratio = `TRUE` / `FALSE`
    ) |>
    select(-`TRUE`, -`FALSE`)

x <- left_join(uniqNodCounts, ratios, by = c("grp", "Attribute"))


x |>
    ggplot(aes(nNodes, TFratio)) +
    geom_point() +
    geom_text(aes(label = Attribute))



sumStats <- tbl |>
    group_by(grp, Attribute, logical) |>
    summarise(
        median = median(distance),
        mean = mean(distance),
        sd = sd(distance),
        n = n()
    )
