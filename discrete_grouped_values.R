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

selectCols <- c(
    "physiology", "attribute", "mcc_mean",
    "ltp_bp_phys"
)


valData <- valFname |>
    read_tsv(show_col_types = FALSE) |>
    filter(rank == "all") |>
    filter(method == "phytools-ltp") |>
    select(all_of(selectCols)) |>
    group_by(physiology) |>
    mutate(mcc_mean = mean(mcc_mean)) |>
    ungroup() |>
    select(-attribute) |>
    distinct()
    # select(where(~ !all(is.na(.x))))

ltp <- ltp()
tr <- ltp$tree

## Drop tips matching more than one taxid
tipData <- ltp$tip_data |>
    group_by(taxid) |>
    mutate(n = n()) |>
    ungroup() |>
    filter(n == 1)

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
attrMult <- c(
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
    "plant pathogenicity",
    "spore formation"
)
attrNames <- c(attrMult, attrBin)

physDatL <- bp[attrNames] |>
    map(~ {
        .x |>
            removeASR() |>
            select(NCBI_ID, Attribute_value) |>
            as_tibble() |>
            group_by(NCBI_ID) |>
            mutate(
                n = n()
            ) |>
            ungroup() |>
            filter(n == 1) |> # unique annotations per taxa only
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

## Small test to make sure no taxid is repeated per attribute/physiology
tipDatL |>
    map_lgl(~ {
        any(duplicated(pull(.x, taxid)))
    }) |>
    sum()


logicalTbl <- tipDatL |>
    imap(~ {
        tipLabs <- pull(.x, tip_label)
        ctSub <- closeTips |>
            filter(node1 %in% tipLabs & node2 %in% tipLabs)
        ctSub |>
            left_join(
                y = select(.x, tip_label, Attribute_value),
                by = c("node1" = "tip_label")
            ) |>
            rename(attrNode1 = Attribute_value) |>
            left_join(
                y = select(.x, tip_label, Attribute_value),
                by = c("node2" = "tip_label")
            ) |>
            rename(attrNode2 = Attribute_value) |>
            mutate_at(
                .vars = c("attrNode1", "attrNode2"), .funs = as.character
            ) |>
            mutate(
                logical = case_when(
                    attrNode1 == attrNode2 ~ TRUE,
                    TRUE ~ FALSE
                )
            )
    }) |>
    bind_rows(.id = "grp")

logicalTbl |>
    count(grp, logical) |>
    pivot_wider(
        names_from = "logical", values_from = "n",
    ) |>
    relocate(`TRUE`, `FALSE`, .after = grp) |>
    mutate_at(
        .vars = c("TRUE", "FALSE"),
        .funs = \(y) ifelse(is.na(y), 1, y)
    ) |>
    mutate(
        ratio = `TRUE` / `FALSE`
    )

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




ratioTbl <- logicalTbl |>
    count(grp, logical) |>
    pivot_wider(
        names_from = "logical", values_from = "n"
    ) |>
    relocate(`TRUE`, `FALSE`, .after = grp) |>
    mutate_at(
        .vars = c("TRUE", "FALSE"),
        .funs = \(y) ifelse(is.na(y), 1, y)
    ) |>
    mutate(
        ratio = `TRUE` / `FALSE`
    )

logicalTbl |>
    ggplot(aes(grp, distance)) +
    geom_boxplot(aes(color = logical)) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
    )

logicalTbl |>
    group_by(grp, logical) |>
    summarise(
        mean = mean(distance)
    )

logicalTbl |>
    group_by(grp) |>
    summarise(
        mean = mean(distance)
    )







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

sumStats <- tbl |>
    group_by(grp, Attribute, logical) |>
    summarise(
        median = median(distance),
        mean = mean(distance),
        sd = sd(distance),
        n = n()
    )

vd <- valData |>
    mutate(
        Attribute = paste0(physiology, "|", attribute)
    )

y <- left_join(x, vd, by = "Attribute") |>
    filter(rank == "all") |>
    left_join(sumStats, by = "Attribute") |>
    filter(!is.na(mcc_mean)) |>
    select(-n, -median, -sd) |>
    pivot_wider(
        names_from = "logical",  values_from = "mean"
    )

y |>
    ggplot(aes(ltp_bp, mcc_mean)) +
    geom_smooth(method = "lm", formula = y ~ x) +
    geom_hline(
        aes(yintercept = 0.7), linetype = 2, size = 0.1, color = "red"
    ) +
    # geom_point(
    #     aes(size = TFratio, color = `TRUE`)
    # ) +
    geom_point(
        aes(size = TFratio, color = TFratio)
    ) +
    ggrepel::geom_text_repel(aes(label = Attribute)) +
    scale_size(name = "True/False ratio") +
    scale_color_viridis_c(option = "C", name = "Mean distance FALSE") +
    labs(
        x = "Number of tips annotated",
        y = "Mean Matthew's correlation coeffcient"
    ) +
    theme_bw()

