library(vroom)
library(dplyr)
library(bugphyzz)
library(taxPPro)
library(purrr)
library(tidyr)
library(ggplot2)

pwd <- vroom(
    file = "pairwise_distances.tsv", delim = "\t", show_col_types = FALSE
)

dat_name <- "aerophilicity"
bp <- importBugphyzz()
dat <- bp[[dat_name]] |>
    dplyr::filter(Evidence != "asr") |>
    mutate(NCBI_ID = as.character(NCBI_ID)) |>
    select(NCBI_ID, Attribute_value)
datl <- split(dat, dat$Attribute_value)

ltp <- ltp()
tip_data <- ltp$tip_data |>
    select(tip_label, taxid)

annotationsL <- datl |>
    map(~{
        select_ncbi <- .x |>
            count(NCBI_ID) |>
            arrange(-n) |>
            filter(n == 1) |>
            pull(NCBI_ID)
        output <- .x |>
            filter(NCBI_ID %in% select_ncbi)
        output <- left_join(output, tip_data, by = c("NCBI_ID" = "taxid"))
        output <- drop_na(select(output, -NCBI_ID))
        output |>
            pull(tip_label)
    })

annotatedTips <- unique(unlist(annotationsL, use.names = FALSE))


pwd2 <- pwd |>
    filter(node1 %in% annotatedTips & node2 %in% annotatedTips)


for (i in seq_along(annotationsL)) {
    annName <- names(annotationsL)[i]

    colName1 <- paste0(annName, "_node1")
    colName2 <- paste0(annName, "_node2")

    pwd2[[colName1]] <- pwd2$node1 %in% annotationsL[[i]]
    pwd2[[colName2]] <- pwd2$node2 %in% annotationsL[[i]]

    pwd2[[annName]] <- paste0(pwd2[[colName1]], "|", pwd2[[colName2]])
}

# tbl <- left_join(pwd, dat, by = c("node1" = "tip_label")) |>
#     rename(node1Val = Attribute_value) |>
#     left_join(dat, by = c("node2" = "tip_label")) |>
#     rename(node2Val = Attribute_value) |>
#     filter(!is.na(node1Val) & !is.na(node2Val))

newDat <- pwd2 |>
    select(-matches("_node\\d$")) |>
    pivot_longer(
        names_to = "Attribute", values_to = "logical",
        cols = aerobic:last_col()
    ) |>
    filter(logical != "FALSE|FALSE")


# newDat |>
#     ggplot(aes(Attribute, distance, group = Attribute, fill = logical)) +
#     geom_boxplot()


x <- newDat |>
    group_by(Attribute, logical) |>
    summarise(
        mean = mean(distance),
        sd = sd(distance)
    ) |>
        ungroup()

x |>
    ggplot(aes(Attribute, mean, fill = logical)) +
    geom_col(position = "dodge")







