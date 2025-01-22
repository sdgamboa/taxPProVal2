library(castor)
library(dplyr)
library(bugphyzz)

dat <- bpSub$`growth temperature`
subTips <- tipDat[which(tipDat$taxid %in% dat$taxid),]
set.seed(2024)
selectTips <- sample(subTips$tip_label, 100)
tipPairs <- as.data.frame(t(combn(selectTips, 2)))
colnames(tipPairs) <- c("tip1", "tip2")
dists <- map2_dbl(tipPairs$tip1, tipPairs$tip2, ~ {
    get_pairwise_distances(tree = tr, A = .x, B = .y)
})
distsDF <- data.frame(dist = dists)
x <- bind_cols(tipPairs, distsDF)
y <- left_join(x, subTips, by = c("tip1" = "tip_label")) |>
    rename(taxid1 = taxid) |>
    left_join(subTips, by = c("tip2" = "tip_label")) |>
    rename(taxid2 = taxid) |>
    left_join(subDat, by = c("taxid1" = "taxid")) |>
    rename(Attribute1 = Attribute, Attribute_value1 = Attribute_value) |>
    left_join(subDat, by = c("taxid2" = "taxid")) |>
    rename(Attribute2 = Attribute, Attribute_value2 = Attribute_value) |>
    mutate(diff = abs(Attribute_value1 - Attribute_value2))

randPair <- function(y) {
    message("Running ", y)
    dat <- bpSub[[y]]
    subTips <- tipDat[which(tipDat$taxid %in% dat$taxid),]
    if (length(subTips$tip_label) < 100) {
        sizeVar <- length(subTips$tip_label)
    } else {
        sizeVar <- 100
    }
    message("Running ", sizeVar)
    set.seed(2024)
    selectTips <- sample(subTips$tip_label, sizeVar)
    tipPairs <- as.data.frame(t(combn(selectTips, 2)))
    colnames(tipPairs) <- c("tip1", "tip2")
    dists <- map2_dbl(tipPairs$tip1, tipPairs$tip2, ~ {
        get_pairwise_distances(tree = tr, A = .x, B = .y)
    })
    distsDF <- data.frame(dist = dists)
    x <- bind_cols(tipPairs, distsDF)
    left_join(x, subTips, by = c("tip1" = "tip_label")) |>
        rename(taxid1 = taxid) |>
        left_join(subTips, by = c("tip2" = "tip_label")) |>
        rename(taxid2 = taxid) |>
        left_join(dat, by = c("taxid1" = "taxid")) |>
        rename(Attribute1 = Attribute, Attribute_value1 = Attribute_value) |>
        left_join(dat, by = c("taxid2" = "taxid")) |>
        rename(Attribute2 = Attribute, Attribute_value2 = Attribute_value) |>
        mutate(diff = abs(Attribute_value1 - Attribute_value2))
}

numAttrs <- c(
    "growth temperature",
    "coding genes",
    "genome size",
    "length",
    "width",
    "optimal ph"
)

res <- map(numAttrs, ~ randPair(.x))
# hh <- randPair("growth temperature")
hh |>
    ggplot(aes(dist, diff)) +
    geom_point()

# randPair("growth temperature")
res |>
    bind_rows() |>
    ggplot(aes(dist, diff)) +
    geom_point(size = 0.1) +
    facet_wrap(~Attribute1, scales = "free_y") +
    labs(
        x = "Distance", y = "Value difference"
    )

# gz <- randPair("genome size") |>
#     ggplot(aes(dist, diff)) +
#     geom_point()
# cg <- randPair("coding genes") |>
#     ggplot(aes(dist, diff)) +
#     geom_point()
# gt <- randPair("growth temperature") |>
#     ggplot(aes(dist, diff)) +
#     geom_point()
# ln <- randPair("length") |>
#     ggplot(aes(dist, diff)) +
#     geom_point()
# wd <- randPair("width") |>
#     ggplot(aes(dist, diff)) +
#     geom_point()
# op <- randPair("optimal ph") |>
#     ggplot(aes(dist, diff)) +
#     geom_point()
#
#









