library(ape)
library(ggtree)
library(stringr)

ar_tree <- read.tree("data/gtdbtk.ar122.classify.tree")

ar_tree$tip.label[!str_detect(ar_tree$tip.label, "matabat2bin")] <- ""

ar_tree <- drop.tip(ar_tree, "")
# plot(ar_tree)

test <- data.frame(TipLabel = gsub("matabat2bin.", "", ar_tree$tip.label))
arto <- data.frame(genome = df$genome,
                   order = gsub("p__", "", df$Phylum),
                   taxonomy = paste(df$Family,
                                    df$Genus,
                                    df$Species,
                                    sep = ";"))

test <- test %>%
  left_join(arto, by = c("TipLabel" = "genome"))
test$order[is.na(test$order)] <- "ZZZ"
test$taxonomy[is.na(test$taxonomy)] <- ""

ar_tree$tip.label <- test$taxonomy
ar_tree$tip.label[ar_tree$tip.label == ""] <- "Other"

groupInfo <- split(ar_tree$tip.label, test$order)
names(groupInfo) <- unique(sort(test$order))
ar_tree <- groupOTU(ar_tree, groupInfo)

ggtree1 <- ggtree(ar_tree, aes(colour = group), layout = "circular") +
  geom_tiplab(size = 1, colour = "black", aes(angle = angle)) +
  labs(colour = "Phylum") +
  scale_colour_manual(values = c("grey", "blue", "red", "green", "orange", "dark grey"),
                      labels = c("", sort(unique(test$order))[1:length(unique(test$order)) - 1], "Other"))

ggsave("ar_tree.pdf",
       plot = ggtree1,
       width = 10,
       height = 10,
       device = "pdf",
       path = "results")

## bac

bac_tree <- read.tree("data/gtdbtk.bac120.classify.tree")

bac_tree$tip.label[!str_detect(bac_tree$tip.label, "matabat2bin")] <- ""

bac_tree <- drop.tip(bac_tree, "")
# plot(bac_tree)

test <- data.frame(TipLabel = gsub("matabat2bin.", "", bac_tree$tip.label))
arto <- data.frame(genome = df$genome,
                   order = gsub("p__", "", df$Phylum),
                   taxonomy = paste(df$Family,
                                    df$Genus,
                                    df$Species,
                                    sep = ";"))

test <- test %>%
  left_join(arto, by = c("TipLabel" = "genome"))
test$order[is.na(test$order)] <- "ZZZ"
test$taxonomy[is.na(test$taxonomy)] <- ""

bac_tree$tip.label <- test$taxonomy
bac_tree$tip.label[bac_tree$tip.label == ""] <- "Other"

groupInfo <- split(bac_tree$tip.label, test$order)
names(groupInfo) <- unique(sort(test$order))
bac_tree <- groupOTU(bac_tree, groupInfo)

ggtree2 <- ggtree(bac_tree, aes(colour = group), layout = "circular") +
  geom_tiplab(size = 0.4, colour = "black", aes(angle = angle)) +
  labs(colour = "Phylum") +
  scale_colour_manual(values = c("grey", rainbow(length(unique(test$order)) - 1), "dark grey"),
                      labels = c("", sort(unique(test$order))[1:length(unique(test$order)) - 1], "Other"))

ggsave("bac_tree.pdf",
       plot = ggtree2,
       width = 15,
       height = 15,
       device = "pdf",
       path = "results")

final_tree <- ggtree1 / ggtree2 +
  plot_annotation(tag_levels = 'A')

ggsave("phylotree.pdf",
       plot = final_tree,
       width = 15,
       height = 30,
       device = "pdf",
       path = "results")
