#Lollipop plot attempt 1

library(ggplot2)
library(ggalt)

#Protein info (approximate)
protein_length <- 460
snps <- data.frame(
  position = c(456, 439, 367, 241, 235),
  effect = c(
    "missense",
    "missense", 
    "synonymous",
    "missense",
    "missense"
  )
)
lollipop_height <- 0.1
ggplot(snps, aes(x = position)) +
  # Protein backbone
  annotate(
    "segment",
    x = 1,
    xend = protein_length,
    y = 0,
    yend = 0,
    linewidth = 2,
    color = "grey40"
  ) +
  # Lollipop stems
  geom_segment(
    aes(xend = position, 
        y = 0, yend = lollipop_height,
        color = effect),
    linewidth = 1.2
  ) +
  # Lollipop heads
  geom_point(
    aes(y = lollipop_height, color = effect),
    size = 4
  ) +
  geom_text(
    aes(
      y = lollipop_height + 0.05,
      label = position,
      color = effect
    ),
    angle = 55,
    hjust = 0,
    vjust = -0.25,
    size = 3,
  )+
  ylim(-0.05, 0.4)+
  labs(
    title = "Synonymous SNP distribution in STM2913 (putative permease)",
    x = "Amino Acid Position",
    y = NULL,
    color = "Mutation type"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "top"
  )

#Figure caption
# "Figure X. Lollipop plot showing the distribution of five high-confidence synonymous SNPs within the STM2913 gene encoding a putative permease protein. The horizontal bar represents the approximate full-length protein (~460 amino acids), while lollipops indicate the positions of synonymous substitutions inferred from gene annotations and manual inspection in IGV."
