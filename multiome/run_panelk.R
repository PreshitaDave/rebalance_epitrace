library(ggplot2)
library(grid)

modeldata <- data.frame(
  stage = c("RG", "Cyc. Prog", "nIPC", "LMO3+", "NR2F1+", "GluN4/5", "GluN2/3"),
  x = c(1, 2, 3, 4, 4, 5, 5),
  y = c(1, 1, 1, 1.3, 0.7, 1.3, 0.7)
)

p_k <- ggplot(modeldata, aes(x = x, y = y, label = stage)) +
  geom_point(size = 6, color = "darkblue") +
  geom_text(vjust = -1, size = 4.5) +
  geom_segment(aes(x = 1, y = 1, xend = 2, yend = 1),
               arrow = arrow(length = unit(0.2, "cm")), linewidth = 1) +
  geom_segment(aes(x = 2, y = 1, xend = 3, yend = 1),
               arrow = arrow(length = unit(0.2, "cm")), linewidth = 1) +
  geom_segment(aes(x = 3, y = 1, xend = 4, yend = 1.3),
               arrow = arrow(length = unit(0.2, "cm")), linewidth = 1) +
  geom_segment(aes(x = 3, y = 1, xend = 4, yend = 0.7),
               arrow = arrow(length = unit(0.2, "cm")), linewidth = 1) +
  geom_segment(aes(x = 4, y = 1.3, xend = 5, yend = 1.3),
               arrow = arrow(length = unit(0.2, "cm")), linewidth = 1) +
  geom_segment(aes(x = 4, y = 0.7, xend = 5, yend = 0.7),
               arrow = arrow(length = unit(0.2, "cm")), linewidth = 1) +
  theme_void() +
  labs(title = "Panel k: EpiTrace model of corticogenesis") +
  theme(plot.title = element_text(hjust = 0.5, size = 16))

ggsave("/projectnb/ds596/projects/Team 3/multiome/panel_k.png", p_k, width = 10, height = 6, dpi = 400)

