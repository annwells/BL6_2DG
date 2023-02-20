# Basic Plot Theme
basic_theme <- theme_bw(base_size=14) +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(size=I(0.25), color="black")
  )
