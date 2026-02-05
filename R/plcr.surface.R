plcr.surface <- function(formula, covar,indata, init, pl, alpha_lo=-2, alpha_hi=2, beta_lo=-2, beta_hi=2, a_step=0.1, b_step=0.1, title_name=NULL) {

library(ggplot2)
library(ggpubr)
library(plyr)
library(ggExtra)
library(ggfortify)
library(tidyverse)
library(png)
library(ggforce)
library(htmlwidgets)
library(plotly)


for (i in seq(alpha_lo, alpha_hi, by = a_step)) {
  for (j in seq(beta_lo, beta_hi, by = b_step)) {

    single.prof.point <- mixcure.penal.profile.point(formula=formula,
                                                     data = indata, k = covar+1,
                                                     kval = c(i,j), init=init, pl=pl)
    outdat <- rbind.data.frame(outdat, single.prof.point$coefficients[[1]][1:4])

  }
}
colnames(outdat) <- c("var", "kval1","kval2","loglik")

p.exlr <- plot_ly(
  outdat,
  x = ~kval1, y = ~kval2, z = ~as.numeric(loglik),
  type = "scatter3d", mode = "markers",
  marker = list(size = 2, color = ~as.numeric(outdat$loglik), colorscale = "Rainbow")
) |>
  layout(
    title = list(
      text = title_name,
      x = 0.5,                # center title
      xanchor = "center",
      yanchor = "top"
    ),
    scene = list(
      aspectmode = "manual",              # allow custom scaling
      aspectratio = list(x = 1, y = 1, z = 1),
      xaxis = list(title = "logOR"),
      yaxis = list(title = "logHR"),
      zaxis = list(title = "Log-Likelihood")
    ),
    coloraxis = list(colorbar = list(title = "Log-Likelihood")),
    uirevision = TRUE                    # keep user camera on re-renders
  ) |>
  config(displayModeBar = TRUE, scrollZoom = TRUE)  # zoom with wheel; toolbar shown


htmlwidgets::saveWidget(p.exlr, "PLCR_surface.html", selfcontained = TRUE)

}
