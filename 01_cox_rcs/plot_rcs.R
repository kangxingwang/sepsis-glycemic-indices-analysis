# Six individual panels, generated only from aggregate spline/histogram data.
plot_rcs_panels <- function(curves, tests, histograms, output_dir) {
  needed <- c("ggplot2", "patchwork", "ragg")
  absent <- needed[!vapply(needed, requireNamespace, logical(1), quietly=TRUE)]
  if (length(absent)) stop("Plotting packages required: ", paste(absent,collapse=", "))
  dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
  format_p <- function(p) ifelse(p < .001, "<0.001", sprintf("%.3f", p))
  labels <- c(GV="Whole-ICU GV (CV, %)", SHR="SHR", HGI="HGI (percentage points)")
  for (exposure in c("GV","SHR","HGI")) {
    same_exposure <- curves[curves$Exposure == exposure, ]
    limits_y <- c(0, max(same_exposure$Upper)*1.05)
    for (database in c("MIMIC-IV","NWICU")) {
      z <- same_exposure[same_exposure$Database == database, ]
      t <- tests[tests$Database == database & tests$Exposure == exposure, ]
      h <- histograms[histograms$Database == database &
                      histograms$Exposure == exposure, ]
      base_theme <- ggplot2::theme_classic(base_size=10, base_family="Times New Roman") +
        ggplot2::theme(
          axis.line=ggplot2::element_line(linewidth=.4),
          axis.ticks=ggplot2::element_line(linewidth=.35),
          axis.text=ggplot2::element_text(size=9, colour="black"),
          axis.title=ggplot2::element_text(size=10),
          plot.title=ggplot2::element_text(size=11, face="bold"),
          plot.subtitle=ggplot2::element_text(size=9),
          plot.margin=ggplot2::margin(5,5,2,5))
      curve <- ggplot2::ggplot(z, ggplot2::aes(Value, HR)) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin=Lower,ymax=Upper),
                             fill="#326C9A", alpha=.16) +
        ggplot2::geom_hline(yintercept=1, colour="grey50", linetype=2, linewidth=.35) +
        ggplot2::annotate("segment", x=t$Reference, xend=t$Reference,
                           y=0, yend=limits_y[2]*.80, colour="grey60",
                           linetype=3, linewidth=.35) +
        ggplot2::geom_line(colour="#326C9A", linewidth=.65) +
        ggplot2::scale_x_continuous(limits=range(z$Value), expand=c(0,0)) +
        ggplot2::coord_cartesian(ylim=limits_y) +
        ggplot2::labs(
          title=database,
          subtitle=sprintf("n = %s; events = %s",format(t$N,big.mark=","), t$Events),
          x=NULL, y="Adjusted HR (95% CI)") +
        ggplot2::annotate(
          "text", x=-Inf, y=Inf,
          label=paste0("Overall Holm P ",format_p(t$P_overall_Holm),
                        "\nNonlinear Holm P ",format_p(t$P_nonlinear_Holm)),
          hjust=-.04, vjust=1.15, size=3, family="Times New Roman") +
        base_theme +
        ggplot2::theme(axis.text.x=ggplot2::element_blank(),
                        axis.ticks.x=ggplot2::element_blank())
      distribution <- ggplot2::ggplot(h) +
        ggplot2::geom_rect(
          ggplot2::aes(xmin=Left, xmax=Right, ymin=0, ymax=Count),
          fill="#CBD1D6", colour="white", linewidth=.12) +
        ggplot2::scale_x_continuous(limits=range(z$Value), expand=c(0,0)) +
        ggplot2::labs(x=labels[exposure], y="Patients") +
        base_theme +
        ggplot2::theme(
          axis.text.y=ggplot2::element_text(size=8),
          plot.margin=ggplot2::margin(0,5,3,5))
      panel <- patchwork::wrap_plots(curve, distribution, ncol=1, heights=c(4,1))
      path <- file.path(output_dir,
                         paste0("RCS_",gsub("-","",database),"_",exposure,".png"))
      ragg::agg_png(path, width=114, height=100, units="mm", res=600,
                     background="white")
      tryCatch(print(panel), finally=grDevices::dev.off())
    }
  }
  invisible(NULL)
}

