get_plot_style <- function(style = c("paper_small", "paper_large", "poster", "paper_smaller"), ...) {
  style <- match.arg(style)
  params <- list(...)
  
  defaults <- switch(style,
                     paper_small = list(
                       # ---- point aesthetics ----
                       line_w = 0.2,
                       line_c = "darkgray",
                       dot_shape = 16,
                       gdot_size = 0.5,
                       gdot_alpha = 0.7,
                       gdot_col = "gray",
                       dot_size = 0.8,
                       dot_alpha = 1,
                       hdot_size = 1.2,
                       # ---- text & labels ----
                       label_size = 2,
                       
                       # ---- plot lines ----
                       line_interc = 0,
                       axis_line_c = "white",
                       axis_txt = 6,
                       title_txt = 8,
                       legend_txt = 7,
                       font = "Helvetica",
                       
                       # ----- plot dimentions ----
                       plot_height = 7,
                       plot_width = 7.8
                     ),
                     paper_large = list(
                       # ---- point aesthetics ----
                       line_w = 0.3,
                       line_c = "darkgray",
                       dot_shape = 16,
                       gdot_size = 0.7,
                       gdot_alpha = 0.7,
                       gdot_col = "gray",
                       dot_size = 1,
                       dot_alpha = 1,
                       hdot_size = 1.2,
                       # ---- text & labels ----
                       label_size = 2,
                       axis_txt = 6,
                       title_txt = 8,
                       legend_txt = 7,
                       font = "Helvetica",
                       
                       # ---- plot lines ----
                       line_interc = 0,
                       axis_line_c = "white",
                       
                       # ----- plot dimentions ----
                       plot_height = 8,
                       plot_width = 9.5
                     ),
                     poster = list(
                       # ---- point aesthetics ----
                       line_w = 1,
                       line_c = "darkgray",
                       dot_shape = 16,
                       gdot_size = 2.2,
                       gdot_alpha = 0.6,
                       gdot_col = "gray",
                       dot_size = 3,
                       dot_alpha = 0.9,
                       hdot_size = 4,
                       
                       # ---- text & labels ----
                       label_size = 6,
                       axis_txt = 16,
                       title_txt = 24,
                       legend_txt = 18,
                       font = "Helvetica",
                       
                       # ---- plot lines ----
                       line_interc = 0,
                       axis_line_c = "white",
                       
                       # ----- plot dimentions ----
                       plot_height = 13,
                       plot_width = 14
                      ), paper_smaller = list(
                        # ---- point aesthetics ----
                        line_w = 0.2,
                        line_c = "darkgray",
                        dot_shape = 16,
                        gdot_size = 0.5,
                        gdot_alpha = 0.7,
                        gdot_col = "gray",
                        dot_size = 0.8,
                        dot_alpha = 1,
                        hdot_size = 1.2,
                        # ---- text & labels ----
                        label_size = 2,
                        
                        # ---- plot lines ----
                        line_interc = 0,
                        axis_line_c = "white",
                        axis_txt = 6,
                        title_txt = 8,
                        legend_txt = 7,
                        font = "Helvetica",
                        
                        # ----- plot dimentions ----
                        plot_height = 5,
                        plot_width = 5.4
                      )
  )
  
  # merge user overrides
  modifyList(defaults, params)
}

#helper function
available_plot_styles <- function() {
  cat("Available styles: paper_small, paper_large, poster\n",
      "  1: Call theme <- get_plot_style('paper_small')\n",
      "  2: To tweak things a bit:\n",
      "     theme <- get_plot_style('paper_small', font='Arial', line_w=3)\n",
      "  3: Use theme$... throughout the plotting code\n")
}

available_plot_styles
