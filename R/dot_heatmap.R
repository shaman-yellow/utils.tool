dot_heatmap <- 
  function(
           ## long data
           df
           ){
    p <- ggplot(df, aes(x = sample, y = no.name)) +
      geom_point(aes(size = abs(value), color = value), shape = 16) +
      theme_minimal() +
      guides(size = "none") +
      scale_color_gradient2(low = "#3182BDFF", high = "#A73030FF") +
      theme(text = element_text(family = "Times"),
            axis.text.x = element_text(angle = 90))
    return(p)
  }
tile_heatmap <- 
  function(
           ## long data
           df
           ){
    p <- ggplot(df, aes(x = sample, y = no.name)) +
      geom_tile(aes(fill = value), color = "white", height = 1, width = 1, size = 0.2) +
      theme_minimal() +
      scale_fill_gradient2(low = "#3182BDFF", high = "#A73030FF",
                           limits = c(min(df$value), max(df$value))) +
      labs(x = "Sample", y = "Feature ID", fill = "log2 (Feature level)") +
      theme(text = element_text(family = "Times", face = "bold"),
            axis.text = element_text(face = "plain"),
            axis.text.x = element_blank()
      )
    return(p)
  }
## ---------------------------------------------------------------------- 
add_tree.heatmap <- 
  function(
           ## wide data
           df,
           phc.height = 0.3,
           p
           ){
    phr <- dist(df) %>% 
      hclust() %>% 
      ggtree::ggtree(layout = "rectangular", branch.length = "branch.length") +
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
    ## ------------------------------------- 
    phc <- t(df) %>% 
      dist() %>% 
      hclust() %>% 
      ggtree::ggtree(layout = "rectangular", branch.length = "branch.length") +
      ggtree::layout_dendrogram() +
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
    ## ------------------------------------- 
    com <- p %>% 
      aplot::insert_left(phr, width = 0.3) %>% 
      aplot::insert_top(phc, height = phc.height)
    return(com)
  }
## ---------------------------------------------------------------------- 
add_xgroup.heatmap <- 
  function(
           df,
           p
           ){
    p.xgroup <- ggplot(df, aes(y = "Group", x = sample)) +
      geom_point(aes(color = group), size = 6) +
      ggsci::scale_color_simpsons() +
      labs(x = "", y = "", fill = "Group") +
      theme_minimal() +
      theme(
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            text = element_text(family = "Times", face = "bold"),
            plot.margin = unit(c(0, 0, 0, 0), "cm")
      )
    com <- p %>% 
      aplot::insert_bottom(p.xgroup, height = 0.05)
    return(com)
  }
## ---------------------------------------------------------------------- 
add_xgroup.tile.heatmap <- 
  function(
           df,
           p,
           pal = NA
           ){
    expr.pal <- ifelse(is.na(pal),
                       'ggsci::scale_fill_simpsons()',
                       'scale_fill_manual(values = pal)')
    p.xgroup <- ggplot(df, aes(y = "Group", x = sample)) +
      geom_tile(aes(fill = group), color = "white", height = 1, width = 1, size = 0.2) +
      eval(parse(text = expr.pal)) +
      labs(x = "", y = "", fill = "Group") +
      theme_minimal() +
      theme(
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            text = element_text(family = "Times", face = "bold"),
            plot.margin = unit(c(0, 0, 0, 0), "cm")
      )
    com <- p %>% 
      aplot::insert_bottom(p.xgroup, height = 0.05)
    return(com)
  }
## ---------------------------------------------------------------------- 
add_ygroup.tile.heatmap <-
  function(
           df,
           p,
           pal = NA
           ){
    expr.pal <- ifelse(is.na(pal),
                       'ggsci::scale_fill_npg()',
                       'scale_fill_manual(values = pal)')
    p.ygroup <- ggplot(df, aes(x = "Class", y = no.name)) +
      geom_tile(aes(fill = class), color = "white", height = 1, width = 1, size = 0.2) +
      labs(x = "", y = "", fill = "From") +
      eval(parse(text = expr.pal)) +
      theme_minimal() +
      theme(
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            text = element_text(family = "Times", face = "bold"),
            plot.margin = unit(c(0, 0, 0, 0), "cm")
      )
    com <- p %>% 
      aplot::insert_left(p.ygroup, width = 0.02) 
    return(com)
  }
