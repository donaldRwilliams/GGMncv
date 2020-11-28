Y <- Sachs
fit <- ggmncv(cor(Y), n = nrow(Sachs[1:3000, ]))
rm(list = ls())


symmetric_mat <- function (x) {
  x[lower.tri(x)] <- t(x)[lower.tri(x)]
  x
}


head.eip <- function(x, n = 5){
  print(x$eip_results[1:n,], row.names = FALSE )
}

head(eip)

get_graph <- function(x){
  returned_object <- list(P = x$P, adj = x$adj)
  class(returned_object) <- "graph"
  return(returned_object)
}

plot.graph <- function(x,
                       layout = "circle",
                       neg_col = "#D55E00",
                       pos_col = "#009E73",
                       edge_magnify = 1,
                       node_size = 10,
                       palette = 2,
                       node_names = NULL,
                       node_groups = NULL){


    x$pcor_adj <- x$P
    p <- ncol(x$P)
    if(is.null(node_names)){
      cn <- 1:p
    } else {
      cn <- node_names
    }

    diag(x$pcor_adj) <- 0

    net <- network::network(x$pcor_adj)

    network::set.edge.value(x = net, attrname = "weights",
                            value = x$pcor_adj)

    network::set.edge.value(x = net, attrname = "abs_weights",
                            value = abs(x$pcor_adj) * edge_magnify)

    network::set.edge.attribute(x = net, attrname = "edge_color",
                                value = ifelse(net %e% "weights" < 0, neg_col,
                                               pos_col))
    e <- abs(as.numeric(x$pcor_adj))

    plt <-GGally::ggnet2(net, edge.alpha = e[e != 0]/max(e),
                         edge.size = "abs_weights", edge.color = "edge_color",
                         node.size = 1, mode = layout)

    if(is.null(node_groups)){
      plt <- plt + geom_point(color = "black",
                 size = node_size + 1) +
      geom_point(size = node_size,
                 color = "white") + guides(color = FALSE) +
      geom_text(label = cn)

  } else {

  plt <-  plt + geom_point(aes(color = node_groups, group = node_groups),
                     size = node_size + 1, alpha = 0.5) +
      geom_point(size = node_size, aes(color = node_groups)) +
      geom_text(label = cn)  +
      scale_color_brewer(palette = palette)


  }
  plt
}

plot(get_graph(fit),
     layout = "circle",
     edge_magnify = 5, node_size = 10,
     node_names = colnames(Sachs),
     node_groups = c(rep("A", 5), rep("B",6)))

