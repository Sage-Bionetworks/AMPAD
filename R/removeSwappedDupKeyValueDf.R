removeSwappedDupKeyValueDf <- function(x){
  foo <- igraph::graph_from_data_frame(x,directed=FALSE)
  foo <- igraph::simplify(foo,edge.attr.comb = list("first"))
  foo <- igraph::as_data_frame(foo,what='edges')
  return(foo)
}
