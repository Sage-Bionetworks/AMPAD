listifyWrapper = function(x,y){
  uniqueKeys <- unique(x)
  res <- lapply(uniqueKeys,
                AMPAD::listify,
                y,
                x)
  names(res) <- uniqueKeys
  return(res)
}
