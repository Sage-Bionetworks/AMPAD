trem2_alt_splicing <- function(){
  synapser::synLogin()
  foo <- synapser::synGet('syn12151859')
  bar <- data.table::fread(foo$path)
  bar <- data.frame(bar,stringsAsFactors=F)
  foobar<-dplyr::filter(bar,V1 == 'ENST00000373113.7' | V1 == 'ENST00000338469.3' | V1 == 'ENST00000373122.8')
}
