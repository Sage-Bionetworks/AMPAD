synapser::synLogin()
foo<-synapser::synGet('syn11914808')
bar <- data.table::fread(foo$path,data.table=F)

foo <- synapser::synGet('syn11914809')
bar <- data.table::fread(foo$path,data.table=F)
