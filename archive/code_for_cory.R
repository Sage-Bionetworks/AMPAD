synapseClient::synapseLogin()

#get module definitions
moduleDefinitions <- synapseClient::synTableQuery('select * from syn11932957')@values

#get differential expression meta analysis
foobarObj <- synapseClient::synGet('syn11914811')
load(foobarObj@filePath)


#
