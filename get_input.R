library(mutationsR)
library(dplyr)

outdir = 'path/to/outdir'
dir.create(outdir)

# signature tables
ref_cosmic = 'path/to/ref/file'
ref = read.table(
  ref_cosmic,
  header = T,
  stringsAsFactors = F,
)
ref$SubType = sapply(ref$Type, get_wt_seq)
ref = select(ref, c(Type, SubType, SBS1, SBS2, SBS5, SBS13, SBS88, SBS89, SBS17a, SBS17b, SBS32, SBS35))


prefix = 'msa' # modify this if necessary
mutation_type='SBS'

refsig_path = paste0(outdir, prefix, '_', mutation_type, '_signatures.csv')
write.table(ref, file = refsig_path, sep = ',', row.names = F)


# input mutation table
input_mutmat = 'path/to/mutation/table'
mutmat = read.table(input_mutmat, header = T, stringsAsFactors = F)
colnames(mutmat)[1] = 'Type'
mutmat_columns = colnames(mutmat)
mutmat$SubType = sapply(mutmat$Type, get_wt_seq)
mutmat = select(mutmat, c(mutmat_columns[1], SubType, mutmat_columns[-1]))

dataset = 'mutmat' # modify this if necessary
context = 96
dir.create(paste0(outdir, dataset))

mutmat_path = paste0(outdir, dataset, '/', 'WGS_', dataset, '.', context, '.csv')
write.table(mutmat, file = mutmat_path, sep = ',', row.names = F)

