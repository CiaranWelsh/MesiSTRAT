library('GEOquery')



WORKING_DIR = 'D:/MesiSTRAT/PublicDataSetSearch'

folders = Sys.glob(file.path(WORKING_DIR, 'GDS*'))

soft_files = Sys.glob(file.path(WORKING_DIR, 'GDS*/*.soft.gz'))
soft_files

d = getGEO(filename = "D:/MesiSTRAT/PublicDataSetSearch/GDS360/GDS360.soft.gz")
boxplot(d)
d
t = Table(d)
boxplot(log2(t[,3:length(colnames(t))]))







GDS360
















