#inputs for reference are in /seq/vgb/swofford/ref/dndscv
#inputs for samples are in /seq/vgb/swofford/R01/convertingToSKLInput/testVCFs/dndscvInput
#outputs for aggregate are in /seq/vgb/swofford/R01/convertingToSKLInput/testVCFs/dndscvOutput


library(dndscv)

buildref(cdsfile='hg19_mart_export_for_dndscv_orthologFiltered.txt', 
         genomefile='/seq/vgb/swofford/ref/hg19.fa', 
         outfile = "hg19_orthologOnly_refcds.rda", 
         excludechrs="MT", useids = T)


df = read.table('/seq/vgb/swofford/R01/convertingToSKLInput/testVCFs/dndscvInput/hg19/dndscv_input_hg19_dndscv.txt', 
                sep = '\t', header = T)

df$chr <- sub("^", "chr", df$chr)

dndsout = dndscv(df, ref = 'hg19_orthologOnly_refcds.rda')

sel_cv = dndsout$sel_cv

write.table(sel_cv, '/seq/vgb/swofford/temp/hg19_sel_cv_031323.txt', 
            sep = '\t', row.names = F, col.names = T, quote = F)