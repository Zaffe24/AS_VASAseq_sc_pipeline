#### Prioritize r library of active conda env ####

r_paths=.libPaths()
.libPaths(tail(r_paths,1))

library(Seurat)
library(stringr)
library(ggplot2)

packageVersion("Seurat")

args = commandArgs(trailingOnly=TRUE)

print(args)

if (length(args)!=7) {
    stop("Wrong number of arguments provided!", call.=FALSE)
} else{
    path_cleaned_files=args[1]
    seurat_object=args[2]
    input_string=args[3]
    Thr=as.numeric(args[4])
    out=args[5]
    label=args[6]
    working_dir=args[7]
}

total.Patient <- readRDS(file = seurat_object)
cell_id <- row.names(total.Patient@meta.data)

if (as.character(label) == 'barcode'){
    trans_cell_id<-sapply(cell_id, function(x) strsplit(x,'-v')[[1]][2])
} else
    {trans_cell_id<-gsub("[^0-9]", "", cell_id)}


un<-rep(1,length(cell_id))

print(paste0('labels: ',input_string,'.'))

# converting string of cluster labels into list 
if (input_string!=''){
    input_string<-gsub(" ","", input_string)
    elements <- strsplit(input_string, ";")
    key_value_pairs <- strsplit(elements[[1]], ":")
    keys <- sapply(key_value_pairs, function(x) as.factor(x[1]))
    values <- sapply(key_value_pairs, function(x) x[2])
    result_list <- setNames(values, keys)
    
    #pools <- rep(1,length(paths))
    group<-total.Patient@meta.data$seurat_clusters
    #condition <- ifelse(group=='2','LSC',ifelse(group=='1','Cluster_1','Cluster_0'))
    
    if (sum(sort(as.character(names(result_list))) == sort(as.character(levels(group))))==length(names(result_list))){
        condition<-lapply(group, function(x){result_list[[as.character(x)]]})
        condition<-unlist(condition)
        print('Generation of new cluster labeles completed')} else 
        {stop(paste0('The clusters specified do not match those present in ',seurat_object))}
    
} else{
    condition<-total.Patient@meta.data$seurat_clusters
}


seurat_df<-data.frame(trf=trans_cell_id,ids=cell_id,unit=un,pools=un,condition=condition)
#head(seurat_df)


fracs <- read.csv(file.path(path_cleaned_files,"fractions.txt"), col.names=c('sample','Mouse_frac','Ambig_frac'), header = 0,sep='\t')

if(any(duplicated(fracs$sample)) | any(is.na(fracs$sample))){
    stop('Duplicated or missing files in fractions.txt')
} else {print('fractions.txt loaded successfully')}


if(any(is.na(fracs$Mouse_frac)) | any(is.na(fracs$Ambig_frac))){
    stop('missing values in fractions.txt')} else 
    { print('All fractions recorded!')}

fracs$sample<-trimws(fracs$sample)

fracs$cell_id<-sapply(strsplit(fracs$sample, ".cle"), function(x) x[1])


if (as.character(label) == 'barcode'){
    barcode_table<-read.csv(file.path(working_dir,'bc_celseq2.tsv'),sep='\t',header = F, col.names =c('barcodes','ids'),colClasses = 'character')
    fracs$trans_fracs_id<-sapply(fracs$cell_id, function(x) paste0(gsub("[^0-9]", "",strsplit(x,'_')[[1]][1]), '_',
                                                                  subset(barcode_table, ids == strsplit(x,'_')[[1]][3], select = 'barcodes')[[1]]))
} else {fracs$trans_fracs_id <-paste0(gsub("[^0-9]", "", 
                                           lapply(fracs$cell_id, function(x) unlist(strsplit(x,'_'))[c(1)])),lapply(fracs$cell_id,function(x) unlist(strsplit(x,'_'))[3]))}


fracs$fq1 <- paste0(file.path(path_cleaned_files,fracs$sample),'.fq.gz')


if (length(unlist(strsplit(fracs$cell_id[1],'_'))) != 3){
    stop("cell barcode is not correct for the fitlering set-up \n and adjust section above")
} else {print('barcodes are of the correct type!')}

merged<-merge(fracs,seurat_df, by.x='trans_fracs_id',by.y='trf',all.x=T)
merged$pass<-ifelse(is.na(merged$pools),'discarded','passed')

if (dim(merged)[1] != dim(fracs)[1] | dim(merged)[1]==0){
    stop("Dimension of merged dataframe incorrect")}

########## PLOTTING ##########################
png(file.path(path_cleaned_files,'Mouse_reads.png'),width=1000,height=1000,res=100)
ggplot(merged, aes(x=pass, y=Mouse_frac)) + geom_boxplot(outlier.shape = NA) +
    geom_jitter(shape=16, position=position_jitter(0.2)) + 
    xlab("single cells") + ylab('Fraction of Mouse reads')
dev.off()

png(file.path(path_cleaned_files,'Ambigous_reads.png'),width=1000,height=1000,res=100)
ggplot(merged, aes(x=pass, y=Ambig_frac)) + geom_boxplot(outlier.shape = NA) +
    geom_jitter(shape=16, position=position_jitter(0.2)) + 
    xlab("single cells") + ylab('Fraction of Ambigous reads')
dev.off()

##########################################################


#### DISCARD CELLS WITH AMBIGOUS CONTENT > threshold ###
merged$pass<-ifelse(merged$Ambig_frac > Thr,'discarded',merged$pass)

################## SAVE FULL METADATA TABLE ################
write.table(merged,file.path(path_cleaned_files,'full_metadata.tsv'),quote = F,sep = '\t',row.names = F,col.names = T)


filtered_df <- subset(merged,pass=='passed')

######################### SAVE SAMPLE TSV FOR MAIN PIPELINE ######################
write.table(filtered_df[,c('sample','condition','pools')], quote=F, sep = "\t",row.names = F, col.names = T, file =  file.path(path_cleaned_files,"samples.tsv"))


######################### SAVE UNITS TSV FOR MAIN PIPELINE #####################
write.table(filtered_df[,c('sample','unit','fq1')], quote=F, sep = "\t",row.names = F, col.names = T, file = file.path(path_cleaned_files,"units.tsv"))


######################### SAVE LOCAL TSV FOR MICROEXON DISCOVERY #####################
filtered_df$path<-filtered_df$fq1
filtered_df$cluster<-filtered_df$condition

write.table(filtered_df[,c('path','sample','cluster')], quote=F, sep = "\t",row.names = F, col.names = T, file = file.path(path_cleaned_files,"local.tsv"))

print('Filtering process completed')
print(paste0('Output files stored in ', out))
