#!/usr/bin/env Rscript

################################################
## Resume                                     ##
################################################

# merge tables from annotation and gene quantification
# USAGE: Rscript mergetables.R -a /srv/www/NOBACKUP/BeiRNA/ParaUnir/BeiRNA_JoinedAnnotations_alex.tsv -c /srv/www/NOBACKUP/BeiRNA/ParaUnir/countsfilteredBeiRNA_ControlDMSO_vs_countsfilteredBeiRNA_DEHP.txt -o /srv/www/NOBACKUP/BeiRNA/ParaUnir -n BeiRNA_merge_table
# argument -a: table with the annotations
# argument -c: table with quantifications
# argument -o: path to output
# argument -n: name for output merge table. Example: merge_table_comparacion1

################################################
## LOAD LIBRARIES                             ##
################################################

# Install and/or load  required packages
required_packages <- c("plyr", "dplyr", "tidyr", "stringr", "optparse")
check_package <- lapply(
  required_packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)) {
        install.packages(x, dependencies = TRUE, repos = "http://cran.us.r-project.org")
        library(x, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
    }
  }
)

################################################
################################################
## PARSE COMMAND-LINE PARAMETERS              ##
################################################
################################################

option_list = list(
    make_option(c("-a", "--path_to_annotation_table"), action="store", default=NA, type='character',
              help="Path to annotation table (Trinity table)"),
    make_option(c("-c", "--path_to_count_table"), action="store", default=NA, type='character',
              help="Path to count table (quantification table)"),
    make_option(c("-o", "--output_path"), action="store", default=NA, type='character',
              help="Output path"),
    make_option(c("-n", "--output_name"), action="store", default=NA, type='character',
              help="Output merge table name. Example: merge_table_comparacion1")
              )

opt_parser<- OptionParser(option_list=option_list)
opt<- parse_args(opt_parser)

################################################
## FUNCTIONS     ###############################
################################################


validate_input<- function(file) {
    if (endsWith(file, "csv")) {
        datos<- read.csv2(file, header = T, sep = ",")
        return(datos)
    } else if (endsWith(file, "tsv")) {
        datos<- read.csv2(file, header = T, sep = "\t")
        return(datos)
    } else if (endsWith(file, "txt")) {
        datos<- read.csv2(file, header = T, sep = "\t")
        return(datos)
    } else {
        stop("Must provide file with extension .csv, .tsv or .txt separate by tab")
    }
}

validate_annotation_table<- function(dataframe) {
    if (ncol(dataframe) == 15) {
        return(dataframe[,c("X.ID", "GENENAME", "DESCRIPTION")])
    } else if (ncol(dataframe) == 35) {
        header<- c("ORF.ID", "Gene.name", "Gene.length", "ORF.length", "ORF.start","ORF.end", "Strand", "Protein.sequence","Pfam", "InterPro", "GENENAME", "DESCRIPTION")
        if (sum(colnames(dataframe) %in% header, na.rm = TRUE) == 12) {
            return(dataframe[,c("ORF.ID", "Gene.name", "Gene.length", "ORF.length", "ORF.start","ORF.end", "Strand", "Protein.sequence","Pfam", "InterPro", "GENENAME", "DESCRIPTION")])
        } else {
            stop("Must provide file with some of these columns: ORF ID, Gene name, Gene length, ORF length, ORF start, ORF end, Strand, Protein sequence, Pfam, InterPro, GENENAME and DESCRIPTION")
        }
    } else {
        print("Must provide count table with 15 columns or 35 columns")
    }
}

validate_count_table<- function(dataframe) {
    if (ncol(dataframe) == 6) {
        colnames(dataframe)<- c("ID", "countsfiltered_ControlDMSO_mean", "countsfiltered_DEHP_mean","theta", "prob", "log2FC")
        return(dataframe)
    } else {
        stop("Must provide count table with 6 columns, like this header: ID, countsfilteredBeiRNA_ControlDMSO_mean, countsfilteredBeiRNA_DEHP_mean, theta, prob, log2FC")
    }
}


merge_table<- function(annotation, count, output, name){
    # Cargamos la tabla de anotacion
    datos<- validate_annotation_table(validate_input(annotation))
    if (ncol(datos) == 3) {
        # Separamos los valores "::"
        newid<- str_split(datos$X.ID, "::", simplify = T)[,c(1:3)]
        # Juntamos los 3 primeros characteres en una nueva variable
        datos$ID<- paste0(newid[,1], "::", newid[,2], "::", newid[,3])
        # Cargar la tabla de expresion
        tabla<- validate_count_table(validate_input(count))
        names(tabla)[names(tabla) == "X"] <- "ID"
        # Mergeamos la tabla
        merge_t<- merge(tabla, datos, by = "ID", all.x = TRUE)
        names(merge_t)[names(merge_t) == "X.ID"] <- "old_ID"
        #Creamos la tabla
        write.table(merge_t, paste0(output, "/", name, ".tsv"), sep = "\t", quote = F, row.names = F)
        Sys.chmod(paste0(output, "/", name, ".tsv"), "775")

    } else if (ncol(datos) == 12) {
        # Cargar la tabla de expresion
        tabla<- validate_count_table(validate_input(count))
        # Separamos los valores "~~" de la 
        tabla$newid<- str_split(tabla$ID, "~~", simplify = T)[,c(2)]
        names(tabla)[names(tabla) == "ID"] <- "old_ID"
        names(tabla)[names(tabla) == "newid"] <- "ID"
        #cambiamos ORF.ID por ID en la anotacion
        names(datos)[names(datos) == "ORF.ID"] <- "ID"
        # Mergeamos la tabla
        merge_t<- merge(tabla, datos, by = "ID", all.x = TRUE)
        #Creamos la tabla
        write.table(merge_t, paste0(output, "/", name, ".tsv"), sep = "\t", quote = F, row.names = F)
        Sys.chmod(paste0(output, "/", name, ".tsv"), "775")
    }
}

################################################
## USE           ###############################
################################################

merge_table(opt$path_to_annotation_table, opt$path_to_count_table, opt$output_path, opt$output_name)
