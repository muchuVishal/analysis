library("RMySQL")
library("reshape2")
source("CodonMap.R")

#***************** General-purpose functions ********************
#****************************************************************

# Connect to database
mydb <- dbConnect(dbDriver("MySQL"), user = "rousniakl", password = "rousn!@k1", dbname = "vest_snvbox", host = "172.16.16.5")

# Get the names of the tables from the db
GetTableList <- function(database) {
	rs <- dbSendQuery(mydb, "show tables;")
	table_name_df <- as.character((fetch(rs, n = -1))[,1])
	dbClearResult(rs)
	table_name_df
}

# Get a df from a table
GetTable <- function(database, tableName) {
	rs <- dbSendQuery(mydb, paste0("
								   select *
								   from ",tableName,";"))
	Table <- fetch(rs, n = -1)
	dbClearResult(rs)
	Table
}

# Select uids
GetCodonPos <- function(database, uid) {
	rs <- dbSendQuery(mydb, paste0("
    select *
    from CodonTable where UID=",uid,";"))
	CodonTable <- fetch(rs, n = -1)
	dbClearResult(rs)
	CodonTable
}

# Get all unique UID's from a chromosome
GetUniqueUID <- function(database,chromosome) {
	rs <- dbSendQuery(mydb, paste0("
								   select distinct UID
								   from CodonTable
								   where chrom=",chromosome,";"))
	UniqueUIDTable <- fetch(rs, n=-1)
	dbClearResult(rs)
	UniqueUIDTable
}

# Split codon base
CodonBaseSplit <- function(df) {
	df$base1 = substr(df$bases,1,1)
	df$base2 = substr(df$bases,2,2)
	df$base3 = substr(df$bases,3,3)
	df$codon = df$bases
	df$bases = NULL
	df
}

# Melt dataframe and merge it
GetGenomicVariants <- function(df) {
	PosCodonTable = df[,-which(names(df) %in% c("base1","base2","base3"))]
	colnames(PosCodonTable)[which(names(PosCodonTable)=="pos1")] = "codpos1"
	colnames(PosCodonTable)[which(names(PosCodonTable)=="pos2")] = "codpos2"
	colnames(PosCodonTable)[which(names(PosCodonTable)=="pos3")] = "codpos3"
	PosCodonTable = melt(PosCodonTable, id.vars = c("UID","chrom","Pos","codon"), value.name = "genpos")

	BaseCodonTable = df[,-which(names(df) %in% c("pos1","pos2","pos3"))]
	colnames(BaseCodonTable)[which(names(BaseCodonTable)=="base1")] = "codpos1"
	colnames(BaseCodonTable)[which(names(BaseCodonTable)=="base2")] = "codpos2"
	colnames(BaseCodonTable)[which(names(BaseCodonTable)=="base3")] = "codpos3"
	BaseCodonTable = melt(BaseCodonTable, id.vars = c("UID","chrom","Pos","codon"), value.name = "base")

	df=merge(PosCodonTable,BaseCodonTable,by = c("UID","chrom","Pos","variable","codon"))
	names(df)[names(df)=="value.x"] = "genpos"
	names(df)[names(df)=="value.y"] = "ref"
	df
}

#***************** Script ***************************************
#****************************************************************

df = CodonTable[,c(4,5,6)]
dfTest = melt(df)
dfTest = dfTest[order(dfTest$value,decreasing=True),]
