## ----first_rse------------------------------------------------
## Lets build our first SummarizedExperiment object
library("SummarizedExperiment") # carga de la libreria
# ?SummarizedExperiment # info del paquete

## De los ejemplos en la ayuda oficial

## Creamos los datos para nuestro objeto de tipo SummarizedExperiment
## para 200 genes a lo largo de 6 muestras
nrows <- 200 # Numero de filas: 200
ncols <- 6 # numero columnas: 6

## Números al azar de cuentas
set.seed(20210223)
counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows) # De 1 a 10000

## Información de nuestros genes
rowRanges <- GRanges( # se crea objeto tipo GRanges
  rep(c("chr1", "chr2"), c(50, 150)),
  IRanges(floor(runif(200, 1e5, 1e6)), width = 100),
  # IRanges: Posiciones a lo largo del cromosoma
  strand = sample(c("+", "-"), 200, TRUE),
  feature_id = sprintf("ID%03d", 1:200) # Id unica para cada gen
)
# Asignar nombres
names(rowRanges) <- paste0("gene_", seq_len(length(rowRanges)))

## Información de nuestras muestras
colData <- DataFrame(
  Treatment = rep(c("ChIP", "Input"), 3),
  row.names = LETTERS[1:6]
)
## Juntamos ahora toda la información en un solo objeto de R
rse <- SummarizedExperiment(
  assays = SimpleList(counts = counts),
  rowRanges = rowRanges,
  colData = colData
)

## Exploremos el objeto resultante
rse

## Número de genes y muestras
dim(rse)

## IDs de nuestros genes y muestras
dimnames(rse)

## Nombres de tablas de cuentas que tenemos (RPKM, CPM, counts, logcounts, etc)
assayNames(rse)

## El inicio de nuestra tabla de cuentas
head(assay(rse))

## Información de los genes en un objeto de Bioconductor
rowRanges(rse)

## Tabla con información de los genes
rowData(rse) # es idéntico a 'mcols(rowRanges(rse))'

## Tabla con información de las muestras
colData(rse)


## ----rse_exercise---------------------------------------------

## Comando 1
rse[1:2, ]
## Repuestas:
## Muestra unicamente dos filas (los genes gene_1 y gene_2)
## y todas sus columnas (A,B,C,D,E,F)
## subconjunto de la matriz

## Comando 2
rse[, c("A", "D", "F")]
## Repuestas:
## Muestra las columas 1,4 y 6, y todas sus filas (gene_1 a gene_200)

## Extra:
identical(rse[,c(1,4,6)],rse[,c("A","D","F")])
## Regresa vector logico que nos indica si los objetos son iguales
# Detiene si la condicón no se cumple:
# stopifnot(identical(rse[,c(1,4,6)],rse[,c("A","D","F")]))

## ----isee_basic, eval = FALSE---------------------------------
# ## Explora el objeto rse de forma interactiva
library("iSEE")
iSEE::iSEE(rse)


## ----download_sce_layer---------------------------------------
## Descarguemos unos datos de spatialLIBD
sce_layer <- spatialLIBD::fetch_data("sce_layer")
sce_layer # filas: 22331 columnas: 76

## Revisemos el tamaño de este objeto
lobstr::obj_size(sce_layer)


## ----explore_sce_layer, eval = FALSE--------------------------
iSEE::iSEE(sce_layer)
