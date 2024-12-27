## 1. Preparación de los datos.

```{r}

# Cargar las librerías necesarias
library(affy)  # Para trabajar con archivos .CEL y ExpressionSet
library(dplyr) # Para manipulación de datos


filter_microarray <- function(allTargets, seed = 123) {
  set.seed(seed)
  filtered <- subset(allTargets, time != "hour 2")
  filtered$group <- interaction(filtered$infection, filtered$agent)
  selected <- do.call(rbind, lapply(split(filtered, filtered$group), function(group_data) {
    if (nrow(group_data) > 4) {
      group_data[sample(1:nrow(group_data), 4), ]
    } else {
      group_data
    }
  }))
  
  original_indices <- match(selected$sample, allTargets$sample)
  
  rownames(selected) <- paste0(selected$sample, ".", original_indices)
  
  selected$group <- NULL
  return(selected)
}

# Cargar el archivo allTargets
allTargets <- read.table("C:/Users/silvi/Desktop/MSc Bioinformática y Bioestadística/202425-1/Anàlisi de dates òmiques/PACS/PAC2/Muestras CEL/allTargets.txt", header = TRUE, sep = " ")

# Verificar que se ha cargado correctamente
head(allTargets)

# Aplicar la función de filtrado
filteredTargets <- filter_microarray(allTargets, seed=53635628)  # Cambia "123" por tu identificador de la UOC

# Verificar que la función se ha aplicado correctamente
head(filteredTargets)

# Cargar los archivos CEL de las muestras seleccionadas
celFiles <- paste0("C:/Users/silvi/Desktop/MSc Bioinformática y Bioestadística/202425-1/Anàlisi de dates òmiques/PACS/PAC2/Muestras CEL/", filteredTargets$sample, ".CEL")

# Leer los archivos CEL para crear un ExpressionSet
rawData <- ReadAffy(filenames = celFiles)

# Ver el ExpressionSet creado
rawData

# Verificar las intensidades de expresión génica
#exprs(rawData)

#Guardar el objeto
saveRDS(rawData, file = "rawData.rds")
```

```{r}

# Revisar la clase del objeto creado
class(rawData)
```

Para crear el objeto ExpressionSet, se deben tener en cuenta los metadatos de filteredTargets, por lo que:

```{r}

library(Biobase)

# Extraer las intensidades de expresión de rawData
expr_matrix <- exprs(rawData)

# Eliminar la extensión .CEL de los nombres de las columnas.
colnames(expr_matrix) <- sub("\\.CEL$", "", colnames(expr_matrix))

# Como las filas de filteredTargets deben coincidir con las muestras de expr_matrix, se ordenan para que coincidan con las muestras seleccionadas
filteredTargets <- filteredTargets[match(filteredTargets$sample, colnames(expr_matrix)), ]

# Eliminar la información adicional después del . en los nombres de las filas
rownames(filteredTargets) <- sub("\\..*","", rownames(filteredTargets))

# Verificar si los nombres coinciden
if(!all(filteredTargets$sample == colnames(expr_matrix))) {
  stop("Los nombres de las muestras no coinciden entre 'filteredTargets' y 'expr_matrix'")
}

# Crear un data frame con una columna que contiene los nombres de las filas de expr_matrix (las sondas)
feature_data <- data.frame(Gene = rownames(expr_matrix))
```

Se genera el ExpressionSet:

```{r}

expression_set <- ExpressionSet(
  # Añadir en assayData la matrix de expresión génica
  assayData = expr_matrix,  
  # Añadir en phenoData los metadatos
  phenoData = new("AnnotatedDataFrame", data = filteredTargets),  
  # featureData debe contener información relacionada con el nombre de las sondas
  featureData = new("AnnotatedDataFrame", data = feature_data)  
)
```

```{r}

# Verificar la clase del objeto
print(class(expression_set))  # Debería devolver "ExpressionSet"
```

Se guarda el objeto ExpressionSet:

```{r}

saveRDS(expression_set, file = "C:/Users/silvi/Desktop/MSc Bioinformática y Bioestadística/202425-1/Anàlisi de dates òmiques/PACS/PAC2/Muestras CEL/expressionSet.rds")
```

## 2. Análisis exploratorio y de control de calidad.

Ahora, se revisa el objeto ExpressionSet para comprobar que se ha generado correctamente:

```{r}

dim(exprs(expression_set))
```

```{r}

head(expression_set)
```

```{r}

dim(pData(expression_set))
head(pData(expression_set))
```

```{r}

dim(fData(expression_set))
head(fData(expression_set))
```

Se realizan aproximaciones de calidad con el paquete arrayQualityMetrics:

```{r}

library(arrayQualityMetrics)
#arrayQualityMetrics(expression_set)
```

Realizar el análisis de componentes principales:

```{r}

library(ggplot2)
library(ggrepel)

# Generar la función para calcular el PCA
plotPCA3 <- function (datos, labels, factor, title, scale,colores, size = 1.5, glineas = 0.25) {
   data <- prcomp(t(datos),scale=scale)
   dataDf <- data.frame(data$x)
   Group <- factor
   loads <- round(data$sdev^2/sum(data$sdev^2)*100,1)
   
   # Creación del gráfico
   p1 <- ggplot(dataDf,aes(x=PC1, y=PC2)) +
     theme_classic() +
     geom_hline(yintercept = 0, color = "gray70") +
     geom_vline(xintercept = 0, color = "gray70") +
     geom_point(aes(color = Group), alpha = 0.55, size = 3) +
     coord_cartesian(xlim = c(min(data$x[,1])-5,max(data$x[,1])+5)) +
     scale_fill_discrete(name = "Group")
   p1 + geom_text_repel(aes(y = PC2 + 0.25, label = labels),segment.size = 0.25, size = size) + 
     labs(x = c(paste("PC1",loads[1],"%")),y=c(paste("PC2",loads[2],"%"))) +  
     ggtitle(paste("Principal Component Analysis for: ",title,sep=" "))+ 
     theme(plot.title = element_text(hjust = 0.5)) +
     scale_color_manual(values=colores)
   }
```

Se revisa el PCA para los diferentes tratamientos y también para el estado de la infección:

```{r}

plotPCA3(exprs(expression_set), labels = expression_set$sample, factor = expression_set$agent, 
          title="Diferentes tratamientos", scale = FALSE, size = 3, 
          colores = c("blue", "red", "green"))
```

```{r}

plotPCA3(exprs(expression_set), labels = expression_set$sample, factor = expression_set$infection, 
          title="Infección", scale = FALSE, size = 3, 
          colores = c("blue", "red"))
```

Finalmente, se realiza el análisis para una combinación de ambos:

```{r}

expression_set$Group <- paste(expression_set$infection, expression_set$agent, sep = "_")
```

```{r}

plotPCA3(exprs(expression_set), labels = expression_set$sample, factor = expression_set$Group, 
          title="Infección y tratamiento", scale = FALSE, size = 3, 
          colores = c("blue", "red", "green", "yellow", "orange", "purple"))
```

Se realiza un boxplot para visualizar la distribución de intensidades de las matrices:

```{r}

boxplot(expression_set, cex.axis=0.5, las=2,  which="all", 
          col = c("blue", "red", "green", "yellow", "orange", "purple")[as.numeric(as.factor(expression_set$Group))],
          main="Distribución de los valores de intensidad crudos")
```

Generar un análisis de clustering jerárquico:

```{r}

clust.euclid.average <- hclust(dist(t(exprs(expression_set))), method = "average")

# Graficar el dendrograma
plot(clust.euclid.average, labels = expression_set$Group, 
     main = "Clustering jerárquico según grupos con datos en crudo", hang = -1, cex = 0.7)
```

Normalización de los datos:

```{r}

eset_rma <- rma(rawData)
```

Control de los datos normalizados:

```{r}

#arrayQualityMetrics(eset_rma)
```

Revisar el análisis de los componentes principales para los datos normalizados:

```{r}

plotPCA3(exprs(eset_rma), labels = eset_rma$sample, factor = eset_rma$infection, 
         title="Datos normalizados", scale = FALSE, size = 3, 
         colores = c("red", "blue"))
```

Comparar la distribución de los datos normalizados:

```{r}

boxplot(eset_rma, cex.axis=0.5, las=2,  which="all", 
          col = c("blue", "red", "green", "yellow", "orange", "purple")[as.numeric(as.factor(filteredTargets$Group))],
          main="Distribución de los valores de intensidad normalizados")
```

Como los datos han mejorado, se completa el ExpressionSet añadiéndole la información faltante:

```{r}

colnames(eset_rma) <- sub("\\.CEL$", "", colnames(expr_matrix))

filteredTargets$Group <- paste(filteredTargets$infection, filteredTargets$agent, sep = "_")

# Asegúrate de que el orden de filteredTargets coincida con las columnas de eset_rma
filteredTargets <- filteredTargets[match(colnames(exprs(eset_rma)), filteredTargets$sample), ]

# Crear el objeto AnnotatedDataFrame
pheno_data <- new("AnnotatedDataFrame", data = filteredTargets)

# Asignar el phenoData al ExpressionSet
phenoData(eset_rma) <- pheno_data

# Verificar que el phenoData ha sido añadido correctamente
pData(eset_rma)
```

Generar un dendograma a partir de clustering jerárquico para revisar la similitud entre las muestras:

```{r}

clust.euclid.average <- hclust(dist(t(exprs(eset_rma))), method = "average")

plot(clust.euclid.average, labels = eset_rma$Group, 
     main = "Clustering jerárquico según grupo con datos normalizados", hang = -1, cex = 0.7)
```

## 3. Filtrado de datos.

```{r}

# Extraer la matriz de expresión normalizada
expr_matrix_norm <- exprs(eset_rma)

# Calcular la desviación estándar de cada sonda
variabilidad <- apply(expr_matrix_norm, 1, sd)

# Ordenar por su SD en orden descendiente, para visualizar los genes con mayor variabilidad primero
variabilidad_orden <- sort(variabilidad, decreasing = TRUE)

# Seleccionar el 10% de los genes con mayor variabilidad
top10 <- variabilidad_orden[1:floor(length(variabilidad_orden)*0.1)]

# Filtrar la matriz de expresión con los top10 genes
expr_matrix_filter <- expr_matrix_norm [rownames(expr_matrix_norm) %in% names(top10), ]

# Revisar la dimensión y la clase del objeto
dim(expr_matrix_filter)
class(expr_matrix_filter)
```

```{r}

# Crear el objeto ExpressionSet
expression_set_filtered <- ExpressionSet(
  assayData = expr_matrix_filter,  
  phenoData = new("AnnotatedDataFrame", data = filteredTargets))  
```

```{r}
# Realizar un análisis de calidad de los datos normalizados y filtrados
#arrayQualityMetrics(expression_set_filtered)
```

```{r}

# Revisar el análisis de componentes principales de los datos filtrados
plotPCA3(exprs(expression_set_filtered), labels = expression_set_filtered$sample, factor = expression_set_filtered$infection, 
         title="Datos normalizados y filtrados", scale = FALSE, size = 3, 
         colores = c("red", "blue"))
```

Comprobar si la simetria de los datos normalizados y filtrados ha mejorado:

```{r}

boxplot(expression_set_filtered, cex.axis=0.5, las=2,  which="all", 
          col = c("blue", "red", "green", "yellow", "orange", "purple")[as.numeric(as.factor(expression_set_filtered$Group))],
          main="Distribución de los valores de intensidad normalizados y filtrados")
```

Revisar también el dendograma de los dastos normalizados y filtrados:

```{r}

clust.euclid.average <- hclust(dist(t(exprs(expression_set_filtered))), method = "average")

plot(clust.euclid.average, labels = expression_set_filtered$Group, 
     main = "Clustering jerárquico según grupo con datos normalizados y filtrados", hang = -1, cex = 0.7)
```

Guardamos los objetos:

```{r}

write.csv(exprs(eset_rma), file="normalized.Data.csv")
write.csv(exprs(expression_set_filtered), file="normalized.Filtered.Data.csv")
save(eset_rma, expression_set_filtered, file="normalized.Data.Rda")
```

## 4. Construcción de las matrices de diseño y de contraste.

Matriz de diseño:

```{r}

library(limma)

# Para crear la matrix de diseño, indicamos '~0+Group' para que no se incluya un intercepto y cada grupo se compare independientemente
designMat <- model.matrix(~0+Group, pData(eset_rma))

# Asignar nombres para cada condición experimental
colnames(designMat) <- c("S. aureus USA300_linezolid", "uninfected_linezolid", "S. aureus USA300_untreated", "uninfected_untreated", "S. aureus USA300_vancomycin", "uninfected_vancomycin")

print(head(designMat))
```

Para la matriz de contraste:

```{r}

# Generar la matriz de contraste
cont.matrix <- makeContrasts(
  # Entre infectados y no infectados sin tratamiento
  Untreated = uninfected_untreated - S_aur_SAUSA300_untreated,
  # Entre infectados y no infectados tratados con Linezolid
  Linezolid = uninfected_linezolid - S_aur_SAUSA300_linezolid,
  # Entre infectados y no infectados tratados con Vancomycin
  Vancomycin = uninfected_vancomycin - S_aur_SAUSA300_vancomycin,
  # Indicar los niveles de cada condición experimental
  levels = c("uninfected_untreated", "S_aur_SAUSA300_untreated", "uninfected_linezolid", 
             "S_aur_SAUSA300_linezolid", "uninfected_vancomycin", "S_aur_SAUSA300_vancomycin")
)

print(cont.matrix)
```

## 5. Obtención del listado de genes diferencialmente expresados para cada comparación.

```{r}

# Ajustamos el modelo considerando los grupos
fit <- lmFit(eset_rma, designMat)

# Renombramos las columnas por cada grupo
colnames(fit$coefficients) <- c("uninfected_untreated", "S_aur_SAUSA300_untreated", 
"uninfected_linezolid", "S_aur_SAUSA300_linezolid", 
"uninfected_vancomycin", "S_aur_SAUSA300_vancomycin")

# Aplicamos los contrastes
fit.main <- contrasts.fit(fit, cont.matrix)

# Ajustar los resultados con un enfoque bayesiano para mejorar las estimaciones
fit.main <- eBayes(fit.main)

# Revisar la clase del objeto creado
class(fit.main)
```

```{r}

# Generar la tabla con los resultados del análisis diferencial de expresión para el contraste entre los infectados y no infectados sin tratamiento
topTab_Untreated <- topTable(fit.main, coef = "Untreated", number = nrow(fit.main), adjust = "fdr")

head(topTab_Untreated)
```

```{r}

# Generar la tabla para la comparación entre los infectados y no infectados tratados con Linezolid
topTab_Linezolid <- topTable(fit.main, coef = "Linezolid", number = nrow(fit.main), adjust = "fdr")
head(topTab_Linezolid)
```

```{r}

# Generar la tabla para la comparación entre infectados y no infectados tratados con Vancomycin
topTab_Vancomycin <- topTable(fit.main, coef = "Vancomycin", number = nrow(fit.main), adjust = "fdr")
head(topTab_Vancomycin)
```

## 6. Anotación de los genes.

```{r}

# Función para añadir anotaciones a la tabla de resultados de los genes diferencialmente expresados 
annotatedTopTable <- function(topTab, anotPackage)
{  
  topTab <- cbind(PROBEID=rownames(topTab), topTab)
myProbes <- rownames(topTab)
thePackage <- eval(parse(text = anotPackage))
geneAnots <- select(thePackage, myProbes, c("SYMBOL", "ENTREZID", "GENENAME"))
annotatedTopTab <- merge(x=geneAnots, y= topTab, by.x="PROBEID", by.y="PROBEID")
return(annotatedTopTab)
}
```

```{r}

library(mouse4302.db, lib.loc = "C:/Users/silvi/AppData/Local/R/win-library/4.4")
library (annotate)

# Utilizar la función para agregar las anotaciones a las tablas de resultados para los tres tratamientos
topAnnotated_Untreated <- annotatedTopTable(topTab_Untreated,                          anotPackage = "mouse4302.db")
topAnnotatedLinezolid <- annotatedTopTable(topTab_Linezolid, anotPackage = "mouse4302.db")
topAnnotatedVancomycin <- annotatedTopTable(topTab_Vancomycin, anotPackage = "mouse4302.db")

# Guardar los resultados
write.csv(topAnnotated_Untreated, file = "topAnnotatedUntreated.csv")
write.csv(topAnnotatedLinezolid, file = "topAnnotatedLinezolid.csv")
write.csv(topAnnotatedVancomycin, file = "topAnnotatedVancomycin.csv")
```

```{r}

# Revisar las primeras líneas
head(topAnnotated_Untreated)
head(topAnnotatedLinezolid)
head(topAnnotatedVancomycin)
```

## 7. Expresión diferencial.

Generar volcano plots para visualizar la expresión diferencial

```{r}

# Para el grupo sin tratamiento, con coefnum=1
coefnum = 1  
# Calcular los p-values
p_values <- fit.main$p.value[, coefnum]
# Ajustar el tamaño de las etiquetas del gráfico
opt <- par(cex.lab = 0.7)

# Crear el gráfico volcano plot, escogiendo con highlight destacar los 5 puntos más significativos
volcanoplot(fit.main, coef = coefnum, highlight = 5, names = topAnnotated_Untreated$SYMBOL,
            main = paste("Genes expresados diferencialmente", colnames(cont.matrix)[coefnum], sep="\n"))

# Asignar color rojo a los puntos con p-value < 0.05, y negro a los demás
points(fit.main$coefficients[, coefnum], -log10(p_values), 
       col = ifelse(p_values < 0.05, "red", "black"), pch = 16, cex = 0.5)

# Añadir línea de corte en el eje X 
abline(v=c(-1,1))

# Restaurar la configuración de los gráficos
par(opt)
```

```{r}

# Aplicar nuevamente el script para el grupo Linezolid
coefnum = 2  #

p_values <- fit.main$p.value[, coefnum]

opt <- par(cex.lab = 0.7)

volcanoplot(fit.main, coef = coefnum, highlight = 5, names = topAnnotatedLinezolid$SYMBOL,
            main = paste("Genes expresados diferencialmente", colnames(cont.matrix)[coefnum], sep="\n"))

points(fit.main$coefficients[, coefnum], -log10(p_values), 
       col = ifelse(p_values < 0.05, "red", "black"), pch = 16, cex = 0.5)
abline(v=c(-1,1))

par(opt)
```

```{r}

# Generar el volcano plot para el grupo Vancomycin
coefnum = 3  

p_values <- fit.main$p.value[, coefnum]

opt <- par(cex.lab = 0.7)

volcanoplot(fit.main, coef = coefnum, highlight = 5, names = topAnnotatedVancomycin$SYMBOL,
            main = paste("Genes expresados diferencialmente", colnames(cont.matrix)[coefnum], sep="\n"))

points(fit.main$coefficients[, coefnum], -log10(p_values), 
       col = ifelse(p_values < 0.05, "red", "black"), pch = 16, cex = 0.5)
abline(v=c(-1,1))

par(opt)
```

## 8. Comparaciones múltiples.

```{r}

# Realizar la prueba de hipótesis en cada coeficiente por separado
res <- decideTests(fit.main, method = "separate", adjust.method = "fdr", p.value = 0.1, lfc = 1)

# Sumar, por gen, cuántos contrastes fueron significativos para las condiciones seleccionadas (p-value < 0.1y logFC >= 1)
sum.res.rows <- apply(abs(res), 1, sum)

# Filtrar los resultados solo para los genes significativos
res.selected <- res[sum.res.rows!=0,]

# Mostrar un resumen de los resultados
print(summary(res))
```

Para mejor visualización, se genera un diagrama de Venn:

```{r}

vennDiagram(res.selected[,1:3], cex = 0.9)
title("Genes en común entre las tres comparaciones, con FDR < 0.1 y logFC >1 ")
```

## 9. Análisis de la significación biológica.

```{r}

# Crear una lista con la comparación de genes por cada tratamiento
listOfTables <- list(Untreated = topTab_Untreated,
                     Linezolid = topTab_Linezolid,
                     Vancomycin = topTab_Vancomycin)

# Generar una función para calcular cuantos genes son significativos por tratamiento y guardarlos con su EntrezID
listOfSelected <- list()
for (i in 1:length(listOfTables)){
  topTab <- listOfTables[[i]]
  whichGenes <- topTab["adj.P.Val"] < 0.15
  selectIDs <- rownames(topTab)[whichGenes]
  EntrezIDs <- select(mouse4302.db, selectIDs, c("ENTREZID"))
  EntrezIDs <- EntrezIDs$ENTREZID
  listOfSelected[[i]] <- EntrezIDs
  names(listOfSelected)[i] <- names(listOfTables)[i]
}

# Aplicar la función sobre nuestras muestras
sapply(listOfSelected, length)
```

```{r}

# Creamos lista de genes con información sobre las funciones biológicas, procesos y componentes celulares
mapped_genes2GO <- mappedkeys(org.Mm.egGO)

# Creamos lista de genes con información sobre las rutas metabólicas y de señalización celular
mapped_genes2KEGG <- mappedkeys(org.Mm.egPATH)

# Combinamos ambas listas, revisando que no hayan duplicados
mapped_genes <- union(mapped_genes2GO , mapped_genes2KEGG)
```

```{r}
library(ReactomePA)
 
# Crear un script para identificar las vía biológicas relevantes asociadas a los genes diferencialmente expresados. Los resultados se guardane en archivos CSV
listOfData <- listOfSelected[1:2]
comparisonsNames <- names(listOfData)
universe <- mapped_genes
 
for (i in 1:length(listOfData)){
   genesIn <- listOfData[[i]]
   comparison <- comparisonsNames[i]
   enrich.result <- enrichPathway(gene = genesIn,
                                  pvalueCutoff = 0.05,
                                  readable = T,
                                  pAdjustMethod = "BH",
                                  organism = "mouse",
                                  universe = universe)
   
   cat("##################################")
   cat("\nComparison: ", comparison,"\n")
   print(head(enrich.result))
 
   if (length(rownames(enrich.result@result)) != 0) {
   write.csv(as.data.frame(enrich.result), 
              file =paste0("ReactomePA.Results.",comparison,".csv"), 
              row.names = FALSE)
   
   pdf(file=paste0("ReactomePABarplot.",comparison,".pdf"))
     print(barplot(enrich.result, showCategory = 15, font.size = 4, 
             title = paste0("Reactome Pathway Analysis for ", comparison,". Barplot")))
   dev.off()
   
   pdf(file = paste0("ReactomePAcnetplot.",comparison,".pdf"))
     print(cnetplot(enrich.result, categorySize = "geneNum", schowCategory = 15, 
          vertex.label.cex = 0.75))
   dev.off()
   }
 }
```

```{r}

# Generar un gráfico de red de las vías biológicas enriquecidas
cnetplot(enrich.result, categorySize = "geneNum", showCategory = 15, 
          vertex.label.cex = 0.75)
