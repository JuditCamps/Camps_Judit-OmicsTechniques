BiocManager::install("Biobase")
library(Biobase)
library(methods)
data(package="Biobase")
data("sample.ExpressionSet")
slotNames(sample.ExpressionSet)
methods(class=class(sample.ExpressionSet))

sample.ExpressionSet
# has 26 rows and 3 columns

dim(exprs(sample.ExpressionSet))
# gives the dimension of the sample

dim(pData(sample.ExpressionSet))
# num of rows and columns

head(pData(sample.ExpressionSet))




samplenames <- letters[1:10]
dataf <- data.frame(treated=sample(c(TRUE, FALSE), 10, replace=TRUE),
                    sex=sample(c("Male", "Female"), 10, replace=TRUE),
                    mood=sample(c("Happy", "Dontt care", "Grumpy"), 10, replace=TRUE),
                    names=samplenames, row.names="names")
dataDesc = data.frame(c("Treated with dark chocolate", "Sex", "Mood while eating"))

pdata <- new("AnnotatedDataFrame", data=dataf,dataDesc)
pdata

pData(pdata)


## -------MIAME---------------------------------
my.desc <- new("MIAME", name="LPS_Experiment",
               lab="National Cancer Institute",
               contact="Lakshman Chelvaraja",
               title="Molecular basis of age associated cytokine dysregulation in LPS stimula"
               url="http://www.jleukbio.org/cgi/content/abstract/79/6/1314")



print(my.desc)

x <- exprs(sample.ExpressionSet)
my.ExpressionSet <- new("ExpressionSet", x, pdata, my.desc)
ExpressionSet(assayData = x, phenoData = pdata)

## Experiment data
##   Experimenter name: LPS_Experiment
##   Laboratory: National Cancer Institute
##   Contact information: Lakshman Chelvaraja
##   Title: Molecular basis of age associated cytokine dysregulation in LPS stimulated mac
##   URL: http://www.jleukbio.org/cgi/content/abstract/79/6/1314
## PMIDs:
##   No abstract available.


data(sample.ExpressionSet)
sample.ExpressionSet
## ExpressionSet (storageMode: lockedEnvironment)
## assayData: 500 features, 26 samples
##   element names: exprs, se.exprs
## protocolData: none
## phenoData
##   sampleNames: A B ... Z (26 total)
##   varLabels: sex type score
##   varMetadata: labelDescription
## featureData: none
## experimentData: use experimentData(object)
## Annotation: hgu95av2

# varmetadata -> explanation if the labels
# samplenames -> rown names

pData()
