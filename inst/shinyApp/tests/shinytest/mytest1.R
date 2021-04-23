app <- ShinyDriver$new("../../", loadTimeout = 1e+05)
app$snapshotInit("mytest1")

app$uploadFile(`matrix1-matrixFile` = "matrix.csv")
app$uploadFile(`genes1-genesetFile` = "genes.csv")
app$setInputs(button = "click")
app$setInputs(`modal.text-cancel` = "click")
app$setInputs(Panels = "GeneSets")
app$snapshot()
