app <- ShinyDriver$new("../../", loadTimeout = 1e+05)
app$snapshotInit("test1")

app$uploadFile(`matrix1-matrixFile` = "leukemia (1).txt") # <-- This should be the path to the file, relative to the app's tests/shinytest directory
app$uploadFile(`genes1-genesetFile` = "c2.all.v7.0.symbols (1).gmt") # <-- This should be the path to the file, relative to the app's tests/shinytest directory
app$setInputs(button = "click")
app$setInputs(`modal.text-cancel` = "click")
app$setInputs(arg = "yes")
app$setInputs(`argumentsInput-minSz` = 450)
app$setInputs(button = "click")
app$snapshot()
