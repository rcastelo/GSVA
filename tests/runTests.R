oldWarnOptionValue <- options(warn=1)
BiocGenerics:::testPackage("GSVA")
options(warn=oldWarnOptionValue[["warn"]])
