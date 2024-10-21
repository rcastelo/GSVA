## Environment that holds global variables and settings
## for GSVA
gsva_global <- new.env(parent=emptyenv())

## whether start and end messages in the gsva*()
## functions should be shown
gsva_global$show_start_and_end_messages <- TRUE
