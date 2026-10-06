## Environment that holds global variables and settings
## for GSVA
gsva_global <- new.env(parent=emptyenv())

## whether start and end messages in the gsva*()
## functions should be shown
gsva_global$show_start_and_end_messages <- TRUE

## whether the gsva*() functions check the memory required by their
## calculations, which gsva() checks for all of them before running them
gsva_global$check_memory <- TRUE
