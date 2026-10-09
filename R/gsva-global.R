## Environment that holds global variables and settings
## for GSVA
gsva_global <- new.env(parent=emptyenv())

## whether start and end messages in the gsva*()
## functions should be shown
gsva_global$show_start_and_end_messages <- TRUE

## whether the gsva*() functions check the memory required by their
## calculations, which gsva() checks for all of them before running them
gsva_global$check_memory <- TRUE

## memory in bytes that other objects keep allocated in the main R process
## while a step processes its input in blocks, such as the input data of
## gsva() while it runs the steps after the row normalization, see .held_mem()
gsva_global$heldmem <- 0
