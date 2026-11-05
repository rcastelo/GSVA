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

## memory in bytes of the matrices of the assays among the objects of
## 'heldmem', which forked workers share with the main R process, see
## .held_shared_mem()
gsva_global$heldshared <- 0

## warnings of blocks narrower than the chunks of the input data on disk, see
## .check_chunk_blocks(), collected while gsva() runs its steps, to give one
## warning for all of them at its end, or NULL when they are not collected,
## see .start_chunk_blocks()
gsva_global$chunkblocks <- NULL
