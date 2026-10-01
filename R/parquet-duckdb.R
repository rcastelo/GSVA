## ----- DuckDB backend to read GSVA Parquet files through HTTP(S) -----
##
## use DuckDB and its extension 'httpfs', which can perform HTTP range requests,
## to download only the footer of the file and the row groups holding the
## requested columns. This requires that the HTTP(S) server supports range
## requests. as with arrow, because the duckdb package is a heavy dependency
## for a functionality that is not always needed, it is a suggested package,
## which is only loaded when reading GSVA output through HTTP(S).

.is_http_url <- function(path) grepl("^https?://", path, ignore.case=TRUE)

#' @importFrom cli cli_abort
.require_duckdb <- function() {
    if (!requireNamespace("duckdb", quietly=TRUE) ||
        !requireNamespace("DBI", quietly=TRUE)) {
        msg <- paste("The 'duckdb' package is required to read GSVA output",
                     "in Apache Parquet format through HTTP(S). Please",
                     "install it with 'install.packages(\"duckdb\")'.")
        cli_abort(c("x"=msg))
    }

    invisible(TRUE)
}

## directory where DuckDB extensions are installed, so that they are kept
## across R sessions, instead of being downloaded again in each one of them
#' @importFrom tools R_user_dir
.duckdb_extension_dir <- function()
    file.path(R_user_dir("GSVA", which="cache"), "duckdb")

## character string as an SQL string literal
.sql_string <- function(x) paste0("'", gsub("'", "''", x, fixed=TRUE), "'")

## DuckDB connections hold pointers to external (C++) objects that cannot be
## serialized to BiocParallel workers, and should not be shared with forked
## processes. For this reason, as with arrow readers, a connection is opened
## on demand and cached by process.
.gsva_duckdb <- new.env(parent=emptyenv())

## connection to an in-memory DuckDB database, with the extension 'httpfs'
## loaded when 'http=TRUE'
#' @importFrom cli cli_abort cli_alert_info
.duckdb_connection <- function(http=FALSE, verbose=FALSE) {
    key <- as.character(Sys.getpid())
    con <- .gsva_duckdb[[key]]
    if (is.null(con) || !DBI::dbIsValid(con)) {
        extdir <- .duckdb_extension_dir()
        dir.create(extdir, recursive=TRUE, showWarnings=FALSE)
        ## the order of the rows read from a file must be the order in which
        ## they are stored, which is DuckDB's default, set here explicitly
        args <- list(config=list(extension_directory=extdir,
                                 preserve_insertion_order="true"))
        ## recent versions of duckdb give messages about its home directory,
        ## used to keep data that GSVA does not need, unless told otherwise
        if ("shared_home" %in% names(formals(duckdb::duckdb)))
            args$shared_home <- FALSE
        con <- DBI::dbConnect(do.call(duckdb::duckdb, args))
        assign(key, con, envir=.gsva_duckdb)
        assign(paste0(key, "|httpfs"), FALSE, envir=.gsva_duckdb)
    }

    if (http) {
        if (!isTRUE(.gsva_duckdb[[paste0(key, "|httpfs")]])) {
            .duckdb_load_httpfs(con, verbose)
            assign(paste0(key, "|httpfs"), TRUE, envir=.gsva_duckdb)
        }
        .duckdb_set_ca_cert_file(con, key)
    }

    con
}

## set in DuckDB the file with certificates of certification authorities to
## verify HTTPS servers given with the option 'GSVA.ca_cert_file', also when
## the option is set, changed or unset after the connection was opened
#' @importFrom cli cli_abort
.duckdb_set_ca_cert_file <- function(con, key) {
    cafile <- getOption("GSVA.ca_cert_file")
    if (!is.null(cafile)) {
        if (!is.character(cafile) || length(cafile) != 1L ||
            !file.exists(cafile))
            cli_abort(c("x"=paste("Cannot find the file {.file {cafile}}",
                                  "given in the option 'GSVA.ca_cert_file'.")))
        cafile <- normalizePath(cafile)
    }

    appliedkey <- paste0(key, "|ca_cert_file")
    if (identical(cafile, .gsva_duckdb[[appliedkey]]))
        return(invisible(TRUE))

    if (is.null(cafile))
        DBI::dbExecute(con, "RESET ca_cert_file")
    else
        DBI::dbExecute(con, paste("SET ca_cert_file =", .sql_string(cafile)))
    assign(appliedkey, cafile, envir=.gsva_duckdb)

    invisible(TRUE)
}

## load the DuckDB extension 'httpfs', which is not part of the duckdb
## package and DuckDB downloads from its own servers the first time it is
## installed
#' @importFrom cli cli_abort cli_alert_info
.duckdb_load_httpfs <- function(con, verbose) {
    installed <- DBI::dbGetQuery(con, paste("SELECT installed FROM",
                                            "duckdb_extensions() WHERE",
                                            "extension_name = 'httpfs'"))
    if (!isTRUE(installed$installed[1])) {
        extdir <- .duckdb_extension_dir()
        if (verbose)
            cli_alert_info(paste("Downloading the DuckDB extension 'httpfs'",
                                 "into {.file {extdir}}"))
        tryCatch(DBI::dbExecute(con, "INSTALL httpfs"), error=function(e) {
            errmsg <- .duckdb_error_message(e)
            cli_abort(c("x"=paste("Cannot install the DuckDB extension",
                                  "'httpfs', necessary to read GSVA output",
                                  "through HTTP(S)."),
                        "i"="{errmsg}",
                        "i"=paste("DuckDB downloads its extensions from",
                                  "'https://extensions.duckdb.org', which",
                                  "requires internet access."),
                        call=NULL))
        })
    }
    DBI::dbExecute(con, "LOAD httpfs")

    invisible(TRUE)
}

## first line of the message of an error given by duckdb, without the
## additional lines with the context of the error
.duckdb_error_message <- function(e) sub("\n.*$", "", conditionMessage(e))

## evaluate a DuckDB query on the file in 'path', giving informative errors
#' @importFrom cli cli_abort
.duckdb_query <- function(path, expr) {
    tryCatch(expr, error=function(e) {
        errmsg <- .hide_credentials(.duckdb_error_message(e), path)
        dpath <- .display_path(path)
        msg <- c("x"="Cannot read {.val {dpath}}.", "i"="{errmsg}")
        if (grepl("SSL|certificate", errmsg, ignore.case=TRUE))
            msg <- c(msg,
                     "i"=paste("The HTTPS server may not be sending the",
                               "intermediate certificates that DuckDB needs",
                               "to verify it, which its administrators can",
                               "fix by configuring the full certificate",
                               "chain. Alternatively, a file with the",
                               "certificates of the certification",
                               "authorities can be given with",
                               "'options(GSVA.ca_cert_file=\"<file>\")'."))
        cli_abort(msg, call=NULL)
    })
}

## key-value metadata and number of rows, (leaf) columns and row groups of the
## Parquet file in 'path', which are cached by process to avoid querying them
## again for the same file
.gsva_duckdb_info <- new.env(parent=emptyenv())

.duckdb_file_info_cached <- function(path)
    !is.null(.gsva_duckdb_info[[paste(Sys.getpid(), path, sep="|")]])

.duckdb_file_info <- function(path, verbose=FALSE) {
    key <- paste(Sys.getpid(), path, sep="|")
    info <- .gsva_duckdb_info[[key]]
    if (!is.null(info))
        return(info)

    con <- .duckdb_connection(http=.is_http_url(path), verbose=verbose)
    src <- .sql_string(path)
    info <- .duckdb_query(path, {
        kv <- DBI::dbGetQuery(con, sprintf(paste(
                  "SELECT decode(key) AS key, decode(value) AS value",
                  "FROM parquet_kv_metadata(%s)",
                  "WHERE starts_with(decode(key), 'gsva_')"), src))
        rg <- DBI::dbGetQuery(con, sprintf(paste(
                  "SELECT count(*) AS n_row_groups, sum(n) AS n_rows FROM",
                  "(SELECT DISTINCT row_group_id, row_group_num_rows AS n",
                  "FROM parquet_metadata(%s))"), src))
        nc <- DBI::dbGetQuery(con, sprintf(paste(
                  "SELECT count(DISTINCT path_in_schema) AS n",
                  "FROM parquet_metadata(%s)"), src))
        md <- as.list(kv$value)
        names(md) <- kv$key
        list(metadata=md,
             num_rows=if (is.na(rg$n_rows)) 0 else as.numeric(rg$n_rows),
             num_columns=as.integer(nc$n),
             num_row_groups=as.integer(rg$n_row_groups))
    })
    assign(key, info, envir=.gsva_duckdb_info)

    info
}

## read the row groups 'gs' (1-based) of the file of the seed 'x', where
## 'what' are the Parquet columns to read, separated by commas. Each run of
## consecutive row groups is read with a query on its range of matrix columns,
## for which DuckDB only reads the row groups that may hold them, using the
## statistics of the column 'col' stored in the footer of the file; results
## are returned in the order of the row groups, as a data frame
.duckdb_read_cols <- function(x, gs, what) {
    con <- .duckdb_connection(http=.is_http_url(x@path))
    lasts <- c(x@firsts[-1] - 1L, x@dim[2])
    runs <- split(gs, cumsum(c(1L, diff(gs) != 1L)))
    res <- lapply(runs, function(rg) {
        q <- sprintf(paste("SELECT %s FROM read_parquet(%s)",
                           "WHERE col BETWEEN %d AND %d"),
                     what, .sql_string(x@path), x@firsts[rg[1]],
                     lasts[rg[length(rg)]])
        .duckdb_query(x@path, DBI::dbGetQuery(con, q))
    })
    if (length(res) == 1L)
        return(res[[1]])
    do.call(rbind, unname(res))
}
