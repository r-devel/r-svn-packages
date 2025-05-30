#  File R/Sweave.R
#  Part of the R package Sweave, http://www.R-project.org
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

Sweave <- function(file, driver=RweaveLatex(),
                   syntax=getOption("SweaveSyntax"), ...)
{
    if(is.character(driver))
        driver <- get(driver, mode="function")()
    else if(is.function(driver))
        driver <- driver()


    if(is.null(syntax))
        syntax <- SweaveGetSyntax(file)
    if(is.character(syntax))
        syntax <- get(syntax, mode="list")

    drobj <- driver$setup(file=file, syntax=syntax, ...)
    on.exit(driver$finish(drobj, error=TRUE))

    text <- SweaveReadFile(file, syntax)
    syntax <- attr(text, "syntax")

    mode <- "doc"
    chunknr <- 0
    chunk <- NULL

    namedchunks <- list()
    for(linenum in seq(along=text)){
    	line <- text[linenum]
        if(any(grep(syntax$doc, line))){
            if(mode=="doc"){
                if(!is.null(chunk))
                    drobj <- driver$writedoc(drobj, chunk)
                mode <- "doc"
            }
            else{
                if(!is.null(chunkopts$label))
                    namedchunks[[chunkopts$label]] <- chunk
                if(!is.null(chunk))
                    drobj <- driver$runcode(drobj, chunk, chunkopts)
                mode <- "doc"
            }
            chunk <- NULL
        }
        else if(any(grep(syntax$code, line))){
            if(mode=="doc"){
                if(!is.null(chunk))
                    drobj <- driver$writedoc(drobj, chunk)
                mode <- "code"
            }
            else{
                if(!is.null(chunkopts$label))
                    namedchunks[[chunkopts$label]] <- chunk
                if(!is.null(chunk))
                    drobj <- driver$runcode(drobj, chunk, chunkopts)
                mode <- "code"
            }
            chunkopts <- sub(syntax$code, "\\1", line)
            chunkopts <- SweaveParseOptions(chunkopts,
                                            drobj$options,
                                            driver$checkopts)
            chunk <- NULL
            chunknr <- chunknr+1
            chunkopts$chunknr <- chunknr
        }
        else{
            if(mode=="code" && any(grep(syntax$coderef, line))){
                chunkref <- sub(syntax$coderef, "\\1", line)
                if(!(chunkref %in% names(namedchunks)))
                    warning(gettextf("reference to unknown chunk '%s'",
                                     chunkref), domain = NA)
                line <- namedchunks[[chunkref]]
            }
            srclines <- c(attr(chunk, "srclines"), rep(linenum, length(line)))
            if(is.null(chunk))
                chunk <- line
            else
                chunk <- c(chunk, line)
            attr(chunk, "srclines") <- srclines
        }
    }
    if(!is.null(chunk)){
        if(mode=="doc") drobj <- driver$writedoc(drobj, chunk)
        else drobj <- driver$runcode(drobj, chunk, chunkopts)
    }

    on.exit()
    driver$finish(drobj)
}

SweaveReadFile <- function(file, syntax)
{
    ## file can be a vector to keep track of recursive calls to
    ## SweaveReadFile.  In this case only the first element is
    ## tried to read in, the rest are forbidden names for further
    ## SweaveInput
    f <- file[1]

    bf <- basename(f)
    df <- dirname(f)
    if(!file.exists(f)){
        f <- list.files(df, full.names=TRUE,
                        pattern=paste(bf, syntax$extension, sep=""))

        if(length(f)==0){
            stop(gettextf("no Sweave file with name '%s' found", file[1]),
                 domain = NA)
        }
        else if(length(f) > 1){
            stop(paste(gettextf("%d Sweave files for basename '%s' found:",
                                length(f), file),
                       paste("\n         ", f, collapse="")),
                 domain = NA)
        }
    }

    text <- readLines(f[1])
    ## <FIXME>
    ## This needs to be more refined eventually ...
    if(any(is.na(nchar(text, "c", TRUE)))) {
        ## Ouch, invalid in the current locale.
        ## (Can only happen in a MBCS locale.)
        ## Try re-encoding from Latin1.
        if(capabilities("iconv"))
            text <- iconv(text, "latin1", "")
        else
            stop("Found invalid multi-byte character data.", "\n",
                 "Cannot re-encode because 'iconv' is not available.", "\n",
                 "Try running R in a single-byte locale.")
    }
    ## </FIXME>

    pos <- grep(syntax$syntaxname, text)

    if(length(pos)>1){
        warning(gettextf("more than one syntax specification found, using the first one"), domain = NA)
    }
    if(length(pos)>0){
        sname <- sub(syntax$syntaxname, "\\1", text[pos[1]])
        syntax <- get(sname, mode = "list")
        if(class(syntax) != "SweaveSyntax")
            stop(gettextf("object '%s' does not have class \"SweaveSyntax\"",
                          sname), domain = NA)

        text <- text[-pos]
    }

    if(!is.null(syntax$input)){
        while(any(pos <- grep(syntax$input, text))){

            pos <- pos[1]
            ifile <- file.path(df, sub(syntax$input, "\\1", text[pos]))
            if(any(ifile==file)){
                stop(paste(gettextf("recursive Sweave input '%s' in stack",
                                    ifile),
                           paste("\n         ", 1:length(file), ": ",
                                 rev(file), collapse="")),
                 domain = NA)
            }
            itext <- SweaveReadFile(c(ifile, file), syntax)

            if(pos==1)
                text <- c(itext, text[-pos])
            else if(pos==length(text))
                text <- c(text[-pos], itext)
            else
                text <- c(text[1:(pos-1)], itext, text[(pos+1):length(text)])
        }
    }

    attr(text, "syntax") <- syntax
    text
}



###**********************************************************

SweaveSyntaxNoweb <-
    list(doc = "^@",
         code = "^<<(.*)>>=.*",
         coderef = "^<<(.*)>>.*",
         docopt = "^[[:space:]]*\\\\SweaveOpts\\{([^\\}]*)\\}",
         docexpr = "\\\\Sexpr\\{([^\\}]*)\\}",
         extension = "\\.[rsRS]?nw$",
         syntaxname = "^[[:space:]]*\\\\SweaveSyntax\\{([^\\}]*)\\}",
         input = "^[[:space:]]*\\\\SweaveInput\\{([^\\}]*)\\}",
         trans = list(
             doc = "@",
             code = "<<\\1>>=",
             coderef = "<<\\1>>",
             docopt = "\\\\SweaveOpts{\\1}",
             docexpr = "\\\\Sexpr{\\1}",
             extension = ".Snw",
             syntaxname = "\\\\SweaveSyntax{SweaveSyntaxNoweb}",
             input = "\\\\SweaveInput{\\1}")
         )

class(SweaveSyntaxNoweb) <- "SweaveSyntax"

SweaveSyntaxLatex <- SweaveSyntaxNoweb
SweaveSyntaxLatex$doc <-  "^[[:space:]]*\\\\end\\{Scode\\}"
SweaveSyntaxLatex$code <- "^[[:space:]]*\\\\begin\\{Scode\\}\\{?([^\\}]*)\\}?.*"
SweaveSyntaxLatex$coderef <- "^[[:space:]]*\\\\Scoderef\\{([^\\}]*)\\}.*"
SweaveSyntaxLatex$extension <- "\\.[rsRS]tex$"

SweaveSyntaxLatex$trans$doc <-  "\\\\end{Scode}"
SweaveSyntaxLatex$trans$code <- "\\\\begin{Scode}{\\1}"
SweaveSyntaxLatex$trans$coderef <- "\\\\Scoderef{\\1}"
SweaveSyntaxLatex$trans$syntaxname <- "\\\\SweaveSyntax{SweaveSyntaxLatex}"
SweaveSyntaxLatex$trans$extension <- ".Stex"

###**********************************************************

SweaveGetSyntax <- function(file){

    synt <- apropos("SweaveSyntax", mode="list")
    for(sname in synt){
        s <- get(sname, mode="list")
        if(class(s) != "SweaveSyntax") next
        if(any(grep(s$extension, file))) return(s)
    }
    return(SweaveSyntaxNoweb)
}


SweaveSyntConv <- function(file, syntax, output=NULL)
{
    if(is.character(syntax))
        syntax <- get(syntax)

    if(class(syntax) != "SweaveSyntax")
        stop("target syntax not of class \"SweaveSyntax\"")

    if(is.null(syntax$trans))
        stop("target syntax contains no translation table")

    insynt <- SweaveGetSyntax(file)
    text = readLines(file)
    if(is.null(output))
        output = sub(insynt$extension, syntax$trans$extension, basename(file))

    TN = names(syntax$trans)

    for(n in TN){
        if(n!="extension")
            text = gsub(insynt[[n]], syntax$trans[[n]], text)
    }

    cat(text, file=output, sep="\n")
    cat("Wrote file", output, "\n")
}




###**********************************************************

SweaveParseOptions <- function(text, defaults=list(), check=NULL)
{
    x <- sub("^[[:space:]]*(.*)", "\\1", text)
    x <- sub("(.*[^[:space:]])[[:space:]]*$", "\\1", x)
    x <- unlist(strsplit(x, "[[:space:]]*,[[:space:]]*"))
    x <- strsplit(x, "[[:space:]]*=[[:space:]]*")

    ## only the first option may have no name: the chunk label
    if(length(x)>0){
        if(length(x[[1]])==1){
            x[[1]] <- c("label", x[[1]])
        }
    }
    else
        return(defaults)

    if(any(sapply(x, length)!=2))
        stop(gettextf("parse error or empty option in\n%s", text), domain = NA)

    options <- defaults

    for(k in 1:length(x))
        options[[ x[[k]][1] ]] <- x[[k]][2]

    if(!is.null(options[["label"]]) && !is.null(options[["engine"]]))
        options[["label"]] <- sub(paste("\\.", options[["engine"]], "$",
                                        sep=""),
                                  "", options[["label"]])

    if(!is.null(check))
        options <- check(options)

    options
}

SweaveHooks <- function(options, run=FALSE, envir=.GlobalEnv)
{
    if(is.null(SweaveHooks <- getOption("SweaveHooks")))
        return(NULL)

    z <- character(0)
    for(k in names(SweaveHooks)){
        if(k != "" && !is.null(options[[k]]) && options[[k]]){
            if(is.function(SweaveHooks[[k]])){
                z <- c(z, k)
                if(run)
                    eval(SweaveHooks[[k]](), envir=envir)
            }
        }
    }
    z
}





###**********************************************************


RweaveLatex <- function()
{
    list(setup = RweaveLatexSetup,
         runcode = RweaveLatexRuncode,
         writedoc = RweaveLatexWritedoc,
         finish = RweaveLatexFinish,
         checkopts = RweaveLatexOptions)
}

RweaveLatexSetup <-
    function(file, syntax,
             output=NULL, quiet=FALSE, debug=FALSE, echo=TRUE,
             eval=TRUE, keep.source=FALSE, split=FALSE, pdf=TRUE, eps=FALSE,
             envir=.GlobalEnv)
{
    if(is.null(output)){
        prefix.string <- basename(sub(syntax$extension, "", file))
        output <- paste(prefix.string, "tex", sep=".")
    }
    else{
        prefix.string <- basename(sub("\\.tex$", "", output))
    }
    if(!quiet) cat("Writing to file ", output, "\n",
                   "Processing code chunks ...\n", sep="")
    output <- file(output, open="w+")

    if(is.character(envir)){
        if(envir=="new") envir <- new.env()
    }

    options <- list(prefix=TRUE, prefix.string=prefix.string,
                    engine="R", print=FALSE, eval=eval,
                    fig=FALSE, pdf=pdf, eps=eps,
                    width=6, height=6, term=TRUE,
                    echo=echo, keep.source=keep.source, results="verbatim",
                    split=split, strip.white="true", include=TRUE,
                    pdf.version="1.1", pdf.encoding="default",
                    concordance=FALSE, expand=TRUE,
                    envir=envir)

    ## to be on the safe side: see if defaults pass the check
    options <- RweaveLatexOptions(options)

    list(output=output, haveconcordance=FALSE,
         debug=debug, quiet=quiet, syntax = syntax,
         options=options, chunkout=list(), srclines=integer(0),
         srcfile=srcfile(file))
}

makeRweaveLatexCodeRunner <- function(evalFunc=RweaveEvalWithOpt)
{
    ## Return a function suitable as the 'runcode' element
    ## of an Sweave driver.  evalFunc will be used for the
    ## actual evaluation of chunk code.
    RweaveLatexRuncode <- function(object, chunk, options)
      {
          if(!(options$engine %in% c("R", "S"))){
              return(object)
          }

          if(!object$quiet){
              cat(formatC(options$chunknr, width=2), ":")
              if(options$echo) cat(" echo")
              if(options$keep.source) cat(" keep.source")
              if(options$eval){
                  if(options$print) cat(" print")
                  if(options$term) cat(" term")
                  cat("", options$results)
                  if(options$fig){
                      if(options$eps) cat(" eps")
                      if(options$pdf) cat(" pdf")
                  }
              }
              if(!is.null(options$label))
                cat(" (label=", options$label, ")", sep="")
              cat("\n")
          }

          chunkprefix <- RweaveChunkPrefix(options)

          if(options$split){
              ## [x][[1]] avoids partial matching of x
              chunkout <- object$chunkout[chunkprefix][[1]]
              if(is.null(chunkout)){
                  chunkout <- file(paste(chunkprefix, "tex", sep="."), "w")
                  if(!is.null(options$label))
                    object$chunkout[[chunkprefix]] <- chunkout
              }
          }
          else
            chunkout <- object$output

	  saveopts <- options(keep.source=options$keep.source)
	  on.exit(options(saveopts))

          SweaveHooks(options, run=TRUE)

          chunkexps <- try(parse(text=chunk), silent=TRUE)
          RweaveTryStop(chunkexps, options)
          openSinput <- FALSE
          openSchunk <- FALSE

          if(length(chunkexps)==0)
            return(object)

          srclines <- attr(chunk, "srclines")
          linesout <- integer(0)
          srcline <- srclines[1]

	  srcrefs <- attr(chunkexps, "srcref")
	  if (options$expand)
	    lastshown <- 0
	  else
	    lastshown <- srcline - 1
	  thisline <- 0
          for(nce in 1:length(chunkexps))
            {
                ce <- chunkexps[[nce]]
                if (nce <= length(srcrefs) && !is.null(srcref <- srcrefs[[nce]])) {
                    if (options$expand) {
                	srcfile <- attr(srcref, "srcfile")
                	showfrom <- srcref[1]
                	showto <- srcref[3]
                    } else {
                    	srcfile <- object$srcfile
                    	showfrom <- srclines[srcref[1]]
                    	showto <- srclines[srcref[3]]
                    }
                    dce <- getSrcLines(srcfile, lastshown+1, showto)
	    	    leading <- showfrom-lastshown
	    	    lastshown <- showto
                    srcline <- srclines[srcref[3]]
                    while (length(dce) && length(grep("^[ \\t]*$", dce[1]))) {
	    		dce <- dce[-1]
	    		leading <- leading - 1
	    	    }
	    	} else {
                    dce <- deparse(ce, width.cutoff=0.75*getOption("width"))
                    leading <- 1
                }
                if(object$debug)
                  cat("\nRnw> ", paste(dce, collapse="\n+  "),"\n")
                if(options$echo && length(dce)){
                    if(!openSinput){
                        if(!openSchunk){
                            cat("\\begin{Schunk}\n",
                                file=chunkout, append=TRUE)
                            linesout[thisline + 1] <- srcline
                            thisline <- thisline + 1
                            openSchunk <- TRUE
                        }
                        cat("\\begin{Sinput}",
                            file=chunkout, append=TRUE)
                        openSinput <- TRUE
                    }
		    cat("\n", paste(getOption("prompt"), dce[1:leading], sep="", collapse="\n"),
		    	file=chunkout, append=TRUE, sep="")
                    if (length(dce) > leading)
                    	cat("\n", paste(getOption("continue"), dce[-(1:leading)], sep="", collapse="\n"),
                    	    file=chunkout, append=TRUE, sep="")
		    linesout[thisline + 1:length(dce)] <- srcline
		    thisline <- thisline + length(dce)
                }

                                        # tmpcon <- textConnection("output", "w")
                                        # avoid the limitations (and overhead) of output text connections
                tmpcon <- file()
                sink(file=tmpcon)
                err <- NULL
                if(options$eval) err <- evalFunc(ce, options)
                cat("\n") # make sure final line is complete
                sink()
                output <- readLines(tmpcon)
                close(tmpcon)
                ## delete empty output
                if(length(output)==1 & output[1]=="") output <- NULL

                RweaveTryStop(err, options)

                if(object$debug)
                  cat(paste(output, collapse="\n"))

                if(length(output)>0 & (options$results != "hide")){

                    if(openSinput){
                        cat("\n\\end{Sinput}\n", file=chunkout, append=TRUE)
                        linesout[thisline + 1:2] <- srcline
                        thisline <- thisline + 2
                        openSinput <- FALSE
                    }
                    if(options$results=="verbatim"){
                        if(!openSchunk){
                            cat("\\begin{Schunk}\n",
                                file=chunkout, append=TRUE)
                            linesout[thisline + 1] <- srcline
                            thisline <- thisline + 1
                            openSchunk <- TRUE
                        }
                        cat("\\begin{Soutput}\n",
                            file=chunkout, append=TRUE)
                        linesout[thisline + 1] <- srcline
                        thisline <- thisline + 1
                    }

                    output <- paste(output,collapse="\n")
                    if(options$strip.white %in% c("all", "true")){
                        output <- sub("^[[:space:]]*\n", "", output)
                        output <- sub("\n[[:space:]]*$", "", output)
                        if(options$strip.white=="all")
                          output <- sub("\n[[:space:]]*\n", "\n", output)
                    }
                    cat(output, file=chunkout, append=TRUE)
                    count <- sum(strsplit(output, NULL)[[1]] == "\n")
                    if (count > 0) {
                    	linesout[thisline + 1:count] <- srcline
                    	thisline <- thisline + count
                    }

                    remove(output)

                    if(options$results=="verbatim"){
                        cat("\n\\end{Soutput}\n", file=chunkout, append=TRUE)
                        linesout[thisline + 1:2] <- srcline
                        thisline <- thisline + 2
                    }
                }
            }

          if(openSinput){
              cat("\n\\end{Sinput}\n", file=chunkout, append=TRUE)
              linesout[thisline + 1:2] <- srcline
              thisline <- thisline + 2
          }

          if(openSchunk){
              cat("\\end{Schunk}\n", file=chunkout, append=TRUE)
              linesout[thisline + 1] <- srcline
              thisline <- thisline + 1
          }

          if(is.null(options$label) & options$split)
            close(chunkout)

          if(options$split & options$include){
              cat("\\input{", chunkprefix, "}\n", sep="",
                file=object$output, append=TRUE)
              linesout[thisline + 1] <- srcline
              thisline <- thisline + 1
          }

          if(options$fig && options$eval){
              if(options$eps){
                  grDevices::postscript(file=paste(chunkprefix, "eps", sep="."),
                                        width=options$width, height=options$height,
                                        paper="special", horizontal=FALSE)

                  err <- try({SweaveHooks(options, run=TRUE)
                              eval(chunkexps, envir=options$envir)})
                  grDevices::dev.off()
                  if(inherits(err, "try-error")) stop(err)
              }
              if(options$pdf){
                  grDevices::pdf(file=paste(chunkprefix, "pdf", sep="."),
                                 width=options$width, height=options$height,
                                 version=options$pdf.version,
                                 encoding=options$pdf.encoding)

                  err <- try({SweaveHooks(options, run=TRUE)
                              eval(chunkexps, envir=options$envir)})
                  grDevices::dev.off()
                  if(inherits(err, "try-error")) stop(err)
              }
              if(options$include) {
                  cat("\\includegraphics{", chunkprefix, "}\n", sep="",
                      file=object$output, append=TRUE)
                  linesout[thisline + 1] <- srcline
                  thisline <- thisline + 1
              }
          }
          object$linesout <- c(object$linesout, linesout)
          return(object)
      }
    RweaveLatexRuncode
}

RweaveLatexRuncode <- makeRweaveLatexCodeRunner()

RweaveLatexWritedoc <- function(object, chunk)
{
    linesout <- attr(chunk, "srclines")
    
    while(any(pos <- grep(object$syntax$docexpr, chunk)))
    {
        cmdloc <- regexpr(object$syntax$docexpr, chunk[pos[1]])
        cmd <- substr(chunk[pos[1]], cmdloc,
                      cmdloc+attr(cmdloc, "match.length")-1)
        cmd <- sub(object$syntax$docexpr, "\\1", cmd)
        if(object$options$eval){
            val <- as.character(eval(parse(text=cmd), envir=object$options$envir))
            ## protect against character(0), because sub() will fail
            if(length(val)==0) val <- ""
        }
        else
            val <- paste("\\\\verb{<<", cmd, ">>{", sep="")

        chunk[pos[1]] <- sub(object$syntax$docexpr, val, chunk[pos[1]])
    }
    while(any(pos <- grep(object$syntax$docopt, chunk)))
    {
        opts <- sub(paste(".*", object$syntax$docopt, ".*", sep=""),
                    "\\1", chunk[pos[1]])
        object$options <- SweaveParseOptions(opts, object$options,
                                             RweaveLatexOptions)
        if (isTRUE(object$options$concordance) 
              && !object$haveconcordance) {
            savelabel <- object$options$label
            object$options$label <- "concordance"
            prefix <- RweaveChunkPrefix(object$options)
            object$options$label <- savelabel
            object$concordfile <- paste(prefix, "tex", sep=".")
            chunk[pos[1]] <- sub(object$syntax$docopt, 
                                 paste("\\\\input{", prefix, "}", sep=""),
                                 chunk[pos[1]])
            object$haveconcordance <- TRUE
        } else
            chunk[pos[1]] <- sub(object$syntax$docopt, "", chunk[pos[1]])
    }
    
    cat(chunk, sep="\n", file=object$output, append=TRUE)
    object$linesout <- c(object$linesout, linesout)
    
    return(object)
}

RweaveLatexFinish <- function(object, error=FALSE)
{
    outputname <- summary(object$output)$description
    inputname <- object$srcfile$filename
    if(!object$quiet && !error){
        if(object$options$eps){
            cat("\n",
                gettextf("You can now run LaTeX on '%s'", outputname),
                "\n", sep = "")
        }
        else{
            cat("\n",
                gettextf("You can now run PDFLaTeX on '%s'", outputname),
                "\n", sep = "")
        }
    }
    close(object$output)
    if(length(object$chunkout) > 0)
        for(con in object$chunkout) close(con)
    if (object$haveconcordance) {
    	# This output format is subject to change.  Currently it contains
    	# three parts, separated by colons:
    	# 1.  The output .tex filename
    	# 2.  The input .Rnw filename
    	# 3.  The input line numbers corresponding to each output line.
    	#     This are compressed using the following simple scheme:
    	#     The first line number, followed by
    	#     a run-length encoded diff of the rest of the line numbers.
        linesout <- object$linesout
        vals <- rle(diff(linesout))
        vals <- c(linesout[1], as.numeric(rbind(vals$lengths, vals$values)))
    	concordance <- paste(strwrap(paste(vals, collapse=" ")), collapse=" %\n")
    	special <- paste("\\special{concordance:", outputname, ":", inputname, ":%\n",
    			 concordance,"}\n", sep="")
    	cat(special, file=object$concordfile)
    }
    invisible(outputname)
}

RweaveLatexOptions <- function(options)
{

    ## ATTENTION: Changes in this function have to be reflected in the
    ## defaults in the init function!

    ## convert a character string to logical
    c2l <- function(x){
        if(is.null(x)) return(FALSE)
        else return(as.logical(toupper(as.character(x))))
    }

    NUMOPTS <- c("width", "height")
    NOLOGOPTS <- c(NUMOPTS, "results", "prefix.string",
                   "engine", "label", "strip.white",
                   "pdf.version", "pdf.encoding", "envir")

    for(opt in names(options)){
        if(! (opt %in% NOLOGOPTS)){
            oldval <- options[[opt]]
            if(!is.logical(options[[opt]])){
                options[[opt]] <- c2l(options[[opt]])
            }
            if(is.na(options[[opt]]))
                stop(gettextf("invalid value for '%s' : %s", opt, oldval),
                     domain = NA)
        }
        else if(opt %in% NUMOPTS){
            options[[opt]] <- as.numeric(options[[opt]])
        }
    }

    if(!is.null(options$results))
        options$results <- tolower(as.character(options$results))
    options$results <- match.arg(options$results,
                                 c("verbatim", "tex", "hide"))

    if(!is.null(options$strip.white))
        options$strip.white <- tolower(as.character(options$strip.white))
    options$strip.white <- match.arg(options$strip.white,
                                     c("true", "false", "all"))

    options
}


RweaveChunkPrefix <- function(options)
{
    if(!is.null(options$label)){
        if(options$prefix)
            chunkprefix <- paste(options$prefix.string, "-",
                                 options$label, sep="")
        else
            chunkprefix <- options$label
    }
    else
        chunkprefix <- paste(options$prefix.string, "-",
                             formatC(options$chunknr, flag="0", width=3),
                             sep="")

    return(chunkprefix)
}

RweaveEvalWithOpt <- function (expr, options){
    if(options$eval){
        res <- try(.Internal(eval.with.vis(expr, options$envir, baseenv())),
                   silent=TRUE)
        if(inherits(res, "try-error")) return(res)
        if(options$print | (options$term & res$visible))
            print(res$value)
    }
    return(res)
}


RweaveTryStop <- function(err, options){

    if(inherits(err, "try-error")){
        cat("\n")
        msg <- paste(" chunk", options$chunknr)
        if(!is.null(options$label))
            msg <- paste(msg, " (label=", options$label, ")", sep="")
        msg <- paste(msg, "\n")
        stop(msg, err, call.=FALSE)
    }
}





###**********************************************************

Stangle <- function(file, driver=Rtangle(),
                    syntax=getOption("SweaveSyntax"), ...)
{
    Sweave(file=file, driver=driver, ...)
}

Rtangle <-  function()
{
    list(setup = RtangleSetup,
         runcode = RtangleRuncode,
         writedoc = RtangleWritedoc,
         finish = RtangleFinish,
         checkopts = RweaveLatexOptions)
}


RtangleSetup <- function(file, syntax,
                         output=NULL, annotate=TRUE, split=FALSE,
                         prefix=TRUE, quiet=FALSE)
{
    if(is.null(output)){
        prefix.string <- basename(sub(syntax$extension, "", file))
        output <- paste(prefix.string, "R", sep=".")
    }
    else{
        prefix.string <- basename(sub("\\.[rsRS]$", "", output))
    }

    if(!split){
        if(!quiet)
            cat("Writing to file", output, "\n")
        output <- file(output, open="w")
    }
    else{
        if(!quiet)
            cat("Writing chunks to files ...\n")
        output <- NULL
    }

    options <- list(split=split, prefix=prefix,
                    prefix.string=prefix.string,
                    engine="R", eval=TRUE)

    list(output=output, annotate=annotate, options=options,
         chunkout=list(), quiet=quiet, syntax=syntax)
}


RtangleRuncode <-  function(object, chunk, options)
{
    if(!(options$engine %in% c("R", "S"))){
        return(object)
    }

    chunkprefix <- RweaveChunkPrefix(options)

    if(options$split){
        outfile <- paste(chunkprefix, options$engine, sep=".")
        if(!object$quiet)
            cat(options$chunknr, ":", outfile,"\n")
        ## [x][[1]] avoids partial matching of x
        chunkout <- object$chunkout[chunkprefix][[1]]
        if(is.null(chunkout)){
            chunkout <- file(outfile, "w")
            if(!is.null(options$label))
                object$chunkout[[chunkprefix]] <- chunkout
        }
    }
    else
        chunkout <- object$output

    if(object$annotate){
        cat("###################################################\n",
            "### chunk number ", options$chunknr,
            ": ", options$label,
            ifelse(options$eval, "", " eval=FALSE"), "\n",
            "###################################################\n",
            file=chunkout, append=TRUE, sep="")
    }

    hooks <- SweaveHooks(options, run=FALSE)
    for(k in hooks)
        cat("getOption(\"SweaveHooks\")[[\"", k, "\"]]()\n",
            file=chunkout, append=TRUE, sep="")

    if(!options$eval)
        chunk <- paste("##", chunk)

    cat(chunk,"\n", file=chunkout, append=TRUE, sep="\n")

    if(is.null(options$label) & options$split)
        close(chunkout)

    return(object)
}

RtangleWritedoc <- function(object, chunk)
{
    while(any(pos <- grep(object$syntax$docopt, chunk)))
    {
        opts <- sub(paste(".*", object$syntax$docopt, ".*", sep=""),
                    "\\1", chunk[pos[1]])
        object$options <- SweaveParseOptions(opts, object$options,
                                             RweaveLatexOptions)
        chunk[pos[1]] <- sub(object$syntax$docopt, "", chunk[pos[1]])
    }
    return(object)
}


RtangleFinish <- function(object, error=FALSE)
{
    if(!is.null(object$output))
        close(object$output)

    if(length(object$chunkout)>0){
        for(con in object$chunkout) close(con)
    }
}

