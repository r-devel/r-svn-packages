\name{RweaveLatex}
\alias{RweaveLatex}
\alias{RweaveLatexSetup}
\title{R/LaTeX Driver for Sweave}
\description{
  A driver for \code{\link{Sweave}} that translates R code chunks in
  LaTeX files.
}
\usage{
RweaveLatex()

RweaveLatexSetup(file, syntax, output = NULL, quiet = FALSE,
                 debug = FALSE, echo = TRUE, eval = TRUE,
                 keep.source = FALSE, split = FALSE,
                 pdf = TRUE, eps = FALSE, envir=.GlobalEnv)
}
\arguments{
  \item{file}{Name of Sweave source file.}
  \item{syntax}{An object of class \code{SweaveSyntax}.}
  \item{output}{Name of output file, default is to remove extension
    \file{.nw}, \file{.Rnw} or \file{.Snw} and to add extension
    \file{.tex}. Any 
    directory names in \code{file} are also removed such that the output
    is created in the current working directory.}
  \item{quiet}{If \code{TRUE} all progress messages are suppressed.}
  \item{debug}{If \code{TRUE}, input and output of all code
    chunks is copied to the console.}
  \item{echo}{set default for option \code{echo}, see details below.}
  \item{eval}{set default for option \code{eval}, see details below.}
  \item{keep.source}{set default for option \code{keep.source}, see
    details below.}
  \item{split}{set default for option \code{split}, see details below.}
  \item{pdf}{set default for option \code{pdf}, see details below.}
  \item{eps}{set default for option \code{eps}, see details below.}
  \item{envir}{Environment where code chunks are evaluated. Either an
    environment, or the character string \code{"new"}, in which case a
    new temprorary environment is used.}
}
\section{Supported Options}{
  RweaveLatex supports the following options for code chunks (the values
  in parentheses show the default values):
  \describe{
    \item{echo:}{logical (\code{TRUE}). Include S code in the
      output file?}
    \item{keep.source:}{logical (\code{FALSE}).  When echoing, if \code{keep.source == TRUE}
      the original source is copied to the file.  Otherwise, deparsed source is echoed.}
    \item{eval:}{logical (\code{TRUE}). If \code{FALSE}, the code chunk is not
      evaluated, and hence no text or graphical output produced.}
    \item{results:}{character string (\code{verbatim}).
      If \code{verbatim}, the output of S commands is
      included in the verbatim-like Soutput environment. If
      \code{tex}, the output is taken to be already proper latex markup
      and included as is. If \code{hide} then all output is
      completely suppressed (but the code executed during the weave).}
    \item{print:}{logical (\code{FALSE})
      If \code{TRUE}, each expression in the
      code chunk is wrapped into a \code{print()} statement before evaluation,
      such that the values of all expressions become visible.}
    \item{term:}{logical (\code{TRUE}). If \code{TRUE}, visibility of values
      emulates an interactive R session: values of assignments are not
      printed, values of single objects are printed. If \code{FALSE},
      output comes only from explicit \code{\link{print}} or
      \code{\link{cat}} statements.}
    \item{split:}{logical (\code{FALSE}). If \code{TRUE}, text output is
      written to separate files 
      for each code chunk.}
    \item{strip.white:}{character string (\code{false}). If \code{true}, blank
      lines at the beginning and end of output are removed. If
      \code{all}, then all blank lines are removed from the output.}
    \item{prefix:}{logical (\code{TRUE}). If \code{TRUE} generated filenames of
      figures and output have a common prefix.}
    \item{prefix.string:}{a character string, default is the name of the
      \file{.Snw} source file.}
    \item{include:}{logical (\code{TRUE}), indicating whether input
      statements for text output and includegraphics statements
      for figures should be auto-generated. Use \code{include = FALSE} if
      the output should appear in a different place than the code chunk
      (by placing the input line manually).}
    \item{fig:}{logical (\code{FALSE}), indicating whether the code
      chunk produces 
      graphical output. Note that only one figure per code chunk can be
      processed this way.}
    \item{eps:}{logical (\code{FALSE}), indicating whether EPS figures shall be
      generated. Ignored if \code{fig = FALSE}.} 
    \item{pdf:}{logical (\code{TRUE}), indicating whether PDF figures shall be
      generated. Ignored if \code{fig = FALSE}.}
    \item{width:}{numeric (6), width of figures in inch.}
    \item{height:}{numeric (6), height of figures in inch.}
    \item{expand:}{logical (\code{TRUE}).  Expand references to other chunks
      so that only R code appears in the output file.  If \code{FALSE}, the
      chunk reference (e.g. \code{<<chunkname>>}) will appear.  The 
      \code{expand=FALSE} option requires \code{keep.source = TRUE} or it 
      will have no effect.}
    \item{concordance:}{logical (\code{FALSE}).  Write a concordance file
      to link the input line numbers to the output line numbers.  This is
      an experimental feature; see the source code for the output format,
      which is subject to change in future releases.}
  }
}
\author{Friedrich Leisch}
\references{
  Friedrich Leisch: Sweave User Manual, 2002\cr
  \url{http://www.ci.tuwien.ac.at/~leisch/Sweave}
}
\seealso{\code{\link{Sweave}}, \code{\link{Rtangle}}}
\keyword{utilities}
