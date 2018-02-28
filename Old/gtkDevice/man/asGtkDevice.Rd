\name{asGtkDevice}
\alias{asGtkDevice}
\title{Use Gtk drawing area widget as new R graphics device}
\description{

  This allows one to create a new R Gtk graphics device that
  draws on the specified Gtk drawing area widget. The widget can be
  created separately from the device and hence embedded within an
  arbitrary graphical interface.
  
}
\usage{
asGtkDevice(widget, width=300, height=300, pointsize=12)
}
\arguments{
  \item{widget}{an object of class \code{GtkDrawingArea}
    typically created using the \link[RGtk]{gtkDrawingArea} function
    in the \url{http://www.omegahat.org/RGtk} package.
    This should have a simple structure, being an external pointer
    object containing the address of the \code{GtkWidget}
    and having a class vector attribute containing \code{"GtkDrawingArea"}.
    }
  \item{width}{a value giving a `hint' for the number of pixels in the
    horizontal dimension.  The device will typically determine this when it first draws on the
    widget. Note this is in pixels, not inches as with the \code{\link{gtk}} device.
  }
  \item{height}{a value giving a `hint' for the number of pixels in the
    vertical dimension.
    The device will typically determine this when it first draws on the
    widget.}
  \item{pointsize}{the default pointsize to be used.}
}
\details{
 
}
\value{
  Currently, a logical value indicating whether the
  creation of the new device was successful or not.
}
\references{\url{http://www.gtk.org},
  \url{http://www.omegahat.org/RGtk}
  \url{http://www.omegahat.org/gtk}
}
\author{Duncan Temple Lang}
\note{
  It may be possible to compile this on Windows
   in the near future.
}


\seealso{
  \code{\link{gtk}}
  \code{\link[RGtk]{gtkDrawingArea}}  
  \code{\link{x11}}  
}

\examples{
\dontrun{
 library(RGtk)
 win <- gtkWindow(show=FALSE)
 d <- gtkDrawingArea()
 asGtkDevice(d)
 win$Add(d)
 plot(rnorm(100))
}
}
\keyword{device}
\keyword{interface}
