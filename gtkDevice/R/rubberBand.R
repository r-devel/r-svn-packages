rubberBand =
function(drawingArea, color = c(65535, 0, 0), lwd = 1, lty = 0, connect = TRUE)
{
  if(!inherits(drawingArea, "GtkDrawingArea"))
    stop("drawingArea must be a GtkDrawingArea widget")

  if(is.character(color))
    color = as.integer(col2rgb(color)*65535/255)

  lwd = as.integer(max(c(lwd, 1)))

  lty = as.integer(min(max(lty, 0), 2))
  
  .Call("R_initRubberBand", drawingArea, lwd, lty, color, as.logical(connect), PACKAGE = "gtkDevice")
}

getRubberBandCoordinates =
function(info)
{
 ans = .Call("R_getRubberBandCoordinates", info, PACKAGE = "gtkDevice")

 ans = matrix(ans, 2, 2, dimnames = list(c("x", "y"), c("start", "end")))
 attr(ans, "unit") <- "pixel"
 
 ans
}

disconnectRubberBand =
function(info)
{
  .Call("R_connectRubberBand", info, FALSE, PACKAGE = "gtkDevice")
}

connectRubberBand =
function(info)
{
  .Call("R_connectRubberBand", info, TRUE, PACKAGE = "gtkDevice")
}  
