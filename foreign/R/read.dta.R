read.dta <- function(file, convert.dates=TRUE,tz="GMT",convert.factors=TRUE){
    rval<-.External("do_readStata", file, PACKAGE = "foreign")
    if (convert.dates){
        ff<-attr(rval,"formats")
        dates<-grep("%_*d",ff)
        for(v in dates)
            rval[[v]]<-ISOdate(1960,1,1,tz=tz)+24*60*60*rval[[v]]
    }
    if (convert.factors){
        ll<-attr(rval,"val.labels")
        tt<-attr(rval,"label.table")
        factors<-which(ll!="")
        for(v in factors)
            rval[[v]]<-factor(rval[[v]],levels=tt[[ll[v]]],labels=names(tt[[ll[v]]]))
    }
        
    rval
}

write.dta <- function(dataframe, file, version = 6,convert.dates=TRUE,tz="GMT") {
    if (convert.dates){
        dates<-which(sapply(data.frame,function(x) inherits(x,"POSIXt")))
        for( v in dates)
            dataframe[[v]]<-round(julian(dataframe[[v]],ISOdate(1960,1,1,tz=tz)))
    }
    if (any(sapply(dataframe, function(x) !is.null(dim(x)))))
        stop("Can't handle multicolumn columns")
    invisible(.External("do_writeStata", file, dataframe, version,
                        PACKAGE = "foreign"))
}
