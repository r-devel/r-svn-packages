## self-starting model function for the Richards growth model
SSRichards <-
    selfStart(~ Asym * (1+exp((xmid - input)/scal))^(-exp(-lpow)),
              function(mCall, data, LHS)
          {
              linit <- unname(attr(SSlogis, "initial")(mCall, data, LHS))
              xy <- data.frame(sortedXyData(mCall[["input"]], LHS, data))
              if (nrow(xy) < 5) {
                  stop("too few distinct input values to fit a logistic model")
              }
              pars <- unname(coef(nls(y ~ (1+exp((xmid - x)/scal))^(-exp(-lpow)),
                                      xy, c(xmid = linit[2], scal = linit[3],
                                            lpow = 0.001), alg = "plinear")))
              value <- pars[c(4,1:3)]
              names(value) <- mCall[c("Asym", "xmid", "scal", "lpow")]
              value
          },
              c("Asym", "xmid", "scal", "lpow"),
              function(input, Asym, xmid, scal, lpow) {})

