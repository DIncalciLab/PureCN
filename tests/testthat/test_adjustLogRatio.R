context("adjustLogRatio")

test_that("Function returns expected values for example coverage", {
     data(purecn.example.output)
     log.ratio <- purecn.example.output$results[[1]]$seg$seg.mean   
     purity <- purecn.example.output$results[[1]]$purity
     ploidy <- purecn.example.output$results[[1]]$ploidy
     log.ratio.adjusted <- adjustLogRatio(log.ratio, purity, ploidy)
     total.ploidy <- 1.73
     p <- 1
     log.ratio.offset <- 0
     opt.C <- (2^(log.ratio.adjusted + log.ratio.offset) *  total.ploidy)/p - ((2 * (1 - p))/p)
     expect_lt(abs(min(log.ratio.adjusted, na.rm=TRUE) + 8), 0.001)
     expect_lt(median(abs(opt.C - purecn.example.output$results[[1]]$seg$C)), 0.1)  
})

