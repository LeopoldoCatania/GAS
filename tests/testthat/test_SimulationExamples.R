test_that('documented simulation examples initialize at their stated targets', {
  for (topic in c('UniGASSim','MultiGASSim')) {
    env <- new.env(parent=globalenv())
    invisible(capture.output(utils::example(topic, package='GAS', local=env,
                                           echo=FALSE, ask=FALSE, character.only=TRUE)))
    target <- if(topic=='UniGASSim') c(0,1.5,7) else env$Theta
    expect_equal(as.numeric(env$Sim@GASDyn$mTheta[,1]),target,tolerance=1e-10)
    expect_true(all(is.finite(env$Sim@GASDyn$mTheta)))
    expect_true(all(is.finite(env$Sim@GASDyn$vLLK)))
    if(topic=='MultiGASSim') {
      expect_true(all(is.finite(env$Sim@GASDyn$mY)))
      # The old example reached the numerical scale clamps and almost
      # singular correlation matrices, failing on Windows R-devel.
      scales <- env$Sim@GASDyn$mTheta[4:6,]
      expect_gt(min(scales),1e-5)
      expect_lt(max(scales),1e5)
    }
  }
})
