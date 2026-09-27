# jamovi serialises every analysis after each run; the headless wrappers never do, so an error
# there (for example a results item named like a jmvcore member, such as `path`) passes every
# other test and only shows in jamovi. Build each analysis as jamovi does and serialise it.
serialise <- function(cls, opts, data = data.frame()) {
    an <- cls$new(options = opts, data = data, datasetId = "d", analysisId = "1", revision = 1)
    an$init(); an$run()
    RProtoBuf::serialize(an$asProtoBuf(final = FALSE, includeState = TRUE), NULL)
}

test_that("every analysis serialises like in jamovi, in design and data modes", {
    jmvcore:::initProtoBuf()
    sq <- read.csv(test_path("..", "..", "data", "sprt_sequential.csv"))
    sq$arm <- factor(sq$arm)
    means <- read.csv(test_path("..", "..", "data", "gsd_interim_means.csv"))
    surv <- read.csv(test_path("..", "..", "data", "maxsprt_surveillance.csv"))
    ok <- function(x) expect_true(is.raw(x) && length(x) > 0)
    ok(serialise(gsdClass, gsdOptions$new()))
    ok(serialise(gsdClass, gsdOptions$new(lookN = "n", estimate = "estimate", se = "se"), means))
    ok(serialise(phase2Class, phase2Options$new(analyse = TRUE, obsN = 11, obsX = 6, boundPlot = TRUE, outcomePlot = TRUE)))
    ok(serialise(phase2Class, phase2Options$new(design = "single", outcomePlot = TRUE)))
    ok(serialise(sprtClass, sprtOptions$new()))
    ok(serialise(sprtClass, sprtOptions$new(dep = "response", pathTable = TRUE, ocPlot = TRUE), sq))
    ok(serialise(maxsprtClass, maxsprtOptions$new()))
    ok(serialise(maxsprtClass, maxsprtOptions$new(cases = "cases", expected = "expected"), surv))
    ok(serialise(bayesmonClass, bayesmonOptions$new(lookN = "n", estimate = "estimate", se = "se", densityPlot = TRUE, ocPlot = TRUE), means))
    ok(serialise(bayesmonClass, bayesmonOptions$new(model = "bf", dep = "score", group = "arm", nsim = 100), sq))
    ok(serialise(anytimeClass, anytimeOptions$new(dep = "score", group = "arm"), sq))
})
