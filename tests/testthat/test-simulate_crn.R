test_that("simple addition works", {
  # setup
  crn <- c("rxn = A + 2 B -> 3 C + D, k = 0.002")
  parm <- parse_reactions(crn)
  concs <- data.frame(species = c("A", "B", "C"),
                      init_conc = c(100, 80, 0))
  results <- simulate_crn(parm, concs, time_step = 10^-3,
               max_time = 6, output_res = 10)
  
  # tests
  expect_equal(as.numeric(results[1, ]),
               c(0, 100, 80, 0, 0, 1280))
  expect_equal(as.numeric(results[11, ]),
               c(6.0000, 60.3357, 0.6714, 118.9928, 39.6643, 0.0544),
               tolerance = 1.0001)
  expect_equal(names(results),
               c("t", "A", "B", "C", "D", "rxn"))
})

test_that("multiple reactions works", {
  # setup
  crn <- c("rxn1 =   A + 2 B -> 3 C + D, k = 0.002",
           "rxn2 = 2 C + D   ->   E,     k = 0.003",
           "rxn3 =   E       ->   F + G, k = 0.001")
  parm <- parse_reactions(crn)
  concs <- data.frame(species = c("A", "B", "C", "E"),
                      init_conc = c(100, 80, 0, 5))
  results <- simulate_crn(parm, concs, time_step = 10^-3,
                          max_time = 6, output_res = 10)
  
  # tests
  expect_equal(as.numeric(results[1, ]),
               c(0, 100, 80, 0, 0, 5, 0, 0, 1280, 0, 0.005))
  expect_equal(as.numeric(results[11, ]),
               c(6.0, 60.3357, 0.6714, 39.6891, 0.0124, 44.3968, 0.2550, 0.2550, 0.0544, 0.05863, 0.0443),
               tolerance = 1.0001)
  expect_equal(names(results),
               c("t", "A", "B", "C", "D", "E", "F",
                 "G", "rxn1", "rxn2", "rxn3"))
})

test_that("differently formatted reactions works", {
  # setup
  crn <- c("rxn.1 =   glucose + 2 B.12      -> 3 ATP      + D, k = 0.002",
           "rxn_2 = 2 ATP     +   D-alanine  ->   E,            k = 10^-3",
           "rxn   =   E                     ->   Factor_X + G, k = 5e-3")
  parm <- parse_reactions(crn)
  concs <- data.frame(species = c("glucose", "B.12", "ATP", "E"),
                      init_conc = c(100, 80, 0, 5))
  results <- simulate_crn(parm, concs, time_step = 10^-3,
                          max_time = 6, output_res = 10)
  
  # tests
  expect_equal(as.numeric(results[1, ]),
               c(0, 0, 80, 0, 0, 5, 0, 0, 100, 1280, 0, 0.025))
  expect_equal(as.numeric(results[11, ]),
               c(6.0000, 118.9928,  0.6714, 39.6643, 0.0000, 4.8522,
                 0.1477,   0.1478, 60.3357,  0.0544, 0.0000, 0.0242),
               tolerance = 1.0001)
  expect_equal(names(results),
               c("t", "ATP", "B.12", "D", "D.alanine", "E", "Factor_X",
                 "G", "glucose", "rxn.1", "rxn_2", "rxn"))
})
