
parsed_simple <- list(
  consts = data.frame(
    rxn = "reaction",
    const = 0.001
  ),
  react_stoich = data.frame(
    rxn = "reaction",
    A = -1,
    B = -1,
    C = 0
  ),
  prod_stoich = data.frame(
    rxn = "reaction",
    A = 0,
    B = 0,
    C = 1
  ),
  partial_orders = data.frame(
    rxn = "reaction",
    A = 1,
    B = 1,
    C = 0
  )
)

parse_reactions(simple)$consts[[1]] == parsed_simple$consts[[1]]

# tests for rxn const
test_that("long format constants work", {
  simple <- c("# A simple addition reaction",
              "reaction = A + B -> C, k = 10^-3")
  results <- parse_reactions(simple)
  
  expect_equal(results$consts[[2]], 0.01),
  
})




