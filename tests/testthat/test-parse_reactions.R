
#------------------tests for right string pattern-------------------------------

test_that("simple reaction works", {
  crn <- c("# A simple addition reaction", "rxn = A + B -> C, k = 0.002")
  expect_no_error(parse_reactions(crn))
})

test_that("catalyzed reaction works", {
  crn <- c("# A simple addition reaction", "rxn = A + B -> B + C, k = 0.002")
  expect_no_error(parse_reactions(crn))
})

test_that("long species names works", {
  crn <- c("# A simple addition reaction",
           "glucoseOxidation = glucose + oxygen + glucoseOxidase -> gluconolactone + glucoseOxidase, k = 0.002")
  expect_no_error(parse_reactions(crn))
})

test_that("special char species names works", {
  crn <- c("# A simple addition reaction",
           "glucose_Oxidation = D-glucose + di.oxygen + glucose_Oxidase -> glucono.lactone + glucose_Oxidase, k = 0.002")
  expect_no_error(parse_reactions(crn))
})

test_that("multiple reactions work", {
  crn <- c("# A simple addition reaction",
           "rxn1 = A + B -> C, k = 0.002",
           "rxn2 = C + D -> E, k = 0.002",
           "rxn3 = A + E -> F, k = 0.002")
  expect_no_error(parse_reactions(crn))
})

#------------------tests for wrong string pattern-------------------------------

# missing parameters
test_that("pattern: missing names", {
  crn <- c("# A simple addition reaction", "A + B -> C, k = 0.002")
  expect_error(parse_reactions(crn))
})

test_that("pattern: missing reactants", {
  crn <- c("# A simple addition reaction", "rxn =  -> C, k = 0.002")
  expect_error(parse_reactions(crn))
})

test_that("pattern: missing products", {
  crn <- c("# A simple addition reaction", "rxn = A + B -> , k = 0.002")
  expect_error(parse_reactions(crn))
})

test_that("pattern: missing const", {
  crn <- c("# A simple addition reaction", "rxn = A + B -> C")
  expect_error(parse_reactions(crn))
})


# parameter separators
test_that("space in reaction name fails", {
  crn <- c("# A simple addition reaction", "best reaction = A + B -> C, k = 0.002")
  expect_error(parse_reactions(crn))
})

test_that("no = after name fails", {
  crn <- c("# A simple addition reaction", "rxn A + B -> C, k = 0.002")
  expect_error(parse_reactions(crn))
})

test_that("no -> between Rs and Ps fails", {
  crn <- c("# A simple addition reaction", "rxn = A + B C, k = 0.002")
  expect_error(parse_reactions(crn))
})

test_that("no , after products fails", {
  crn <- c("# A simple addition reaction", "rxn = A + B -> C k = 0.002")
  expect_error(parse_reactions(crn))
})

# pluses and stoichiometry
test_that("hanging + in reactants fails", {
  crn <- c("# A simple addition reaction", "rxn = A + B + -> C, k = 0.002")
  expect_error(parse_reactions(crn))
})

test_that("hanging + in products fails", {
  crn <- c("# A simple addition reaction", "rxn = A + B -> C +, k = 0.002")
  expect_error(parse_reactions(crn))
})

test_that("letter as reactant+ stoich fails", {
  crn <- c("# A simple addition reaction", "rxn = a A + B -> C, k = 0.002")
  expect_error(parse_reactions(crn))
})

test_that("letter as product+ stoich fails", {
  crn <- c("# A simple addition reaction", "rxn = A + B -> c C + D, k = 0.002")
  expect_error(parse_reactions(crn))
})

test_that("letter as reactant-> stoich fails", {
  crn <- c("# A simple addition reaction", "rxn = A + b B -> C, k = 0.002")
  expect_error(parse_reactions(crn))
})

test_that("letter as product-> stoich fails", {
  crn <- c("# A simple addition reaction", "rxn = A + B -> C + d D, k = 0.002")
  expect_error(parse_reactions(crn))
})

test_that("extra plusses fails", {
  crn <- c("# A simple addition reaction", "rxn = A + + B -> C, k = 0.002")
  expect_error(parse_reactions(crn))
})


#------------------tests for rate constant--------------------------------------

test_that("long format constants work", {
  crn <- c("# A simple addition reaction", "rxn = A + B -> C, k = 0.002")
  expect_equal(parse_reactions(crn)$consts[[2]], 0.002)
})

test_that("10^Y format constants work", {
  crn <- c("# A simple addition reaction", "rxn = A + B -> C, k = 0.005")
  expect_equal(parse_reactions(crn)$consts[[2]], 0.005)
})

test_that("1eY format constants work", {
  crn <- c("# A simple addition reaction", "rxn = A + B -> C, k = 0.006")
  expect_equal(parse_reactions(crn)$consts[[2]], 0.006)
})


#------------------tests for stoich---------------------------------------------

test_that("omitted reactant stoich = 1", {
  crn <- c("# A simple addition reaction", "rxn = A + 2 B -> 3 C + D, k = 0.002")
  expect_equal(parse_reactions(crn)$react_stoich[[2]], -1)
})

test_that("supplied reactant stoich is retained", {
  crn <- c("# A simple addition reaction", "rxn = A + 2 B -> 3 C + D, k = 0.002")
  expect_equal(parse_reactions(crn)$react_stoich[[3]], -2)
})

test_that("react_stoich for non-reactant = 0", {
  crn <- c("# A simple addition reaction", "rxn = A + 2 B -> 3 C + D, k = 0.002")
  expect_equal(parse_reactions(crn)$react_stoich[[4]], 0)
})

test_that("omitted product stoich = 1", {
  crn <- c("# A simple addition reaction", "rxn = A + 2 B -> 3 C + D, k = 0.002")
  expect_equal(parse_reactions(crn)$prod_stoich[[5]], 1)
})

test_that("supplied product stoich is retained", {
  crn <- c("# A simple addition reaction", "rxn = 1 A + 2 B -> 3 C + D, k = 0.002")
  expect_equal(parse_reactions(crn)$prod_stoich[[4]], 3)
})

test_that("product_stoich for non-product = 0", {
  crn <- c("# A simple addition reaction", "rxn = A + 2 B -> 3 C + D, k = 0.002")
  expect_equal(parse_reactions(crn)$prod_stoich[[2]], 0)
})


#------------------tests for partial orders-------------------------------------

test_that("partial_orders is the negative of react_stoich", {
  crn <- c("# A simple addition reaction", "rxn = A + 2 B -> 3 C + D, k = 0.002")
  expect_equal(as.numeric(parse_reactions(crn)$react_stoich[1, 2:4]),
               -1 * as.numeric(parse_reactions(crn)$partial_orders[1, 2:4]))
})
