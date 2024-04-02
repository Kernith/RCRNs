### Simulate the progression of species concentrations in a system of reactions
### these systems are commonly referred to as chemical reaction networks

#------------------Known Issues and updates-------------------------------------

# add a check for in a concentration went negative

# when reaction and species have the same name, creates new names for everything
# New names:
# • `r.1` -> `r.1...5`
# • `r.2` -> `r.2...6`
# • `r.3` -> `r.3...7`
# • `r.1` -> `r.1...8`
# • `r.2` -> `r.2...9`
# • `r.3` -> `r.3...10`

#------------------Functions----------------------------------------------------

#' Parse a Chemical Reaction Network
#'
#' Parse a set of reactions written as equations to return the reaction
#' parameters (reactants and products, stoichiometry, and rate constant).
#'
#' Each reaction string must follow a particular format to be parsed correctly.
#' The format is "*name* = *a* *A* + *b* *B* + *...* -> *c* *C* + *d* *D* +
#' *...*, k = *const*".
#' \itemize{
#'   \item *name* is an identifier for the reaction and must be unique for each
#'     reaction in the network.
#'   \item *A*, *B*, *C*, *D*, *...* are names of chemical species acting as
#'     reactants or products in the reaction. Each species must be separated by
#'     a plus, "+", and reactants and products must be separated by a reaction
#'     arrow, "->". At least one reactant and one product must be supplied for
#'     each reaction. A species can appear on both sides of the equation as both
#'     a reactant and product (e.g. a catalyst).
#'   \item *a*, *b*, *c*, *d*, *...* are the reaction stoichiometries for
#'     species *A*, *B*, *C*, *D*, *...*, respectively. Reaction stoichiometries
#'     must be separated from the chemical species name by white space.
#'     Different stoichiometries can be supplied for a species present as both a
#'     reactant and product.
#'   \item *const* is the reaction rate constant. It can be supplied as any
#'     number that is interpretable by `as.numeric`, e.g. "50000" or "5e4", as
#'     well as scientific format, e.g. "5*10^4".
#' }
#' The units for rate constant are unspecified, but should be consistent with
#' the units of the time step and concentrations used in subsequent simulations.
#' `simulate_reaction`.
#'
#' @param rxn_strings A character vector containing the reaction equations
#'   formatted as described below
#' @param partial_orders A data frame of partial orders for the reaction. By
#'   default the reactant stoichiometries are used for the partial orders.
#'
#' @return A list of four data frames with the parameters of the reactions. The
#'   data frames returned contain:
#' \enumerate{
#'   \item rate constants
#'   \item reactant stoichiometries
#'   \item product stoichiometries
#'   \item partial orders of reaction
#' }
#'
#' @export
#'
#' @examples
#' data(CRNs)
#' parse_reactions(crns$rps$rxns)
parse_reactions <- function(rxn_strings, partial_orders = NULL) {
  
  # Check for valid supplied arguments arguments
  arg_rxns_valid <- is.character(rxn_strings)
  stopifnot("`rxn_strings` must be a character vecor." = arg_rxns_valid)
  arg_orders_valid <- is.data.frame(partial_orders) | is.null(partial_orders)
  stopifnot("`partial_orders` must be NULL or a data frame" = arg_orders_valid)
  
  
  # Filter comments
  rxns <- rxn_strings[!startsWith(rxn_strings, "#")]
  rxns <- as.data.frame(rxns)
  names(rxns) <- "rxn"
  
  # check format of rxn_strings and warn if not matching
  rxn_pattern <-
    paste0("^\\s*[[:graph:]]+\\s*=\\s*",  # reaction name
           "(",                           # 0 or more .....
           "([[:digit:]]\\s+)?",          # 0 or 1 stoich number
           "[[:graph:]]+\\s*\\+\\s*",     # reactant followed by +
           ")*",
           "([[:digit:]]\\s+)?",          # 0 or 1 stoich number
           "[[:graph:]]+\\s*->\\s*",      # reactant followed by ->
           "(",                           # 0 or more .....
           "([[:digit:]]\\s+)?",          # 0 or 1 stoich number
           "[[:graph:]]+\\s*\\+\\s*",     # products followed by +
           ")*", 
           "([[:digit:]]\\s+)?",          # 0 or 1 stoich number
           "[[:graph:]]+\\s*,\\s*",       # product followed by ,
           "k\\s*=\\s*[[:graph:]]+\\s*$", # k = const
           collapse = "")
  stopifnot("rxn_strings (excluding comments) do not match pattern.  Check strings" = grepl(rxn_pattern, rxn_strings))
  
  # separate reaction names, reactants, products, and constants
  rxns$rxn <- gsub("k = ", "", rxns$rxn)
  rxns$rxn <- gsub("=", ",", rxns$rxn)
  rxns$rxn <- gsub("->", ",", rxns$rxn)
  
  rxns <- as.data.frame(do.call(rbind, strsplit(rxns$rxn, ",")))
  rxns <- lapply(rxns, trimws)
  rxns <- data.frame(rxns)
  colnames(rxns) <- c("rxn", "react", "prod", "const")
  
  # parse k statements to get numeric rate constants
  rxns$const <- trimws(rxns$const)
  rxns$const <- gsub("10\\^", "1e", rxns$const)
  rxns$const <- as.numeric(rxns$const)
  consts <- rxns[, c(1, 4)]
  
  # Get reactant stoichiometries
  # Separate reactants and rxn names, add leading 1's, and make stoichs negative
  reacts <- strsplit(rxns$react, "\\+")
  names(reacts) <- rxns$rxn
  reacts <- stack(reacts)
  names(reacts) <- c("react", "rxn")
  reacts$react <- trimws(reacts$react)
  reacts$react <- gsub("^(\\D)", "1 \\1", reacts$react)
  reacts$react <- gsub("(^[0-9]+)", "-\\1", reacts$react)
  # separate reactant name and reactant stoich, and join to rxn name
  react_names <- as.character(reacts$rxn)
  reacts <- do.call(rbind, strsplit(reacts$react, " "))
  reacts <- cbind(reacts, react_names)
  reacts <- `colnames<-`(reacts, c("stoich", "react", "rxn"))
  reacts <- as.data.frame(reacts)
  # pivot reactant stoichs wider
  reacts <- reshape(data = reacts,
                    direction = "wide",
                    idvar = "rxn",
                    timevar = "react",
                    v.names = "stoich")
  colnames(reacts) <- sub("^stoich\\.", "", colnames(reacts))
  # replace NAs
  reacts <- as.data.frame(lapply(reacts, function(x) replace(x, is.na(x), 0)))
  # rearrange columns
  reacts <- reacts[c(1, order(colnames(reacts))[order(colnames(reacts)) != 1])]
  
  # Get product stoichiometries
  # Separate products and rxn names and add leading 1's
  prods <- strsplit(rxns$prod, "\\+")
  prods <- lapply(prods, trimws)
  prods <- lapply(prods, function (x) gsub("^(\\D)", "1 \\1", x))
  prods <- `names<-`(prods, rxns$rxn)
  prods <- stack(prods)
  prods <- `names<-`(prods, c("prod", "rxn"))
  # separate product name and product stoich, and join to rxn name
  prod_names <- as.character(prods$rxn)
  prods <- do.call(rbind, strsplit(prods$prod, " "))
  prods <- cbind(prods, prod_names)
  prods <- `colnames<-`(prods, c("stoich", "prod", "rxn"))
  prods <- as.data.frame(prods)
  # pivot product stoichs wider
  prods <- reshape(data = prods,
                   direction = "wide",
                   idvar = "rxn",
                   timevar = "prod",
                   v.names = "stoich")
  colnames(prods) <- sub("^stoich\\.", "", colnames(prods))
  # replace NAs
  prods <- as.data.frame(lapply(prods, function(x) replace(x, is.na(x), 0)))
  # rearrange columns with rxn first and species in alphabetical order
  prods <- prods[c(1, order(colnames(prods))[order(colnames(prods)) != 1])]
  
  # Add columns to reacts with species missing from table, giving them 0 stoich
  react_missing <- setdiff(colnames(prods),  colnames(reacts))
  react_missing <- prods[react_missing]
  react_missing[] <- 0
  reacts <- cbind(reacts, react_missing)
  reacts <- reacts[c(1, order(colnames(reacts))[order(colnames(reacts)) != 1])]
  reacts[2:length(reacts)] <- lapply(reacts[2:length(reacts)], as.numeric)
  # Add columns to prods with species missing from table, giving them 0 stoich
  prod_missing <- setdiff(colnames(reacts), colnames(prods))
  prod_missing <- reacts[prod_missing]
  prod_missing[] <- 0
  prods <- cbind(prods, prod_missing)
  prods <- prods[c(1, order(colnames(prods))[order(colnames(prods)) != 1])]
  prods[2:length(prods)] <- lapply(prods[2:length(prods)], as.numeric)
  
  # get partial orders from react_stoichs (by default) or input file
  orders <-
    if (is.null(partial_orders)) {
      cbind(reacts[1], lapply(reacts[2:length(reacts)], function (x) x * -1))
    } else {
      partial_orders
    }
  
  # return a named list of the reaction parameters (const, react, prod, partial)
  return(list(consts         = consts,
              react_stoich   = reacts,
              prod_stoich    = prods,
              partial_orders = orders))
  
}

#' Simulate a Chemical Reaction Network
#'
#' Simulates the progression of a network of chemical reactions from a set of
#' initial concentrations.
#'
#' Calculates reaction rates and changes in concentration for a chemical
#' reaction network described by `rxn_params`. Calculations are done using
#' classic mass-action rate equations, and are performed at a series of time
#' points incremented by `time_step` from 0 to `max time` . \cr \cr The units of
#' measure (e.g. molar or seconds) for `time_step` and `max_time` and
#' `rxn_params` (rates, rate constants, and concentrations) are arbitrary, but
#' are assumed to be consistent (e.g. if concentration is supplied in nanomolar
#' and `time_step` in seconds, rate is assumed to be in nanomolar per second).
#'
#' @param rxn_params A character vector containing strings with the reaction
#'   equations. This set of parameters can be generated from a human-friendly
#'   strings via the reaction equations with the `parse_reactions` function.
#' @param init_concs A data frame containing initial concentrations of chemical
#'   species in the reaction network. Species not supplied with an initial
#'   concentration are assumed to have an initial concentration of 0.
#' @param time_step The amount of time that elapses in each simulation cycle
#' @param max_time The full amount of time the simulation will run
#' @param output_res The number of time points that will be returned from the
#'   simulation
#'
#' @return Returns a data frame with the concentrations and the reaction rates
#'   at each output time point
#' @export
#'
#' @examples
#' data(CRNs)
#' results <- simulate_crn(
#'   rxn_params = parse_reactions(crns$rps$rxns),
#'   init_concs = crns$rps$concs
#' )
#' 
simulate_crn <- function(rxn_params, init_concs, time_step = 10^-3,
                              max_time = 60, output_res = 1000) {
  
  # print number of points in simulation and output - warn if beyond thresholds
  cat("\ntotal steps = ", max_time / time_step, ". It should be <= 10^6",
      sep = "")
  cat("\noutput a data point every ", max_time / time_step / output_res,
      " steps. It must be >= 1. It should be >=100 for smooth graph\n",
      sep = "")
  
  # read in reaction equations and initial concentrations
  rxn_params   <- rxn_params
  init_concs$init_conc <- as.double(init_concs$init_conc)
  
  # add missing species to init_concs and order species
  species <- names(rxn_params$react_stoich)[-1]
  # warn of species in supplied concs missing from rxn_params
  species_extra <- !(init_concs[[1]] %in% species)
  if(any(species_extra)) {
    stop(paste0("Chemical species {",
               paste0(init_concs[[1]][species_extra], collapse = ", "),
               "} supplied in inital concentrations but ommited from reactions\n",
               "simulation terminated"))
  }
  
  species_miss <- !(species %in% init_concs[[1]])
  if (any(species_miss)) {
    concs_miss <- data.frame("species" = species[species_miss],
                             "init_conc" = 0)
    init_concs <- rbind(init_concs, concs_miss)
    # warn of species missing from supplied concs
    cat(paste0("Chemical species {",
               paste0(species[species_miss], collapse = ", "),
               "} omitted from inital concentrations\n",
               "these species are assumed to have conc = 0."))
  }
  init_concs <- init_concs[order(init_concs$species), ]
  
  # make matrices of reaction parameters, and vector of initial concentrations
  consts_mat <- as.matrix(rxn_params$consts[-1])
  consts_mat <- `rownames<-`(consts_mat, rxn_params$consts[[1]])
  reacts_mat <- as.matrix(rxn_params$react_stoich[-1])
  reacts_mat <- `rownames<-`(reacts_mat, rxn_params$react_stoich[[1]])
  prods_mat <- as.matrix(rxn_params$prod_stoich[-1])
  prods_mat <- `rownames<-`(prods_mat, rxn_params$prod_stoich[[1]])
  orders_mat <- as.matrix(rxn_params$partial_orders[-1])
  orders_mat <- `rownames<-`(orders_mat, rxn_params$partial_orders[[1]])
  concs_vec <- init_concs[, 2]
  names(concs_vec) <- init_concs[, 1]
  
  # make transposed matrices to reduce transposition in the loop
  t_consts_mat <- t(consts_mat)
  t_prods_mat  <- t(prods_mat)
  t_reacts_mat <- t(reacts_mat)
  t_orders_mat <- t(orders_mat)
  
  # Set up the simulation loop's initial parameters
  results_colnames <- c("t", names(concs_vec), colnames(t_consts_mat))
  results <- data.frame(matrix(nrow = 0, ncol = length(results_colnames))) 
  results <- `colnames<-`(results, results_colnames)
  max_i    <- max_time / time_step
  out_step <- as.integer(max_i / output_res)
  i        <- 0
  pb       <- utils::txtProgressBar(min = 0, max = max_i, style = 3, width = 50)
  rxn_ones <- matrix(rep(1, length(consts_mat)))
  
  # Loop the steps for the duration of the simulation
  while (i <= max_i + 1) {
    
    # The reaction rate "v" for a reaction "a A + b B + ... -> P" with rate
    # constant "k" is calculated as "v = k * [A]^a' * [B]^b' * ...", where [A],
    # [B], ..., are the concentrations of the reactants and a', b', ... are the
    # partial orders, taken by default to be the reactant stoichiometries a' = a
    
    # make a matrix of values for the terms A^a, B^b, ... for each reaction
    rate_terms <- t(concs_vec ^ t_orders_mat)
    
    # calculate "v = k * A^a * B^b * ... = k * rate_terms" for each reaction
    rxn_rates <- apply(rate_terms, 1, prod) * t_consts_mat[,]
    
    # append results to output df if this time point is one of the output points
    if (i %% out_step == 0) {
      results <- rbind(results, c(t = i * time_step, concs_vec, rxn_rates))
      utils::setTxtProgressBar(pb, i)
    }
    
    # calculate new concentrations
    prods_change  <- prods_mat * rxn_rates * time_step
    reacts_change <- reacts_mat * rxn_rates * time_step
    conc_change <- t(t(prods_change + reacts_change) %*% rxn_ones)
    concs_vec <- concs_vec + as.numeric(conc_change)
    
    # add a check for in a concentration went negative
    # check check check
    stopifnot(
      "Concentrations went below 0. Consider reducing the timestep." = 
      concs_vec >= 0
    )
    
    # iterate for next loop
    i <- i + 1
    
  }
  
  results <- `colnames<-`(results, results_colnames)
  return(results)
  
}
