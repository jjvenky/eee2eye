#' Calculate the E:I for a Steady-State Lake using 18O-H2O
#'
#' This function allows you to calculate the evaporation to inflow (E:I) ratio of a lake based on a few assumptions and its water isotope value.
#' @param df dataframe with your values
#' @param dL steady-state lake isotope value (‰)
#' @param dI input/source water isotope value (‰)
#' @param dP precipitation isotope value (‰)
#' @param temp temperature at surface of lake (C)
#' @param humid relative humidity (decimal)
#' @param k influence of seasonality (0.5: highly seasonal; 1: no seasonality)
#' @keywords water isotopes, E:I, isotope mass balance
#' @export
#' @examples
#' ei_input <- data.frame(dL_permille = c(-11.77, -15.67, -18.23),dI_permille = c(-20.7, -18.2, -20.2),dP_permille = c(-23, -28, -32),temp_C = c(14.3, 12.1, 8.9),h_dec = c(0.68, 0.71, 0.58),k = c(0.7, 0.72, 0.65))
#' eee2eye(ei_input, 'dL_permille', 'dI_permille', 'dP_permille', 'temp_C', 'h_dec', 'k')


eee2eye <- function(df, dL, dI, dP, temp, humid, k){

  #check if using values or specific columns for each parameter
  dL_input <- if(is.numeric(dL) == TRUE) {
    dL
  } else {
    df[[dL]]
  }

  dI_input <- if(is.numeric(dI) == TRUE) {
    dI
  } else {
    df[[dI]]
  }

  dP_input <- if(is.numeric(dP) == TRUE) {
    dP
  } else {
    df[[dP]]
  }

  temp_input <- if(is.numeric(temp) == TRUE) {
    temp
  } else {
    df[[temp]]
  }

  humid_input <- if(is.numeric(humid) == TRUE) {
    humid
  } else {
    df[[humid]]
  }

  k_input <- if(is.numeric(k) == TRUE) {
    k
  } else {
    df[[k]]
  }

  #Calculate isotope fractionation and separation factors:
  alpha_plus <- exp((-7.685/(10^3)) + (6.7123/(273.15 + temp_input)) - (1666.4/((273.15 + temp_input)^2)) + (350410/((273.15 + temp_input)^3)))
  e_plus_permille <- (alpha_plus - 1) * 1000

  theta <- 1
  Ck_permille <- 14.2
  ek_permille <- theta * Ck_permille * (1 - humid_input)


  #Calculate dA
  dA_permille <- (dP_input - (k_input * e_plus_permille)) / (1 + ((10^-3) * k_input * e_plus_permille))

  #Calculate temporal enrichment slope (m)
  m <- (humid_input - (10^-3) * (ek_permille + (e_plus_permille / alpha_plus))) / (1 - humid_input + ((10^-3) * ek_permille))

  #Calculate limiting isotopic composition
  dstar_permille <- ((humid_input * dA_permille) + ek_permille + (e_plus_permille/alpha_plus)) / (humid_input - ((10^-3) * (ek_permille + (e_plus_permille/alpha_plus))))

  #Calculate our E:I based on lake sig
  df$E.I <- (dL_input - dI_input) / (m * (dstar_permille - dL_input))

  return(df)

}
