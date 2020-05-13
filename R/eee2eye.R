#' Calculate the E:I for a Steady-State Lake using 18O-H2O
#'
#' This function allows you to calculate the evaporation to inflow (E:I) ratio of a lake based on a few assumptions and its water isotope value.
#' @param df dataframe with your values
#' @param dL steady-state lake isotope value (‰)
#' @param dI input/source water isotope value (‰)
#' @param dP precipitation isotope value (‰)
#' @param temp temperature at surface of lake (C)
#' @param humid relative humidity (decimal)
#' @param evap annual evaporation rate of area (m/yr)
#' @param k influence of seasonality (0.5: highly seasonal; 1: no seasonality)
#' @keywords water isotopes, E:I, isotope mass balance
#' @export
#' @examples
#' eee2eye(ei_input, 'd18O_H2O_permille', 'dI_permille', 'dP_permille', 'temp_C', 'h_dec', 'evap_m.yr', 'k')


eee2eye <- function(df, dL, dI, dP, temp, humid, evap, k){

  #Calculate isotope fractionation and separation factors:
  alpha_plus <- exp((-7.685/(10^3)) + (6.7123/(273.15 + df[[temp]])) - (1666.4/((273.15 + df[[temp]])^2)) + (350410/((273.15 + df[[temp]])^3)))
  e_plus_permille <- (alpha_plus - 1) * 1000

  theta <- 1
  Ck_permille <- 14.2
  ek_permille <- theta * Ck_permille * (1 - df[[humid]])

  #Calculate dA
  dA_permille <- (df[[dP]] - (df[[k]] * e_plus_permille)) / (1 + ((10^-3) * df[[k]] * e_plus_permille))

  #Calculate temporal enrichment slope (m)
  m <- (df[[humid]] - (10^-3) * (ek_permille + (e_plus_permille / alpha_plus))) / (1 - df[[humid]] + ((10^-3) * ek_permille))

  #Calculate limiting isotopic composition
  dstar_permille <- ((df[[humid]] * dA_permille) + ek_permille + (e_plus_permille/alpha_plus)) / (df[[humid]] - ((10^-3) * (ek_permille + (e_plus_permille/alpha_plus))))

  #Calculate our E:I based on lake sig
  df$E.I <- (df[[dL]] - df[[dI]]) / (m * (dstar_permille - df[[dL]]))

  return(df)

}
