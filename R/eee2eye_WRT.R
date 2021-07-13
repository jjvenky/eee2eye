#' Calculate the Water Residence Time from E:I for a Steady-State Lake
#'
#' This function allows you to calculate the water residence time of a lake based on an evaporation to inflow (E:I) ratio.
#' @param df dataframe with your values
#' @param E.I E:I ratio
#' @param e annual evaporation depth (m/yr)
#' @param SA lake surface area (m2)
#' @param V lake volume (m3)
#' @keywords water residence time, water isotopes, E:I, isotope mass balance
#' @export
#' @examples
#' ei_input <- data.frame(E.I = c(0.2042, 0.3138, 0.1838),
#'                        e_myr = c(0.3965, 0.3965, 0.3965),
#'                        SA_m2 = c(315900, 300825, 589950),
#'                        V_m3 = c(2466000, 3004064, 5712829))
#' eee2eye_WRT(ei_input, 'E.I', 'e_myr', 'SA_m2', 'V_m3')


eee2eye_WRT <- function(df, E.I, e, SA, V){

  #check if using values or specific columns for each parameter
  E.I_input <- if(is.numeric(E.I) == TRUE) {
    E.I
  } else {
    df[[E.I]]
  }

  e_input <- if(is.numeric(e) == TRUE) {
    e
  } else {
    df[[e]]
  }

  SA_input <- if(is.numeric(SA) == TRUE) {
    SA
  } else {
    df[[SA]]
  }

  V_input <- if(is.numeric(V) == TRUE) {
    V
  } else {
    df[[V]]
  }

  #calculate evaporation rate (m3/yr)
  E <- e_input * SA_input

  #calculate WRT
  df$WRT_yr <- E.I_input * V_input / E

  return(df)

  }
