#
# Read and clean data from Stata files
#
clean_data <- function(){
  #
  # read data from Stata files
  data_1980s <- read.dta("data/data1980s_share.dta")
  data_regulation_share <- read.dta("data/data_regulation_share.dta")
  #
  data_1980s$islands <- ifelse(data_1980s$n_islands == 2, 0, 1)
  data_1980s$landlock <- ifelse(data_1980s$n_landlock == 2, 0, 1)
  data_1980s$trade <- ifelse(is.na(data_1980s$ln_trade), 0, 1)
  #
  # Define continents
  #
  countries_africa  <- read.csv(file="data/countries/Africa.csv")
  countries_asia    <- read.csv(file="data/countries/Asia.csv")
  countries_europe  <- read.csv(file="data/countries/Europe.csv")
  countries_NA      <- read.csv(file="data/countries/North_America.csv")
  countries_SA      <- read.csv(file="data/countries/South_America.csv")
  countries_oceania <- read.csv(file="data/countries/Oceania.csv")
  #
  unique_countries <- as.character(unique(data_1980s[,"impcode"]))
  n_countries <- length(unique_countries)
  which_continent <- rep("unassigned",n_countries)
  #
  for(i in 1:n_countries){
    country_i <- unique_countries[i]
    #
    test_africa <- any(country_i == countries_africa)
    if(test_africa==TRUE){
      which_continent[i] <- "AFRICA"
    }else{
      test_asia <- any(country_i == countries_asia)
      if(test_asia==TRUE){
        which_continent[i] <- "ASIA"
      }else{
        test_europe <- any(country_i == countries_europe)
        if(test_europe==TRUE){
          which_continent[i] <- "EUROPE"
        }else{
          test_NA <- any(country_i == countries_NA)
          if(test_NA==TRUE){
            which_continent[i] <- "NORTH_AMERICA"
          }else{
            test_SA <- any(country_i == countries_SA)
            if(test_SA==TRUE){
              which_continent[i] <- "SOUTH_AMERICA"
            }else{
              test_oceania <- any(country_i == countries_oceania)
              if(test_oceania==TRUE){
                which_continent[i] <- "OCEANIA"
              }
            }
          }
        }
      }
    }
  }
  #
  # Deal with unassigned values
  #
  which_continent[unique_countries=="LIBY ARAB JM"] <- "AFRICA"
  which_continent[unique_countries=="CENTRAL AFR. REP."] <- "AFRICA"
  which_continent[unique_countries=="ZAIRE"] <- "AFRICA"
  which_continent[unique_countries=="EQ. GUINEA"] <- "AFRICA"
  which_continent[unique_countries=="COTE D'IVOIRE"] <- "AFRICA"
  which_continent[unique_countries=="GUINEA-BISSAU (includes CAPE VERDE)"] <- "AFRICA"
  which_continent[unique_countries=="EQ. GUINEA"] <- "AFRICA"
  which_continent[unique_countries=="UNTD RP TANZANIA"] <- "AFRICA"
  which_continent[unique_countries=="WESTERN SAHARA"] <- "AFRICA"
  #
  which_continent[unique_countries=="SYRN ARAB RP"] <- "ASIA"
  which_continent[unique_countries=="UNTD ARAB EM"] <- "ASIA"
  which_continent[unique_countries=="MYANMAR (BURMA)"] <- "ASIA"
  which_continent[unique_countries=="INDONESIA (includes MACAU)"] <- "ASIA"
  which_continent[unique_countries=="KOREA RP (SOUTH)"] <- "ASIA"
  which_continent[unique_countries=="LAOS P.DEM.R"] <- "ASIA"
  which_continent[unique_countries=="KOREA D P RP (NORTH)"] <- "ASIA"
  which_continent[unique_countries=="FM USSR"] <- "ASIA"
  #
  which_continent[unique_countries=="BELGIUM-LUX."] <- "EUROPE"
  which_continent[unique_countries=="CZECHOSLOVAKIA"] <- "EUROPE"
  which_continent[unique_countries=="DENMARK (includes FAROE ISLDS)"] <- "EUROPE"
  which_continent[unique_countries=="FM YUGOSLAVIA (includes CROATIA, SLOVENIA)"] <- "EUROPE"
  #
  which_continent[unique_countries=="SOLOMON ISLDS"] <- "OCEANIA"
  which_continent[unique_countries=="KIRIBATI (includes SOLOMON ISLDS, TONGA, TUVALU)"] <- "OCEANIA"
  which_continent[unique_countries=="NEW CALEDONIA (includes FR POLYNESIA, VANUATA)"] <- "OCEANIA"
  which_continent[unique_countries=="PAPUA N.GUINEA"] <- "OCEANIA"
  #
  which_continent[unique_countries=="USA"] <- "NORTH_AMERICA"
  which_continent[unique_countries=="CAYMAN ISLDS"] <- "NORTH_AMERICA"
  which_continent[unique_countries=="DOMINICAN RP"] <- "NORTH_AMERICA"
  which_continent[unique_countries=="GUADELOUPE (includes MARTINIQUE)"] <- "NORTH_AMERICA"
  which_continent[unique_countries=="NETH ANTILLES"] <- "NORTH_AMERICA"
  which_continent[unique_countries=="ST KITTS NEV (includes DOMINICA, MONTSERRAT, ST LUCA,ST VINCT, GRENADA)"] <- "NORTH_AMERICA"
  which_continent[unique_countries=="TRINIDAD-TOBAGO"] <- "NORTH_AMERICA"
  which_continent[unique_countries=="TURKS CAICOS ISL"] <- "NORTH_AMERICA"
  #
  which_continent[unique_countries=="SURINAM"] <- "SOUTH_AMERICA"
  #
  data_1980s$exp_continent <- "NA" # exporter continent
  data_1980s$imp_continent <- "NA" # importer continent
  for(i in 1:n_countries){
    country_i <- unique_countries[i]
    exp_inds <- country_i == as.character(data_1980s$expcode)
    imp_inds <- country_i == as.character(data_1980s$impcode)
    #    
    data_1980s$exp_continent[exp_inds] <- which_continent[i]
  data_1980s$imp_continent[imp_inds] <- which_continent[i]
  }
  #
  # create continent dummies
  #
  continent_exp_africa  <- ifelse(data_1980s$exp_continent=="AFRICA",1,0)
  continent_exp_asia    <- ifelse(data_1980s$exp_continent=="ASIA",1,0)
  continent_exp_europe  <- ifelse(data_1980s$exp_continent=="EUROPE",1,0)
  continent_exp_oceania <- ifelse(data_1980s$exp_continent=="OCEANIA",1,0)
  continent_exp_NA      <- ifelse(data_1980s$exp_continent=="NORTH_AMERICA",1,0)
  continent_exp_SA      <- ifelse(data_1980s$exp_continent=="SOUTH_AMERICA",1,0)
  #
  continent_imp_africa  <- ifelse(data_1980s$imp_continent=="AFRICA",1,0)
  continent_imp_asia    <- ifelse(data_1980s$imp_continent=="ASIA",1,0)
  continent_imp_europe  <- ifelse(data_1980s$imp_continent=="EUROPE",1,0)
  continent_imp_oceania <- ifelse(data_1980s$imp_continent=="OCEANIA",1,0)
  continent_imp_NA      <- ifelse(data_1980s$imp_continent=="NORTH_AMERICA",1,0)
  continent_imp_SA      <- ifelse(data_1980s$imp_continent=="SOUTH_AMERICA",1,0)
  #
  # setup regression matrix
  #
  reg_1980s <- data.frame(data_1980s$trade,
                          data_1980s$ln_trade,
                          data_1980s$year,
                          data_1980s$ln_distance,
                          data_1980s$border,
                          data_1980s$islands,
                          data_1980s$landlock,
                          data_1980s$legalsystem_same,
                          data_1980s$common_lang,
                          data_1980s$colonial,
                          data_1980s$cu,
                          data_1980s$fta,
                          data_1980s$religion_same,
                          data_1980s$impcode,
                          data_1980s$expcode,
                          data_1980s$imp_continent,
                          data_1980s$exp_continent)
  colnames(reg_1980s) <- c("trade",
                           "ln_trade",
                           "year",
                           "ln_distance",
                           "border",
                           "islands",
                           "landlock",
                           "legal",
                           "language",
                           "colonial",
                           "currency",
                           "fta",
                           "religion",
                           "impcode",
                           "expcode",
                           "imp_continent",
                           "exp_continent")
  #
  # log trade and regressors for log trade (outcome equation)
  #
  Y1 <- c(reg_1980s$ln_trade)
  X  <- cbind(1,
              -reg_1980s$ln_distance,
              -reg_1980s$border,
              -reg_1980s$islands,
              -reg_1980s$landlock,
              -reg_1980s$legal,
              -reg_1980s$language,
              -reg_1980s$colonial,
              -reg_1980s$currency,
              -reg_1980s$fta,
              -reg_1980s$religion,
              continent_exp_africa,
              continent_exp_asia,
              continent_exp_europe,
              continent_exp_NA,
              continent_exp_SA,
              continent_imp_africa,
              continent_imp_asia,
              continent_imp_europe,
              continent_imp_NA,
              continent_imp_SA)
  #
  # trade indicator and regressors for trade (selection equation)
  #
  Y2 <- matrix(c(reg_1980s$trade))
  Z  <- cbind(1,
              -reg_1980s$ln_distance,
              -reg_1980s$border,
              -reg_1980s$islands,
              -reg_1980s$landlock,
              -reg_1980s$legal,
              -reg_1980s$language,
              -reg_1980s$colonial,
              -reg_1980s$currency,
              -reg_1980s$fta,
              -reg_1980s$religion,
              continent_exp_africa,
              continent_exp_asia,
              continent_exp_europe,
              continent_exp_NA,
              continent_exp_SA,
              continent_imp_africa,
              continent_imp_asia,
              continent_imp_europe,
              continent_imp_NA,
              continent_imp_SA)
  # drop NAs from no trade in Y1
  which_na <- is.na(Y1)
  X <- X[(is.na(Y1)==FALSE),]
  Z0 <- Z[(is.na(Y1)==TRUE),]
  Z1 <- Z[(is.na(Y1)==FALSE),]
  Y1 <- matrix(Y1[(is.na(Y1)==FALSE)])
  #
  reg_data_full <- list(Y1=Y1,X=X,Z=Z,Z0=Z0,Z1=Z1)
  #
  # focus on 1986
  #
  is_1986 <- data_1980s$year==1986
  reg_1986 <- reg_1980s[is_1986,]
  #
  Y1 <- c(reg_1986$ln_trade)
  X <- cbind(1,
             -reg_1986$ln_distance,
             -reg_1986$border,
             -reg_1986$islands,
             -reg_1986$landlock,
             -reg_1986$legal,
             -reg_1986$language,
             -reg_1986$colonial,
             -reg_1986$currency,
             -reg_1986$fta,
             -reg_1986$religion,
             continent_exp_africa[is_1986],
             continent_exp_asia[is_1986],
             continent_exp_europe[is_1986],
             continent_exp_NA[is_1986],
             continent_exp_SA[is_1986],
             continent_imp_africa[is_1986],
             continent_imp_asia[is_1986],
             continent_imp_europe[is_1986],
             continent_imp_NA[is_1986],
             continent_imp_SA[is_1986])
  #
  # trade indicator and regressors for trade
  Y2 <- matrix(c(reg_1986$trade))
  Z <- cbind(1,
             -reg_1986$ln_distance,
             -reg_1986$border,
             -reg_1986$islands,
             -reg_1986$landlock,
             -reg_1986$legal,
             -reg_1986$language,
             -reg_1986$colonial,
             -reg_1986$currency,
             -reg_1986$fta,
             -reg_1986$religion,
             continent_exp_africa[is_1986],
             continent_exp_asia[is_1986],
             continent_exp_europe[is_1986],
             continent_exp_NA[is_1986],
             continent_exp_SA[is_1986],
             continent_imp_africa[is_1986],
             continent_imp_asia[is_1986],
             continent_imp_europe[is_1986],
             continent_imp_NA[is_1986],
             continent_imp_SA[is_1986])
  #
  # drop NAs from no trade in Y1
  which_na <- is.na(Y1)
  X <- X[(is.na(Y1)==FALSE),]
  Z0 <- Z[(is.na(Y1)==TRUE),]
  Y1 <- matrix(Y1[(is.na(Y1)==FALSE)])
  #
  reg_data_1986 <- list(Y1=Y1,X=X,Z0=Z0)
  #
  return(list(data_1980s=reg_data_full,data_1986=reg_data_1986))
}
#
# END