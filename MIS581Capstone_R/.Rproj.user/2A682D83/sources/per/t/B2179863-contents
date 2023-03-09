## Libraries
library(fields)

## Working Files
Assessor_Master_Res_Data <- paste(getwd(), "Working/Assessor_Master_Res_Data_2022.csv", sep = "/")
Well_Master_Data <- paste(getwd(), "Working/Master_Well_Data_2022.csv", sep = "/")

## Read in .csv file
Assessor_Data <- read.csv(Assessor_Master_Res_Data, header = TRUE)
Well_Data <- read.csv(Well_Master_Data, header = TRUE)

## FYI
# Earth Equatorial Radius in km = 6,378.388
# Earth Polar minimum Radius in km = 6,357
# Earth Mean Radius in km = 6,371

## Calculating distances between to points on Earth using Radian Lat and Long coordinates

## Convernt decimal degrees (DD) to radians
deg2rad <- function(deg) return(deg*pi/180)

### The following code for calculating distances between two sets of coordinates came from:
### Pineda-Krch, M. (2010, November 23). Great-circle distance calculations in R | R-bloggers. 
### https://www.r-bloggers.com/2010/11/great-circle-distance-calculations-in-r/#:~:text=%23%20Calculates%20the%20geodesic%20distance%20between%20two%20points


## Calculates the geodesic distance between two points specified by radian latitude/longitude using the
# Spherical Law of Cosines (slc)
gcd.slc <- function(long1, lat1, long2, lat2) {
  R <- 6371 # Earth mean radius [km]
  d <- acos(sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2) * cos(long2-long1)) * R
  return(d) # Distance in km
}

# Calculates the geodesic distance between two points specified by radian latitude/longitude using the
# Haversine formula (hf)
gcd.hf <- function(long1, lat1, long2, lat2) {
  R <- 6371 # Earth mean radius [km]
  delta.long <- (long2 - long1)
  delta.lat <- (lat2 - lat1)
  a <- sin(delta.lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta.long/2)^2
  c <- 2 * asin(min(1,sqrt(a)))
  d = R * c
  return(d) # Distance in km
}

# Calculates the geodesic distance between two points specified by radian latitude/longitude using
# Vincenty inverse formula for ellipsoids (vif)
gcd.vif <- function(long1, lat1, long2, lat2) {
  
  # WGS-84 ellipsoid parameters
  a <- 6378137         # length of major axis of the ellipsoid (radius at equator)
  b <- 6356752.314245  # length of minor axis of the ellipsoid (radius at the poles)
  f <- 1/298.257223563 # flattening of the ellipsoid
  
  L <- long2-long1 # difference in longitude
  U1 <- atan((1-f) * tan(lat1)) # reduced latitude
  U2 <- atan((1-f) * tan(lat2)) # reduced latitude
  sinU1 <- sin(U1)
  cosU1 <- cos(U1)
  sinU2 <- sin(U2)
  cosU2 <- cos(U2)
  
  cosSqAlpha <- NULL
  sinSigma <- NULL
  cosSigma <- NULL
  cos2SigmaM <- NULL
  sigma <- NULL
  
  lambda <- L
  lambdaP <- 0
  iterLimit <- 100
  while (abs(lambda-lambdaP) > 1e-12 & iterLimit>0) {
    sinLambda <- sin(lambda)
    cosLambda <- cos(lambda)
    sinSigma <- sqrt( (cosU2*sinLambda) * (cosU2*sinLambda) +
                        (cosU1*sinU2-sinU1*cosU2*cosLambda) * (cosU1*sinU2-sinU1*cosU2*cosLambda) )
    if (sinSigma==0) return(0)  # Co-incident points
    cosSigma <- sinU1*sinU2 + cosU1*cosU2*cosLambda
    sigma <- atan2(sinSigma, cosSigma)
    sinAlpha <- cosU1 * cosU2 * sinLambda / sinSigma
    cosSqAlpha <- 1 - sinAlpha*sinAlpha
    cos2SigmaM <- cosSigma - 2*sinU1*sinU2/cosSqAlpha
    if (is.na(cos2SigmaM)) cos2SigmaM <- 0  # Equatorial line: cosSqAlpha=0
    C <- f/16*cosSqAlpha*(4+f*(4-3*cosSqAlpha))
    lambdaP <- lambda
    lambda <- L + (1-C) * f * sinAlpha *
      (sigma + C*sinSigma*(cos2SigmaM+C*cosSigma*(-1+2*cos2SigmaM*cos2SigmaM)))
    iterLimit <- iterLimit - 1
  }
  if (iterLimit==0) return(NA)  # formula failed to converge
  uSq <- cosSqAlpha * (a*a - b*b) / (b*b)
  A <- 1 + uSq/16384*(4096+uSq*(-768+uSq*(320-175*uSq)))
  B <- uSq/1024 * (256+uSq*(-128+uSq*(74-47*uSq)))
  deltaSigma = B*sinSigma*(cos2SigmaM+B/4*(cosSigma*(-1+2*cos2SigmaM^2) -
                                             B/6*cos2SigmaM*(-3+4*sinSigma^2)*(-3+4*cos2SigmaM^2)))
  s <- b*A*(sigma-deltaSigma) / 1000
  
  return(s) # Distance in km
}

# Calculates the geodesic distance between two points specified by degrees (DD) latitude/longitude using
# Haversine formula (hf), Spherical Law of Cosines (slc) and Vincenty inverse formula for ellipsoids (vif)
gcd <- function(long1, lat1, long2, lat2) {
  
  # Convert degrees to radians
  long1 <- deg2rad(long1)
  lat1 <- deg2rad(lat1)
  long2 <- deg2rad(long2)
  lat2 <- deg2rad(lat2)
  
  return(list(haversine = gcd.hf(long1, lat1, long2, lat2),
              sphere = gcd.slc(long1, lat1, long2, lat2),
              vincenty = gcd.vif(long1, lat1, long2, lat2)) )
}

## R function for determining distance between to points on a sphere
## Uses coordinates in decimal degree.
# Well 1 - 06061
long1 <- -105.1659
lat1 <- 40.05646
#Well 2 - 06073
long1 <- -105.1667
lat1 <- 40.08523
#House 1
long2 <- -105.1749
lat2 <- 40.06061

dist_km <- rdist.earth(matrix(c(long1,lat1), ncol=2), matrix(c(long2,lat2), ncol=2), miles = FALSE, R = 6371)

# 1,000 ft = 0.3048 km
# 1 mile = 1.609344 km

## Stacked for loop
n1 <- nrow(Well_Data)
n2 <- nrow(Assessor_Data)

df <- data.frame(API=character(),
                 H1000=integer(), 
                 H5280=integer(), 
                 Well_PTax=numeric(), 
                 Sum_Home_PTax=numeric(),
                 stringsAsFactors = FALSE)
colnames(df) <- c("API", "H1000", "H5280", "Well_PTax", "Sum_Home_PTax")

# Test
new_row <- c(Well_Data$API_Label[1], 1, 2, Well_Data$Tax_Est[1], 2000)
df1 <- rbind(df, new_row)


for (i in 1:n1){
  tax_sum = 0
  h_count = 0
  for (j in 1:n2){
    dist_km <- rdist.earth(matrix(c(Well_Data$long1[i],Well_Data$lat1[i]), ncol=2), matrix(c(Assessor_Data$long2[j],Assessor_Data$lat2[j]), ncol=2), miles = FALSE, R = 6371)
    if (dist_km > 1.609344){
      next
    }else if (dist_km > 0.3048){
      H1K = 0
      h_count_new = 1
      h_tax = Assessor_Data$Calc_Taxes[j]
    }else if (dist_km <= 0.3048){
      H1K = 1
      h_count_new = 1
      h_tax = Assessor_Data$Calc_Taxes[j]
    }
    h_count = h_count + h_count_new
    tax_sum = tax_sum + h_tax
  }
  new_row <- c(Well_Data$API_Label[i], H1K, h_count, Well_Data$Tax_Est[i], tax_sum)
  df <- rbind(df, new_row)
}
  
View(df)


