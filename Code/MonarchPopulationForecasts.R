#####################################################################################################
#Combining results from a retrospective population model for monarch butterflies in eastern 
#North America (2004-2018) with projections of climate on the spring and summer breeding grounds 
#to forecast the population's response to future climate change

#Last updated: 21 Oct 2021
#####################################################################################################

#-----------------------------------------------------------------------------------------------------#
# Set working directory and load packages
#-----------------------------------------------------------------------------------------------------#

# setwd()

library(plyr)
library(reshape2)
library(rstan)
rstan_options(auto_write = TRUE)
rstan_options(javascript = FALSE)

# rm(list=ls())

#-----------------------------------------------------------------------------------------------------#
# Read in data for the retrospective population model (assuming these files are in a "Data" subfolder)
#-----------------------------------------------------------------------------------------------------#

#Summer monarch data
  summer <- read.csv('Data/Monarchs_summer.csv',header=TRUE,stringsAsFactors=FALSE)

#Winter monarch data and annual covariate data
  dat.y <- read.csv('Data/YearlyData.csv',header=TRUE,stringsAsFactors=FALSE)

#Covariates: county
  cov.c <- read.table('Data/Covariates_County.txt',sep='\t',header=TRUE,quote='\"',
                      colClasses=c(rep('numeric',2),rep('character',4),rep('numeric',8)))

#Covariates: county * year
  cov.cy <- read.csv('Data/Covariates_CountyYear.csv',header=TRUE)

#Covariates: county * year * week
  cov.cw <- read.csv('Data/Covariates_CountyWeek.csv',header=TRUE)

#-----------------------------------------------------------------------------------------------------#
# Retrospective population model: format summer survey data (2004-2018)
#-----------------------------------------------------------------------------------------------------#   

#Specify number of weeks, years, counties, sites, etc.
  uyears <- sort(unique(summer$yr))
  n_years <- length(uyears)
  uweeks <- sort(unique(summer$wk))
  n_weeks <- length(uweeks)
  ucounties <- cov.c$county.ind
  n_counties <- length(ucounties)
  usites <- sort(unique(summer$site.ind))
  n_sites <- length(usites)

#Add indicators for state monitoring programs (using NABA as a reference level)
  summer$ia.ind <- ifelse(summer$program=='Iowa',1,0)
  summer$il.ind <- ifelse(summer$program=='Illinois',1,0)
  summer$mi.ind <- ifelse(summer$program=='Michigan',1,0)
  summer$oh.ind <- ifelse(summer$program=='Ohio',1,0)

#Standardize estimates of %open (ie, %unforested within 2.5 or 12.5 km of a NABA or BMN survey, respectively)
  openS.m <- mean(summer$perc.open)
  openS.sd <- sd(summer$perc.open)
  summer$openS.st <- (summer$perc.open - openS.m)/openS.sd

#-----------------------------------------------------------------------------------------------------#
# Retrospective population model: standardize annual covariates
#-----------------------------------------------------------------------------------------------------#  

#Size of monarch population in late February
  feb.m <- mean(dat.y$area.feb)
  feb.sd <- sd(dat.y$area.feb)
  dat.y$feb.st <- (dat.y$area.feb - feb.m)/feb.sd
  
#GDD in spring, eastern Texas
  spGDD.m <- mean(dat.y$spGDD.east)
  spGDD.sd <- sd(dat.y$spGDD.east)
  dat.y$spGDD.st <- (dat.y$spGDD.east - spGDD.m)/spGDD.sd
  #Quadratic
  dat.y$spGDD.st2 <- dat.y$spGDD.st*dat.y$spGDD.st
  
#Precipitation in spring (FMA), eastern Texas
  spPCP.m <- mean(dat.y$spPCP.east)
  spPCP.sd <- sd(dat.y$spPCP.east)
  dat.y$spPCP.st <- (dat.y$spPCP.east - spPCP.m)/spPCP.sd
  #Quadratic
  dat.y$spPCP.st2 <- dat.y$spPCP.st*dat.y$spPCP.st

#Autumn nectar availability (NDVI, along first half of migration route)
  nectar.m <- mean(dat.y$NDVI)
  nectar.sd <- sd(dat.y$NDVI)
  dat.y$nectar.st <- (dat.y$NDVI - nectar.m)/nectar.sd  

#Dense forest cover at overwintering sites
  forest.m <- mean(dat.y$denseforest)
  forest.sd <- sd(dat.y$denseforest)
  dat.y$forest.st <- (dat.y$denseforest-forest.m)/forest.sd

#-----------------------------------------------------------------------------------------------------#  
# Retrospective population model: standardize county covariates
#-----------------------------------------------------------------------------------------------------#  

#Sort by county index
  cov.c <- cov.c[with(cov.c,order(county.ind)),]
  
#Average GDD for weeks 10-24, 2004-2018
  avgGDD.m <- mean(cov.c$avgGDD)
  avgGDD.sd <- sd(cov.c$avgGDD)
  cov.c$avgGDD.st <- (cov.c$avgGDD - avgGDD.m)/avgGDD.sd

#Average summer precipitation (AMJJA), 2004-2018
  avgPCP.m <- mean(cov.c$avgPCP)
  avgPCP.sd <- sd(cov.c$avgPCP)
  cov.c$avgPCP.st <- (cov.c$avgPCP - avgPCP.m)/avgPCP.sd  
  
#%open (ie, %unforested in each county)
  openC.m <- mean(cov.c$perc.open)
  openC.sd <- sd(cov.c$perc.open)
  cov.c$openC.st <- (cov.c$perc.open - openC.m)/openC.sd
  
#%crop
  cropC.m <- mean(cov.c$perc.crop)
  cropC.sd <- sd(cov.c$perc.crop)
  cov.c$cropC.st <- (cov.c$perc.crop - cropC.m)/cropC.sd 

#-----------------------------------------------------------------------------------------------------#  
# Retrospective population model: standardize county*year covariates
#-----------------------------------------------------------------------------------------------------#  

#Sort by county and year index
  cov.cy <- cov.cy[with(cov.cy,order(yr,county.ind)),]  
  
#Glyphosate use (proportion of corn and soy crops sprayed)
  #First, calculate county averages over 2004-2018
    gly.c <- ddply(cov.cy,.(county.ind,state.county),summarize,gly.avg=mean(glyphosate),gly.min=min(glyphosate))
    # sum(gly.c$gly.avg==0) #28 counties
    # summary(cov.c$perc.crop[cov.c$county.ind %in% c(gly.c$county.ind[gly.c$gly.avg==0])]) #all with <9% crop cover
  #Calculate minimum (non-zero) county average
    gly.c.min <- min(gly.c$gly.avg[gly.c$gly.avg!=0])
  #Impute the minimum county average in all years for those counties with glyphosate use = 0
    cov.cy$glyphosate[cov.cy$state.county %in% gly.c$state.county[gly.c$gly.avg==0]] <- gly.c.min
  
  gly.m <- mean(cov.cy$glyphosate)
  gly.sd <- sd(cov.cy$glyphosate)
  cov.cy$gly.st <- (cov.cy$glyphosate-gly.m)/gly.sd
  
#Summer precipitation (annual deviations from 2004-2018 means for AMJJA)
  diffPCP.m <- mean(cov.cy$diffPCP)
  diffPCP.sd <- sd(cov.cy$diffPCP)
  cov.cy$diffPCP.st <- (cov.cy$diffPCP - diffPCP.m)/diffPCP.sd

#-----------------------------------------------------------------------------------------------------#  
# Retrospective population model: standardize county*week covariates
#-----------------------------------------------------------------------------------------------------#  

#Sort by county, week, and year index
  cov.cw <- cov.cw[with(cov.cw,order(yr,wk,county.ind)),]   
  
#Difference between GDD and average GDD for that week and county across all years of the study
  diffGDD.m <- mean(cov.cw$diffGDD)
  diffGDD.sd <- sd(cov.cw$diffGDD)
  cov.cw$diffGDD.st <- (cov.cw$diffGDD-diffGDD.m)/diffGDD.sd  
  
#-----------------------------------------------------------------------------------------------------# 
# Retrospective population model: bundle covariates for the summer submodel
#-----------------------------------------------------------------------------------------------------#  

  cyw <- expand.grid(wk=uweeks,county.ind=ucounties,yr=uyears)
  cyw$yr.ind <- cyw$yr-2003
    #Standarize week
    wk.m <- mean(cyw$wk)
    wk.sd <- sd(cyw$wk)
    cyw$wk.st <- (cyw$wk-wk.m)/wk.sd
  cyw <- join(cyw,dat.y[,c('yr','feb.st','spGDD.st','spGDD.st2','spPCP.st','spPCP.st2')],by='yr',type='left')
  cyw <- join(cyw,cov.c[,c('county.ind','avgGDD.st','avgPCP.st','cropC.st')],by='county.ind',type='left')
  cyw <- join(cyw,cov.cy[,c('county.ind','yr','diffPCP.st','gly.st')],by=c('county.ind','yr'),type='left')
  cyw <- join(cyw,cov.cw[,c('county.ind','yr','wk','diffGDD.st')],by=c('county.ind','yr','wk'),type='left')
  
  cyw$diffavgGDD <- cyw$diffGDD.st*cyw$avgGDD.st
  cyw$diffGDD2 <- cyw$diffGDD.st*cyw$diffGDD.st
  cyw$diffavgPCP <- cyw$diffPCP.st*cyw$avgPCP.st
  cyw$diffPCP2 <- cyw$diffPCP.st*cyw$diffPCP.st
  cyw$glycrop <- cyw$gly.st*cyw$cropC.st
  
  cyw <- cyw[,c('county.ind','yr','yr.ind','wk','wk.st','feb.st','spGDD.st','spGDD.st2','spPCP.st','spPCP.st2',
                'avgGDD.st','diffGDD.st','diffGDD2','diffavgGDD',
                'avgPCP.st','diffPCP.st','diffPCP2','diffavgPCP',
                'gly.st','cropC.st','glycrop')]
  names(cyw)[6:21] <- c('feb','spGDD','spGDD2','spPCP','spPCP2',
                        'avgGDD','diffGDD','diffGDD2','diffavgGDD',
                        'avgPCP','diffPCP','diffPCP2','diffavgPCP',
                        'gly','crop','glycrop')
  
  #Need to sort in a particular way to create the model-based annual index of peak summer population size in Stan
    #Put first 5 weeks (16-20) on top
    #Then sort week-fastest, year-slowest: first 4 rows would be county=1, yr=2004, wk=21-24, next 4 rows = county=2, yr=2004, wk=21-24
    cyw15 <- cyw[cyw$wk %in% 16:20,]
    cyw69 <- cyw[cyw$wk %in% 21:24,]
    # head(cyw69); tail(cyw69)
    cyw <- rbind(cyw15,cyw69)
  
  #Finally, create an index for unique combinations of county, year, and week
  cyw$cyw.ind <- 1:nrow(cyw)
  #And attach these indices to the survey data
  summer <- join(summer,cyw[,c('county.ind','wk','yr','cyw.ind')],by=c('county.ind','wk','yr'),type='left')

#-----------------------------------------------------------------------------------------------------# 
# Retrospective population model: county weights based on area unforested land
#-----------------------------------------------------------------------------------------------------# 

  area.open <- as.vector(cov.c$area.land.sqmi*cov.c$perc.open/100)
  weights.open <- area.open/sum(area.open)

#-----------------------------------------------------------------------------------------------------# 
# Retrospective population model: package data, initial values, and parameters for Stan
#-----------------------------------------------------------------------------------------------------# 

#Survey-level covariates in matrix
  X_survey <- as.matrix(summer[,c('ia.ind','il.ind','mi.ind','oh.ind','openS.st')]) 
  
#County-year-week covariates 
  X_county <- as.matrix(cyw[,c('feb','spGDD','spGDD2','spPCP','spPCP2',
                               'avgGDD','diffGDD','diffGDD2','diffavgGDD',
                               'avgPCP','diffPCP','diffPCP2','diffavgPCP',
                               'gly','crop','glycrop')]) 

#Covariates in winter submodel  
  X_winter <- as.matrix(dat.y[,c('forest.st','nectar.st')])

#Bundle data
  standata <- list(n_years=n_years,
                   n_counties=n_counties,
                   n_cyw=nrow(cyw),
                   n_sites=n_sites,
                   n_surveys=nrow(summer),
                   year_id=cyw$yr.ind,
                   county_id=cyw$county.ind,
                   site_id=summer$site.ind,
                   cyw_id=summer$cyw.ind,
                   n_cov_alpha=ncol(X_county),
                   n_cov_beta=ncol(X_survey),
                   n_cov_gamma=ncol(X_winter),
                   X_county=X_county,
                   week_st=cyw$wk.st,
                   X_survey=X_survey,
                   effort=summer$duration,
                   X_winter=X_winter,
                   y_count=summer$monarch,
                   area=dat.y$area.dec,
                   weights=weights.open,
                   ind1=seq(1,n_counties*n_years*4,by=4),
                   ind2=seq(1,n_counties*n_years,by=n_counties),
                   start_peak=which(cyw$wk %in% 21:24)[1],
                   n_peak=sum(cyw$wk %in% 21:24),
                   n_cy=n_counties*n_years)

#MCMC parameters
  ni <- 4000       # No. iterations (including warmup)
  nb <- 3000       # No. burn-in iterations to discard (ie, warmup)
  nt <- 1          # Thin rate
  nc <- 3          # No. chains
  
#Initial values
  set.seed(126)
  inits <- lapply(1:nc, function(i)
           list(alpha0=runif(1,2,3),
                alphaFE=runif(ncol(X_county),-0.5,0.5),
                alphaRE_week=runif(1,-0.5,0.5),
                alphaRE_week2=runif(1,-0.5,0.5),
                betaFE=runif(ncol(X_survey),-0.5,0.5),
                gamma0=runif(1,0,2),
                gamma_sum=runif(1,-0.5,0.5),
                gammaFE=runif(ncol(X_winter),-0.5,0.5),
                sd_county=runif(1,0,1),
                sd_week=runif(1,0,1),
                sd_week2=runif(1,0,1),
                sd_site=runif(1,0,1),
                sd_year=runif(1,0,1),
                r_count=runif(1,0,2),
                shape=runif(1,0,2)))

#Parameters to monitor
  params <- c('alpha0','alphaFE','alphaRE_week','alphaRE_week2','betaFE',
              'gamma0','gamma_sum','gammaFE','sd_county','sd_week','sd_week2',
              'sd_site','sd_year','r_count','shape','pred_orig_log','pred_orig_exp','pred_sum','mu_win',
              'randweek','randweek2','randyear','randcounty')

#-----------------------------------------------------------------------------------------------------# 
# Retrospective population model: Run model in Stan using rstan package and save posterior samples
#-----------------------------------------------------------------------------------------------------# 

  out <- stan('FACM_2004-2018.stan',
              control=list(adapt_delta=0.8),
              data=standata, init=inits, pars=params,
              chains=nc, iter=ni, warmup=nb, thin=nt,
              seed=1,cores=3,open_progress=FALSE)

  #If run previously, load workspace with stanfit object
  # load('...Rdata')

  posterior <- as.matrix(out)
  #Extract 1000 samples from posterior distribution (thin by 3)
  iter <- seq(1,3000,by=3)
  posterior <- posterior[iter,]
  
  #Rename matrix of covariate values used in retrospective model
  cyw.obs <- cyw
  
  #Remove objects from workspace that aren't needed anymore
  rm(standata)
  rm(out)
  rm(X_survey)
  rm(X_county)
  rm(X_winter)

#-----------------------------------------------------------------------------------------------------# 
# Climate projections: import and format climate data
#-----------------------------------------------------------------------------------------------------# 

#Import spring data 
  springclim <- read.csv('Data/SpringClimateProjections.csv',header=TRUE)
  names(springclim)[1:3] <- c('scenarioname','yr','modelname')

#Import summer data
  summerclim1 <- read.csv('Data/SummerClimateProjections1.csv',header=TRUE)
  summerclim2 <- read.csv('Data/SummerClimateProjections2.csv',header=TRUE)
  summerclim3 <- read.csv('Data/SummerClimateProjections3.csv',header=TRUE)
  summerclim <- rbind(summerclim1,summerclim2,summerclim3)
  names(summerclim)[1:3] <- c('scenarioname','yr','modelname')
  rm(list=c('summerclim1','summerclim2','summerclim3'))
  
#Time periods
  periods <- 1:3
  
#Years
  years1 <- 2023:2043
  years2 <- 2050:2070
  years3 <- 2080:2100
  
#Emissions scenarios: 
  #1 = low, SSP1-2.6 (SSP126); 2 = moderate, SSP2-4.5 (SSP245); 3 = high, SSP3-7.0 (SSP370); 4 = worst-case, SSP5-8.5 (SSP585) 
  sc <- data.frame(scenario=1:4,scenarioname=c('SSP126','SSP245','SSP370','SSP585'))
  
#GCMs
  mod <- data.frame(model=1:5,modelname=sort(unique(springclim$modelname)))

#Adding time period, model number, and scenario number to spring and summer data  
  springclim$period <- ifelse(springclim$yr %in% years1,1,ifelse(springclim$yr %in% years2,2,3))
  summerclim$period <- ifelse(summerclim$yr %in% years1,1,ifelse(summerclim$yr %in% years2,2,3))
  springclim$model <- mod$model[match(springclim$modelname,mod$modelname)]
  summerclim$model <- mod$model[match(summerclim$modelname,mod$modelname)]
  springclim$scenario <- sc$scenario[match(springclim$scenarioname,sc$scenarioname)]
  summerclim$scenario <- sc$scenario[match(summerclim$scenarioname,sc$scenarioname)]
  
#Adding county index (1:545) to summer data 
  summerclim$county.ind <- cov.c$county.ind[match(summerclim$state.county,cov.c$state.county)]
  
#Sort dataframes
  springclim <- springclim[with(springclim,order(scenario,model,yr)),]
  summerclim <- summerclim[with(summerclim,order(scenario,model,yr,county.ind)),]
  
#-----------------------------------------------------------------------------------------------------# 
# Climate projections: create climate covariates and standardize
#-----------------------------------------------------------------------------------------------------# 
#Standardizing all covariates using the mean/SD of 2004-2018 data  

#Spring GDD (spGDD)
  springclim$spGDD.st <- (springclim$spGDD - spGDD.m)/spGDD.sd
  #Quadratic
  springclim$spGDD.st2 <- springclim$spGDD.st*springclim$spGDD.st
  
#Spring precipitation (spPCP)
  springclim$spPCP <- rowSums(springclim[,c('pcp.Feb','pcp.Mar','pcp.Apr')])
  springclim$spPCP.st <- (springclim$spPCP - spPCP.m)/spPCP.sd
  #Quadratic
  springclim$spPCP.st2 <- springclim$spPCP.st*springclim$spPCP.st

#Summer GDD (avgGDD and diffGDD)
  avgGDD <- ddply(summerclim,.(scenario,model,period,county.ind,state.county),summarize,avgGDD=mean(gdd.wk10.24))
  avgGDD$avgGDD.st <- (avgGDD$avgGDD - avgGDD.m)/avgGDD.sd
  
  #Calculate difference from average for that time period, for each county, year, and week
  GDDlong <- melt(summerclim,id.vars=c('scenario','model','period','yr','county.ind','state.county'),
                  measure.vars=paste('gdd.wk10',21:24,sep='.'),variable.name='wk',value.name='GDD')
  GDDlong$wk <- as.numeric(substr(GDDlong$wk,10,11))
  wkGDDavg <- ddply(GDDlong,.(scenario,model,period,county.ind,state.county,wk),summarize,
                    GDD.mn=mean(GDD))  
  diffGDD <- join(GDDlong,wkGDDavg,by=c('scenario','model','period','wk','county.ind','state.county'),type='left')
  diffGDD$diffGDD <- diffGDD$GDD - diffGDD$GDD.mn
  diffGDD$diffGDD.st <- (diffGDD$diffGDD - diffGDD.m)/diffGDD.sd
  #Quadratic
  diffGDD$diffGDD.st2 <- diffGDD$diffGDD.st*diffGDD$diffGDD.st
  
#Summer PCP (avgPCP and diffGDD)
  summerclim$suPCP <- rowSums(summerclim[,c('pcp.Apr','pcp.May','pcp.Jun','pcp.Jul','pcp.Aug')])
  avgPCP <- ddply(summerclim,.(scenario,model,period,county.ind,state.county),summarize,avgPCP=mean(suPCP))
  avgPCP$avgPCP.st <- (avgPCP$avgPCP - avgPCP.m)/avgPCP.sd
  
  #Calculate difference from average for that time period, for each county and year
  diffPCP <- join(summerclim[,c('scenario','model','period','yr','county.ind','state.county','suPCP')],
                  avgPCP[,c('scenario','model','period','county.ind','state.county','avgPCP')],
                  by=c('scenario','model','period','county.ind','state.county'),type='left')
  diffPCP$diffPCP <- diffPCP$suPCP - diffPCP$avgPCP
  diffPCP$diffPCP.st <- (diffPCP$diffPCP - diffPCP.m)/diffPCP.sd
  #Quadratic
  diffPCP$diffPCP.st2 <- diffPCP$diffPCP.st*diffPCP$diffPCP.st

#-----------------------------------------------------------------------------------------------------# 
# Population forecasts: format and standardize other (non-climate) covariates 
#----------------------------------------------------------------------------------------------------# 
#February population size (feb) assuming value equal to 2004-2018 mean; standardized = 0
#Forest cover (forest) assuming value equal to 2004-2018 mean; standardized = 0
#Nectar availability (nectar) assuming value equal to 2004-2018 mean; standardized = 0
#Crop cover: using 2004-2018 county-level values, which were time-invariant (cov.c$cropC.st)

#Glyphosate: using 2004-2018 county-level means (gly.c)
  #Impute the minimum non-zero value (0.62) for counties with no data
  gly.c$gly.avg[gly.c$gly.avg==0] <- gly.c.min
  gly.c$gly.st <- (gly.c$gly.avg - gly.m)/gly.sd
  
#-----------------------------------------------------------------------------------------------------# 
# Population forecasts: create vectors/matrices with all posterior samples or the median posterior value from retrospective model
#----------------------------------------------------------------------------------------------------# 

  ucounties <- cov.c$county.ind
  n_counties <- length(ucounties) 
  uweeks <- sort(unique(diffGDD$wk))
  n_weeks <- length(uweeks)
  n_years <- length(years1)
  n_iter <- nrow(posterior)
  
#Matrix with posterior samples for fixed effects in summer submodel
    alpha <- as.matrix(posterior[,grep('alpha0|alphaFE',colnames(posterior))])
    #Remove samples for regression coefficient associated with February population size
    alpha <- alpha[,colnames(alpha)!='alphaFE[1]']
  #Median values from posterior distribution
    alpha.md <- apply(alpha,2,median)

#Vectors with posterior samples for week effects
    alpha.wk <- posterior[,'alphaRE_week']
    alpha.wk2 <- posterior[,'alphaRE_week2']
    sd.wk <- posterior[,'sd_week']
    sd.wk2 <- posterior[,'sd_week2']
  #Median values of week effects
    alpha.wk.md <- median(alpha.wk)
    alpha.wk2.md <- median(alpha.wk2)

#Matrix with posterior samples for county random effects          
    RE.counties <- as.matrix(posterior[,grep('randcounty',colnames(posterior))])
    RE.counties <- t(RE.counties)
    RE.county1yr <- RE.counties[rep(1:nrow(RE.counties),each=n_weeks),]
  #Median values from posterior samples
    RE.county1yr.md <- apply(RE.county1yr,1,median)
    RE.countyallyr.md <- rep(RE.county1yr.md,n_years)
    
#Matrix with posterior samples from intercept in winter submodel
    gamma0 <- as.matrix(posterior[,'gamma0'])
    WinInt <- as.matrix(rep(1,n_years),nrow=n_years)
    Wgamma0 <- WinInt %*% t(gamma0) 
  #Median of those values
    Wgamma0.md <- apply(Wgamma0,1,median)

#Matrix with posterior samples from summer effect in winter submodel
    gamma.sum <- as.matrix(posterior[,'gamma_sum'])
  #Median value
    gamma.sum.md <- apply(gamma.sum,2,median)
    
#Vector with posterior samples from SD parameter associated with random yearly effects in winter submodel
    sd.yr <- posterior[,'sd_year']
    
#Matrix with posterior samples from shape parameter (for gamma distribution of winter areas)
    shape <- as.matrix(posterior[,'shape'])

#-----------------------------------------------------------------------------------------------------# 
# Population forecasts: forecasting summer counts and winter population size for each scenario, GCM, year
# Accounting for parameter uncertainty, climate uncertainty, environmental stochasticity
#-----------------------------------------------------------------------------------------------------# 

#Forecasts will be compiled into lists:  
  win.list <- list()
  sum.list <- list()
  
  i <- 0
  set.seed(44)
  
#Loop through emissions scenarios, GCMs, and time periods:
  
  for(ss in 1:nrow(sc)){ #scenario

    for(mm in 1:nrow(mod)){ #GCM

      for(pp in 1:3){ #time period
        uyears <- get(paste0('years',pp))
        i <- i+1
        
        #Combine ecological covariates in summer submodel in dataframe
          cyw <- expand.grid(wk=uweeks,county.ind=ucounties,yr=uyears)
          cyw$yr.ind <- cyw$yr - min(cyw$yr) + 1
            #Standardize week (want same values as those used in 2004-2018 analysis)
            cyw$wk.st <- (cyw$wk-wk.m)/wk.sd
          cyw <- join(cyw,springclim[springclim$scenario==ss & springclim$model==mm & springclim$period==pp,
                                     c('yr','spGDD.st','spGDD.st2','spPCP.st','spPCP.st2')],by='yr',type='left')
          cyw <- join_all(list(cyw,cov.c[,c('county.ind','cropC.st')],gly.c[,c('county.ind','gly.st')],
                               avgGDD[avgGDD$scenario==ss & avgGDD$model==mm & avgGDD$period==pp,c('county.ind','avgGDD.st')],
                               avgPCP[avgPCP$scenario==ss & avgPCP$model==mm & avgPCP$period==pp,c('county.ind','avgPCP.st')]),
                          by='county.ind',type='left')
          cyw <- join(cyw,diffPCP[diffPCP$scenario==ss & diffPCP$model==mm & diffPCP$period==pp,c('county.ind','yr','diffPCP.st','diffPCP.st2')],
                      by=c('county.ind','yr'),type='left')
          cyw <- join(cyw,diffGDD[diffGDD$scenario==ss & diffGDD$model==mm & diffGDD$period==pp,c('county.ind','yr','wk','diffGDD.st','diffGDD.st2')],
                      by=c('county.ind','yr','wk'),type='left')
          cyw$diffavgGDD <- cyw$diffGDD.st*cyw$avgGDD.st
          cyw$diffavgPCP <- cyw$diffPCP.st*cyw$avgPCP.st
          cyw$glycrop <- cyw$gly.st*cyw$cropC.st
          cyw <- cyw[,c('county.ind','yr','yr.ind','wk','wk.st','spGDD.st','spGDD.st2','spPCP.st','spPCP.st2',
                        'avgGDD.st','diffGDD.st','diffGDD.st2','diffavgGDD',
                        'avgPCP.st','diffPCP.st','diffPCP.st2','diffavgPCP',
                        'gly.st','cropC.st','glycrop')]
          names(cyw)[6:ncol(cyw)] <- c('spGDD','spGDD2','spPCP','spPCP2',
                                       'avgGDD','diffGDD','diffGDD2','diffavgGDD',
                                       'avgPCP','diffPCP','diffPCP2','diffavgPCP',
                                       'gly','crop','glycrop')
         
        #Calculate summer expected counts from posterior samples
          #Fixed effects
          Xcounty <- as.matrix(cbind(rep(1,nrow(cyw)),
                                     cyw[,c('spGDD','spGDD2','spPCP','spPCP2',
                                            'avgGDD','diffGDD','diffGDD2','diffavgGDD',
                                            'avgPCP','diffPCP','diffPCP2','diffavgPCP',
                                            'gly','crop','glycrop')]))
          Xalpha <- Xcounty %*% t(alpha)

          #Effect of week
          Xweeks <- matrix(c(cyw$wk.st,cyw$wk.st*cyw$wk.st),ncol=2,nrow=nrow(cyw),byrow=FALSE)  
          alphaRE.wk <- alphaRE.wk2 <- matrix(NA,nrow=n_years,ncol=n_iter)
          for(s in 1:n_years){
            alphaRE.wk[s,] <- alpha.wk + rnorm(n_iter,0,sd.wk)
            alphaRE.wk2[s,] <- alpha.wk2 + rnorm(n_iter,0,sd.wk2)
          } #s
          Xalpha.wk <- Xalpha.wk2 <- matrix(NA,nrow=nrow(Xweeks),ncol=n_iter)
          for(s in 1:nrow(cyw)){
            Xalpha.wk[s,] <- Xweeks[s,1]*alphaRE.wk[cyw$yr.ind[s],] 
            Xalpha.wk2[s,] <- Xweeks[s,2]*alphaRE.wk2[cyw$yr.ind[s],] 
          } #s
          
          #Random county effects
          eps.county <- RE.county1yr[rep(1:nrow(RE.county1yr),times=n_years),]

          #Combine components in summer submodel
          mu.county.log <- Xalpha + Xalpha.wk + Xalpha.wk2 + eps.county
    
          #Calculate means for each county-year combination
          split.muc <- lapply(split(mu.county.log,rep(1:(n_counties*n_years),each=n_weeks)),matrix,ncol=ncol(mu.county.log))
          mu.cy.log <- t(sapply(split.muc,colMeans))
          
        #Add matrix of forecasted counts (in each county and year) to sum.list
          mu.cy <- exp(mu.cy.log)
          params.cy <- matrix(as.numeric(c(ss,mm,pp)),nrow=n_years*n_counties,ncol=3,byrow=TRUE)
          params.cy <- cbind(params.cy,unique(cyw[,c('county.ind','yr')]))
          sum.list[[i]] <- cbind(params.cy,mu.cy)
          colnames(sum.list[[i]]) <- c('scenario','model','period','county.ind','yr',paste0('i',1:n_iter))
      
        #Calculate weighted mean of expected counts in each year (mean across counties)
          weights.mat <- matrix(weights.open,nrow=nrow(mu.cy.log),ncol=ncol(mu.cy.log),byrow=FALSE)  
          mu.cy.weights <- mu.cy.log*weights.mat
          split.mucy <- lapply(split(mu.cy.weights,rep(1:n_years,each=n_counties)),matrix,ncol=ncol(mu.cy.weights))
          mu.y.log <- t(sapply(split.mucy,colSums))
          
          #Standardize annual population index by fixed mean, SD
          mu.y.z <- (mu.y.log - 1.12)/0.59
          
        #Calculate expected area occupied from posterior samples    
          Wsum <- matrix(NA,nrow=n_years,ncol=n_iter)
          Wyr <- matrix(NA,nrow=n_years,ncol=n_iter)          
          
          #Summer effect 
            for(ii in 1:n_iter){
              Wsum[,ii] <- mu.y.z[,ii]*gamma.sum[ii,1]
          
          #Random yearly effects
              Wyr[,ii] <- rnorm(n_years,0,sd.yr[ii])
            } #ii

          #Combine components and exponentiate
            mu.win <- exp(Wgamma0 + Wsum + Wyr)
            
          #Calculate rate parameter
            rate <- matrix(NA,nrow=n_years,ncol=n_iter)
            pred.win <- matrix(NA,nrow=n_years,ncol=n_iter)
            for(ii in 1:n_iter){
              rate[,ii] <- shape[ii,1]*(1/mu.win[,ii])
            
          #Generate value from gamma distribution with specified shape, rate parameters
              pred.win[,ii] <- rgamma(n_years,shape=shape[ii,1],rate=rate[,ii])   
            } #ii           
        
        #Add matrix of annual forecasted values of the area occupied to win.list  
          params.y <- matrix(as.numeric(c(ss,mm,pp)),nrow=n_years,ncol=3,byrow=TRUE)
          params.y <- cbind(params.y,uyears)
          win.list[[i]] <- cbind(params.y,pred.win)
          colnames(win.list[[i]]) <- c('scenario','model','period','yr',paste0('i',1:n_iter)) 
          
        #See progress
          print(i)

      } #pp
    } #mm
  } #ss
  
  preds.sum <- do.call(rbind,sum.list)

#-----------------------------------------------------------------------------------------------------# 
# Summarize forecasted monarch counts on the summer breeding grounds during each future time period
#-----------------------------------------------------------------------------------------------------#  

#-----------------------------------------------------------------------------#
#For comparison, first calculate the expected means in each county and year between 2004-2018
  cyw.obs <- cyw.obs[cyw.obs$wk %in% 21:24,]
  Xcounty <- as.matrix(cbind(rep(1,nrow(cyw.obs)),                                     
                             cyw.obs[,c('feb','spGDD','spGDD2','spPCP','spPCP2',
                                        'avgGDD','diffGDD','diffGDD2','diffavgGDD',
                                        'avgPCP','diffPCP','diffPCP2','diffavgPCP',
                                        'gly','crop','glycrop')]))
  alpha.full <- as.matrix(posterior[,grep('alpha0|alphaFE',colnames(posterior))])
  Xalpha <- Xcounty %*% t(alpha.full)
  
  Xweeks <- matrix(c(cyw.obs$wk.st,cyw.obs$wk.st*cyw.obs$wk.st),ncol=2,byrow=FALSE)
  randweekboth <- posterior[,grep('randweek',colnames(posterior))]
  randweek <- randweekboth[,1:15]
  randweek2 <- randweekboth[,16:30]
  alphaRE.wk <- alphaRE.wk2 <- matrix(NA,nrow=15,ncol=n_iter) 
  for(s in 1:15){
    alphaRE.wk[s,] <- alpha.wk + randweek[,s]
    alphaRE.wk2[s,] <- alpha.wk2 + randweek2[,s]
  }
  Xalpha.wk <- Xalpha.wk2 <- matrix(NA,nrow=nrow(Xweeks),ncol=n_iter)  
  for(s in 1:nrow(cyw.obs)){
    Xalpha.wk[s,] <- Xweeks[s,1]*alphaRE.wk[cyw.obs$yr.ind[s],]
    Xalpha.wk2[s,] <- Xweeks[s,2]*alphaRE.wk2[cyw.obs$yr.ind[s],]
  }
  
  eps.county <- RE.county1yr[rep(1:nrow(RE.county1yr),times=15),]
  
  #Combine components
  obs.log <- Xalpha + Xalpha.wk + Xalpha.wk2 + eps.county

  #Calculate means (across weeks) for each county-year-iteration combination
  split.cy <- lapply(split(obs.log,rep(1:(n_counties*15),each=4)),matrix,ncol=ncol(obs.log))
  cy.log <- t(sapply(split.cy,colMeans))
  cy.real <- as.data.frame(exp(cy.log))
  
  #For each county, calculate the 15-year mean (for each iteration)
  cy.obs <- data.frame(unique(cyw.obs[,c('yr','county.ind')]))
  c.obs.list <- split(cy.real,f=cy.obs$county.ind)
  c.obs.byperiod <- lapply(c.obs.list,colMeans)
  c.obs.byperiod <- do.call(rbind,c.obs.byperiod)
  
  #Then summarize across iterations for each county
  c.obs.summary <- data.frame(county.ind=unique(cy.obs$county.ind))
  c.obs.summary$obs.mn <- apply(c.obs.byperiod,1,mean) 
  c.obs.summary$obs.sd <- apply(c.obs.byperiod,1,sd)
  c.obs.summary$obs.q0.05 <- apply(c.obs.byperiod,1,quantile,0.05) 
  c.obs.summary$obs.q0.5 <- apply(c.obs.byperiod,1,quantile,0.5) 
  c.obs.summary$obs.q0.95 <- apply(c.obs.byperiod,1,quantile,0.95) 
  c.obs.summary$state.county <- cov.c$state.county[match(c.obs.summary$county.ind,cov.c$county.ind)]
  
  #Summary of median counts
  summary(c.obs.summary$obs.q0.5)
  
#-----------------------------------------------------------------------------#
#Summarizing forecasts for a low emissions scenario, SSP126 (across GCMs)
  preds.sum1 <- preds.sum[preds.sum$scenario==1,]
  index.cpm <- paste(preds.sum1$county.ind,preds.sum1$period,preds.sum1$model,sep='-')
  preds.list.cpm1 <- split(preds.sum1,f=index.cpm)
  #For each county and GCM, calculate the 21-year mean (for each iteration)
  preds.bycpm1 <- lapply(preds.list.cpm1,colMeans)
  preds.bycpm1a <- do.call(rbind,preds.bycpm1)
  preds.bycpm1a <- as.data.frame(preds.bycpm1a)
  
  #Summarize period means across 5000 iteration-GCM combinations
  index.cp <- paste(preds.bycpm1a$county.ind,preds.bycpm1a$period,sep='-')
  preds.bycp1 <- split(preds.bycpm1a,f=index.cp) 
  preds.summary.cp1 <- unique(preds.bycpm1a[,c('period','county.ind')])
  preds.summary.cp1$mn <- sapply(preds.bycp1,function(x) mean(as.matrix(x[,6:(n_iter+5)])))
  preds.summary.cp1$sd <- sapply(preds.bycp1,function(x) sd(as.matrix(x[,6:(n_iter+5)])))
  preds.summary.cp1$q0.05 <- sapply(preds.bycp1,function(x) quantile(as.matrix(x[,6:(n_iter+5)]),0.05))
  preds.summary.cp1$q0.5 <- sapply(preds.bycp1,function(x) quantile(as.matrix(x[,6:(n_iter+5)]),0.5))
  preds.summary.cp1$q0.95 <- sapply(preds.bycp1,function(x) quantile(as.matrix(x[,6:(n_iter+5)]),0.95))
  
  #Attach 2004-2018 mean/median, and county covariates
  preds.summary.cp1 <- join(preds.summary.cp1,c.obs.summary[,c('county.ind','obs.mn','obs.q0.5')],
                            by='county.ind',type='left')
  preds.summary.cp1 <- join(preds.summary.cp1,cov.c[,c('county.ind','state.county')],
                            by='county.ind',type='left')
  
  #Calculate absolute and percent difference from 2004-2018 means
  #Using median of distributions of means
    preds.summary.cp1$diff <- preds.summary.cp1$q0.5 - preds.summary.cp1$obs.q0.5
    preds.summary.cp1$perc.diff <- preds.summary.cp1$diff/preds.summary.cp1$obs.q0.5*100  
    
#-----------------------------------------------------------------------------#
#Summarizing forecasts for a high emissions scenario, SSP585 (across GCMs)
  preds.sum4 <- preds.sum[preds.sum$scenario==4,]
  index.cpm <- paste(preds.sum4$county.ind,preds.sum4$period,preds.sum4$model,sep='-')
  preds.list.cpm4 <- split(preds.sum4,f=index.cpm)
  #For each county and GCM, calculate the 21-year mean (for each iteration)
  preds.bycpm4 <- lapply(preds.list.cpm4,colMeans)
  preds.bycpm4a <- do.call(rbind,preds.bycpm4)
  preds.bycpm4a <- as.data.frame(preds.bycpm4a)
  
  #Summarize period means across 5000 iteration-GCM combinations
  index.cp <- paste(preds.bycpm4a$county.ind,preds.bycpm4a$period,sep='-')
  preds.bycp4 <- split(preds.bycpm4a,f=index.cp) 
  preds.summary.cp4 <- unique(preds.bycpm4a[,c('period','county.ind')])
  preds.summary.cp4$mn <- sapply(preds.bycp4,function(x) mean(as.matrix(x[,6:(n_iter+5)])))
  preds.summary.cp4$sd <- sapply(preds.bycp4,function(x) sd(as.matrix(x[,6:(n_iter+5)])))
  preds.summary.cp4$q0.05 <- sapply(preds.bycp4,function(x) quantile(as.matrix(x[,6:(n_iter+5)]),0.05))
  preds.summary.cp4$q0.5 <- sapply(preds.bycp4,function(x) quantile(as.matrix(x[,6:(n_iter+5)]),0.5))
  preds.summary.cp4$q0.95 <- sapply(preds.bycp4,function(x) quantile(as.matrix(x[,6:(n_iter+5)]),0.95))
  
  #Attach 2004-2018 mean/median, and county covariates
  preds.summary.cp4 <- join(preds.summary.cp4,c.obs.summary[,c('county.ind','obs.mn','obs.q0.5')],
                            by='county.ind',type='left')
  preds.summary.cp4 <- join(preds.summary.cp4,cov.c[,c('county.ind','state.county')],
                            by='county.ind',type='left')
  
  #Calculate absolute and percent difference from 2004-2018 means
  #Using median of distributions of means
    preds.summary.cp4$diff <- preds.summary.cp4$q0.5 - preds.summary.cp4$obs.q0.5
    preds.summary.cp4$perc.diff <- preds.summary.cp4$diff/preds.summary.cp4$obs.q0.5*100  
    
#-----------------------------------------------------------------------------------------------------# 
# Summarize forecasted area occupied on the overwintering grounds during each future time period
#-----------------------------------------------------------------------------------------------------#  

#-----------------------------------------------------------------------------#
#For comparison, summarize observed data, 2004-2018
  obsmn.0418 <- mean(dat.y$area.dec)
  #Standard error and 90% CI
  se.obs0418 <- sd(dat.y$area.dec)/sqrt(15)
  ci90.obs0418 <- obsmn.0418 + c(-1,1)*1.761*se.obs0418  #(t-stat for 90% CI, 14 df is 1.761)
  
#-----------------------------------------------------------------------------#
#Summarize forecasts for each GCM-emission scenario combination
  #First, average forecasts over years for each iteration-period-GCM-scenario
  preds.byperiod <- lapply(win.list,colMeans)
  #Then summarize period means (mean, sd, quants) over iterations for each GCM-scenario combination
  preds.byperiod2 <- do.call(rbind,preds.byperiod)
  preds.summary <- as.data.frame(preds.byperiod2[,1:3])
  preds.summary$mn <- apply(preds.byperiod2[,5:(n_iter+4)],1,mean) 
  preds.summary$sd <- apply(preds.byperiod2[,5:(n_iter+4)],1,sd) 
  preds.summary$q0.05 <- apply(preds.byperiod2[,5:(n_iter+4)],1,quantile,0.05) 
  preds.summary$q0.25 <- apply(preds.byperiod2[,5:(n_iter+4)],1,quantile,0.25) 
  preds.summary$q0.5 <- apply(preds.byperiod2[,5:(n_iter+4)],1,quantile,0.5) 
  preds.summary$q0.75 <- apply(preds.byperiod2[,5:(n_iter+4)],1,quantile,0.75) 
  preds.summary$q0.95 <- apply(preds.byperiod2[,5:(n_iter+4)],1,quantile,0.95) 

#-----------------------------------------------------------------------------#  
#Summarize forecasts for each emission scenario (across GCMs)
  #Summarize period means (mean, sd, quants) across (5000) iteration-GCM combinations for each scenario
  preds.byperiod.df <- as.data.frame(preds.byperiod2)
  index.s <- paste(preds.byperiod.df$scenario,preds.byperiod.df$period,sep='-')
  preds.list.s <- split(preds.byperiod.df,f=index.s)
  preds.summary.s <- unique(preds.byperiod.df[,c('scenario','period')])
  preds.summary.s$mn <- sapply(preds.list.s,function(x) mean(as.matrix(x[,5:(n_iter+4)])))
  preds.summary.s$sd <- sapply(preds.list.s,function(x) sd(as.matrix(x[,5:(n_iter+4)])))
  preds.summary.s$q0.05 <- sapply(preds.list.s,function(x) quantile(as.matrix(x[,5:(n_iter+4)]),0.05))
  preds.summary.s$q0.25 <- sapply(preds.list.s,function(x) quantile(as.matrix(x[,5:(n_iter+4)]),0.25))
  preds.summary.s$q0.5 <- sapply(preds.list.s,function(x) quantile(as.matrix(x[,5:(n_iter+4)]),0.5))
  preds.summary.s$q0.75 <- sapply(preds.list.s,function(x) quantile(as.matrix(x[,5:(n_iter+4)]),0.75))
  preds.summary.s$q0.95 <- sapply(preds.list.s,function(x) quantile(as.matrix(x[,5:(n_iter+4)]),0.95))
  pss <- preds.summary.s

#-----------------------------------------------------------------------------#  
#Plot with forecasts (medians, 50 & 90% CIs), averaged across GCMs for each emission scenario 
#(Figs. 2a-c in paper)

  figcol4 <- c('dodgerblue4','darkseagreen4','lightsalmon2','tomato4')  
  
  par(mfrow=c(1,3),mar=c(0,0,0.3,0),oma=c(1.4,3.5,0.2,0.5),cex=0.9)
  plot(0,0,type='n',xlim=c(0,max(pss$scenario)+1),ylim=c(0,9.6),axes=FALSE,xlab='',ylab='')
    polygon(x=c(par('usr')[1],par('usr')[2],par('usr')[2],par('usr')[1],par('usr')[1]),
            y=c(ci90.obs0418[1],ci90.obs0418[1],ci90.obs0418[2],ci90.obs0418[2],ci90.obs0418[1]),
            col='gray95',border=NA)  
    arrows(x0=par('usr')[1],x1=par('usr')[2],y0=obsmn.0418,y1=obsmn.0418,length=0,col='gray40',lty=2)
    axis(1,at=par('usr')[1:2],tcl=0,labels=FALSE)
    axis(2,at=par('usr')[3:4],tcl=0,labels=FALSE)
    axis(2,at=seq(0,9,by=2),tcl=-0.25,labels=seq(0,9,by=2),las=1,mgp=c(1.5,0.5,0))
    axis(3,at=par('usr')[1:2],tcl=0,labels=FALSE)
    axis(4,at=par('usr')[3:4],tcl=0,labels=FALSE)
    arrows(x0=pss$scenario[pss$period==1],x1=pss$scenario[pss$period==1],y0=pss$q0.05[pss$period==1],y1=pss$q0.95[pss$period==1],
           length=0,col=figcol4,lwd=1.2)
    arrows(x0=pss$scenario[pss$period==1],x1=pss$scenario[pss$period==1],y0=pss$q0.25[pss$period==1],y1=pss$q0.75[pss$period==1],
           length=0,col=figcol4,lwd=2.5)
    arrows(x0=par('usr')[1],x1=par('usr')[2],y0=obsmn.0418,y1=obsmn.0418,length=0,col='gray40',lty=2)
    points(q0.5~scenario,data=pss[pss$period==1,],pch=21,col=figcol4,bg=figcol4)
    legend(x=-0.2,y=10.1,c('Low','Mod-Low','Mod-High','High'),col=figcol4,lwd=2.5,bty='n',cex=0.8,x.intersp=0.7)
    mtext('2023-2043',side=1,line=0.2,cex=0.9)
  plot(0,0,type='n',xlim=c(0,max(pss$scenario)+1),ylim=c(0,9.6),axes=FALSE,xlab='',ylab='')
    polygon(x=c(par('usr')[1],par('usr')[2],par('usr')[2],par('usr')[1],par('usr')[1]),
            y=c(ci90.obs0418[1],ci90.obs0418[1],ci90.obs0418[2],ci90.obs0418[2],ci90.obs0418[1]),
            col='gray95',border=NA)  
    arrows(x0=par('usr')[1],x1=par('usr')[2],y0=obsmn.0418,y1=obsmn.0418,length=0,col='gray40',lty=2)
    axis(1,at=par('usr')[1:2],tcl=0,labels=FALSE)
    axis(2,at=par('usr')[3:4],tcl=0,labels=FALSE)
    axis(3,at=par('usr')[1:2],tcl=0,labels=FALSE)
    axis(4,at=par('usr')[3:4],tcl=0,labels=FALSE)
    arrows(x0=pss$scenario[pss$period==2],x1=pss$scenario[pss$period==2],y0=pss$q0.05[pss$period==2],y1=pss$q0.95[pss$period==2],
           length=0,col=figcol4,lwd=1.2)
    arrows(x0=pss$scenario[pss$period==2],x1=pss$scenario[pss$period==2],y0=pss$q0.25[pss$period==2],y1=pss$q0.75[pss$period==2],
           length=0,col=figcol4,lwd=2.5)
    points(q0.5~scenario,data=pss[pss$period==2,],pch=21,col=figcol4,bg=figcol4)
    mtext('2050-2070',side=1,line=0.2,cex=0.9)
  plot(0,0,type='n',xlim=c(0,max(pss$scenario)+1),ylim=c(0,9.6),axes=FALSE,xlab='',ylab='')
    polygon(x=c(par('usr')[1],par('usr')[2],par('usr')[2],par('usr')[1],par('usr')[1]),
            y=c(ci90.obs0418[1],ci90.obs0418[1],ci90.obs0418[2],ci90.obs0418[2],ci90.obs0418[1]),
            col='gray95',border=NA)  
    arrows(x0=par('usr')[1],x1=par('usr')[2],y0=obsmn.0418,y1=obsmn.0418,length=0,col='gray40',lty=2)
    axis(1,at=par('usr')[1:2],tcl=0,labels=FALSE)
    axis(2,at=par('usr')[3:4],tcl=0,labels=FALSE)
    axis(3,at=par('usr')[1:2],tcl=0,labels=FALSE)
    axis(4,at=par('usr')[3:4],tcl=0,labels=FALSE)
    arrows(x0=pss$scenario[pss$period==3],x1=pss$scenario[pss$period==3],y0=pss$q0.05[pss$period==3],y1=pss$q0.95[pss$period==3],
           length=0,col=figcol4,lwd=1.2)
    arrows(x0=pss$scenario[pss$period==3],x1=pss$scenario[pss$period==3],y0=pss$q0.25[pss$period==3],y1=pss$q0.75[pss$period==3],
           length=0,col=figcol4,lwd=2.5)
    points(q0.5~scenario,data=pss[pss$period==3,],pch=21,col=figcol4,bg=figcol4)
    mtext('2080-2100',side=1,line=0.2,cex=0.9)
  mtext('Total area occupied (ha)',side=2,line=2.0,outer=TRUE,cex=0.9)

#-----------------------------------------------------------------------------# 
#Summarizing annual variation in the forecasted area occupied

  match <- data.frame(summary=c('min','max','mn','md','sd'),summary.no=1:5,stringsAsFactors=FALSE)
  preds.annmin <- lapply(win.list,function(x) c(x[1,1:3],apply(x[,5:(n_iter+4)],2,min)))
  preds.annmax <- lapply(win.list,function(x) c(x[1,1:3],apply(x[,5:(n_iter+4)],2,max)))
  preds.annmn <- lapply(win.list,function(x) c(x[1,1:3],apply(x[,5:(n_iter+4)],2,mean)))
  preds.annmd <- lapply(win.list,function(x) c(x[1,1:3],apply(x[,5:(n_iter+4)],2,median)))
  preds.annsd <- lapply(win.list,function(x) c(x[1,1:3],apply(x[,5:(n_iter+4)],2,sd)))
  
  preds.annmin <- lapply(1:length(preds.annmin), function(x) c(summary.no=1,preds.annmin[[x]]))
  preds.annmax <- lapply(1:length(preds.annmax), function(x) c(summary.no=2,preds.annmax[[x]]))
  preds.annmn <- lapply(1:length(preds.annmn), function(x) c(summary.no=3,preds.annmn[[x]]))
  preds.annmd <- lapply(1:length(preds.annmd), function(x) c(summary.no=4,preds.annmd[[x]]))
  preds.annsd <- lapply(1:length(preds.annsd), function(x) c(summary.no=5,preds.annsd[[x]]))

#Dataframe with the min/max/mn/md/sd of the forecasted annual area occupied (preds.ann)
#(min/max/mn/md/sd among years for each scenario, GCM, period and iteration)
  preds.ann <- rbind(do.call(rbind,preds.annmin),
                     do.call(rbind,preds.annmax),
                     do.call(rbind,preds.annmn),
                     do.call(rbind,preds.annmd),
                     do.call(rbind,preds.annsd))
  preds.ann <- as.data.frame(preds.ann)
  preds.ann <- join(preds.ann,match,by='summary.no',type='left')
  preds.ann <- preds.ann[,c(2:4,ncol(preds.ann),5:(ncol(preds.ann)-1))]
 
#Then, calculate the median of those values across iterations
#(ie, on average, the minimum/maximum/mean/median/sd of annual forecasted values across 21 years for each scenario, GCM, period)
  preds.annvar.l <- preds.ann[,1:4]
  preds.annvar.l$md <- apply(preds.ann[,5:ncol(preds.ann)],1,median)
  preds.annvar <- dcast(preds.annvar.l,scenario + model + period ~ summary,value.var='md')
  preds.annvar <- preds.annvar[,c(1:3,6,4,7,5,8)]

#-----------------------------------------------------------------------------#   
#Calculate the probability that one or more years in each period falls below a given threshold
  th <- 0.67  #0.67 = minimum observed value
  #Extracting minimum areas for each scenario, GCM, period, and iteration
  mins <- preds.ann[preds.ann$summary=='min',]
  #Create table with median probability that at least one year per period is below threshold
  #with min and max probabilities among GCMs
  mins$pbelow <- apply(mins[,5:1004],1,function(x) sum(x<th)/length(x))
  probbelow <- ddply(mins,.(period,scenario),summarize,prob.md=round(median(pbelow),2),
                     prob.min=round(min(pbelow),2),prob.max=round(max(pbelow),2))
  probbelow$scenarioname <- sc$scenarioname[match(probbelow$scenario,sc$scenario)]  

#-----------------------------------------------------------------------------# 
#Forecast the number of years per period that fall below a given threshold
  #First, create dataframe with number of years below threshold for each scenario, GCM, period, and iteration
  preds.low <- lapply(win.list,function(x) c(x[1,1:3],apply(x[,5:(n_iter+4)],2,function(xx) sum(xx<th))))
  preds.lowyrs <- as.data.frame(do.call(rbind,preds.low))
  #Then calculate the probable number of years for each scenario and period (averaged across GCMs)
  index.s <- paste(preds.lowyrs$scenario,preds.lowyrs$period,sep='-')
  lowyrs.list.s <- split(preds.lowyrs,f=index.s)
  lowyrs.summary.s <- unique(preds.lowyrs[,c('scenario','period')])
  lowyrs.matlist.s <- lapply(lowyrs.list.s,'[',4:(n_iter+3))
  lowyrs.veclist.s <- lapply(lowyrs.matlist.s,unlist)
  lowyrs.summary.s$q0.5 <- sapply(lowyrs.veclist.s,median)
  lowyrs.summary.s$q0.05 <- sapply(lowyrs.veclist.s,quantile,0.05)
  lowyrs.summary.s$q0.25 <- sapply(lowyrs.veclist.s,quantile,0.25)
  lowyrs.summary.s$q0.75 <- sapply(lowyrs.veclist.s,quantile,0.75)
  lowyrs.summary.s$q0.95 <- sapply(lowyrs.veclist.s,quantile,0.95)
  pls <- lowyrs.summary.s  

#-----------------------------------------------------------------------------#  
#Plot with forecasted number of years below threshold (medians, 50 & 90% CIs), averaged across GCMs for each emission scenario
#(Figs. 2d-f in paper)

  par(mfrow=c(1,3),mar=c(0,0,0.3,0),oma=c(1.4,3.5,0.2,0.5),cex=0.9)
  plot(0,0,type='n',xlim=c(0,max(pls$scenario)+1),ylim=c(0,21),axes=FALSE,xlab='',ylab='')
    axis(1,at=par('usr')[1:2],tcl=0,labels=FALSE)
    axis(2,at=par('usr')[3:4],tcl=0,labels=FALSE)
    axis(2,at=seq(0,20,by=5),tcl=-0.25,labels=seq(0,20,by=5),las=1,mgp=c(1.5,0.5,0))
    axis(3,at=par('usr')[1:2],tcl=0,labels=FALSE)
    axis(4,at=par('usr')[3:4],tcl=0,labels=FALSE)
    arrows(x0=pls$scenario[pls$period==1],x1=pls$scenario[pls$period==1],y0=pls$q0.05[pls$period==1],y1=pls$q0.95[pls$period==1],
           length=0,col=figcol4,lwd=1.2)
    arrows(x0=pls$scenario[pls$period==1],x1=pls$scenario[pls$period==1],y0=pls$q0.25[pls$period==1],y1=pls$q0.75[pls$period==1],
           length=0,col=figcol4,lwd=2.5)    
    points(q0.5~scenario,data=pls[pls$period==1,],pch=21,col=figcol4,bg=figcol4)
    mtext('2023-2043',side=1,line=0.2,cex=0.9)
    legend(x=0,y=21,c('Low','Mod-Low','Mod-High','High'),pch=21,col=figcol4,pt.bg=figcol4,bty='n',cex=0.8,)
  plot(0,0,type='n',xlim=c(0,max(pls$scenario)+1),ylim=c(0,21),axes=FALSE,xlab='',ylab='')
    axis(1,at=par('usr')[1:2],tcl=0,labels=FALSE)
    axis(2,at=par('usr')[3:4],tcl=0,labels=FALSE)
    axis(3,at=par('usr')[1:2],tcl=0,labels=FALSE)
    axis(4,at=par('usr')[3:4],tcl=0,labels=FALSE)
    arrows(x0=pls$scenario[pls$period==2],x1=pls$scenario[pls$period==2],y0=pls$q0.05[pls$period==2],y1=pls$q0.95[pls$period==2],
           length=0,col=figcol4,lwd=1.2)
    arrows(x0=pls$scenario[pls$period==2],x1=pls$scenario[pls$period==2],y0=pls$q0.25[pls$period==2],y1=pls$q0.75[pls$period==2],
           length=0,col=figcol4,lwd=2.5)
    points(q0.5~scenario,data=pls[pls$period==2,],pch=21,col=figcol4,bg=figcol4)
    mtext('2050-2070',side=1,line=0.2,cex=0.9)
  plot(0,0,type='n',xlim=c(0,max(pls$scenario)+1),ylim=c(0,21),axes=FALSE,xlab='',ylab='')
    axis(1,at=par('usr')[1:2],tcl=0,labels=FALSE)
    axis(2,at=par('usr')[3:4],tcl=0,labels=FALSE)
    axis(3,at=par('usr')[1:2],tcl=0,labels=FALSE)
    axis(4,at=par('usr')[3:4],tcl=0,labels=FALSE)
    arrows(x0=pls$scenario[pls$period==3],x1=pls$scenario[pls$period==3],y0=pls$q0.05[pls$period==3],y1=pls$q0.95[pls$period==3],
           length=0,col=figcol4,lwd=1.2)
    arrows(x0=pls$scenario[pls$period==3],x1=pls$scenario[pls$period==3],y0=pls$q0.25[pls$period==3],y1=pls$q0.75[pls$period==3],
           length=0,col=figcol4,lwd=2.5)
    points(q0.5~scenario,data=pls[pls$period==3,],pch=21,col=figcol4,bg=figcol4)
    mtext('2080-2100',side=1,line=0.2,cex=0.9)
  mtext('Number of years below 0.67 ha',side=2,line=2.0,outer=TRUE,cex=0.9)

#---------------------------------------------------------------------------------------------------------#
# Partitioning uncertainty
#---------------------------------------------------------------------------------------------------------#  

#-----------------------------------------------------------------------------#    
#Forecasts that don't incorporate environmental stochasticity
#(but do account for parameter uncertainty and climate uncertainty)

  win.list.noes <- list()
  i <- 0
  set.seed(44)
  
  for(ss in 1:nrow(sc)){ #scenario

    for(mm in 1:nrow(mod)){ #GCM

      for(pp in 1:3){ #time period
        uyears <- get(paste0('years',pp))
        i <- i+1
        
        #Combine ecological covariates in summer submodel in dataframe
          cyw <- expand.grid(wk=uweeks,county.ind=ucounties,yr=uyears)
          cyw$yr.ind <- cyw$yr - min(cyw$yr) + 1
            #Standardize week (want same values as those used in 2004-2018 analysis)
            cyw$wk.st <- (cyw$wk-wk.m)/wk.sd
          cyw <- join(cyw,springclim[springclim$scenario==ss & springclim$model==mm & springclim$period==pp,
                                     c('yr','spGDD.st','spGDD.st2','spPCP.st','spPCP.st2')],by='yr',type='left')
          cyw <- join_all(list(cyw,cov.c[,c('county.ind','cropC.st')],gly.c[,c('county.ind','gly.st')],
                               avgGDD[avgGDD$scenario==ss & avgGDD$model==mm & avgGDD$period==pp,c('county.ind','avgGDD.st')],
                               avgPCP[avgPCP$scenario==ss & avgPCP$model==mm & avgPCP$period==pp,c('county.ind','avgPCP.st')]),
                          by='county.ind',type='left')
          cyw <- join(cyw,diffPCP[diffPCP$scenario==ss & diffPCP$model==mm & diffPCP$period==pp,c('county.ind','yr','diffPCP.st','diffPCP.st2')],
                      by=c('county.ind','yr'),type='left')
          cyw <- join(cyw,diffGDD[diffGDD$scenario==ss & diffGDD$model==mm & diffGDD$period==pp,c('county.ind','yr','wk','diffGDD.st','diffGDD.st2')],
                      by=c('county.ind','yr','wk'),type='left')
          cyw$diffavgGDD <- cyw$diffGDD.st*cyw$avgGDD.st
          cyw$diffavgPCP <- cyw$diffPCP.st*cyw$avgPCP.st
          cyw$glycrop <- cyw$gly.st*cyw$cropC.st
          cyw <- cyw[,c('county.ind','yr','yr.ind','wk','wk.st','spGDD.st','spGDD.st2','spPCP.st','spPCP.st2',
                        'avgGDD.st','diffGDD.st','diffGDD.st2','diffavgGDD',
                        'avgPCP.st','diffPCP.st','diffPCP.st2','diffavgPCP',
                        'gly.st','cropC.st','glycrop')]
          names(cyw)[6:ncol(cyw)] <- c('spGDD','spGDD2','spPCP','spPCP2',
                                       'avgGDD','diffGDD','diffGDD2','diffavgGDD',
                                       'avgPCP','diffPCP','diffPCP2','diffavgPCP',
                                       'gly','crop','glycrop')

        #Calculate summer population index from posterior samples
          #Fixed effects
          Xcounty <- as.matrix(cbind(rep(1,nrow(cyw)),
                                     cyw[,c('spGDD','spGDD2','spPCP','spPCP2',
                                            'avgGDD','diffGDD','diffGDD2','diffavgGDD',
                                            'avgPCP','diffPCP','diffPCP2','diffavgPCP',
                                            'gly','crop','glycrop')]))
          Xalpha <- Xcounty %*% t(alpha)

          #Effect of week
          Xweeks <- matrix(c(cyw$wk.st,cyw$wk.st*cyw$wk.st),ncol=2,nrow=nrow(cyw),byrow=FALSE)  
          alphaREweeks <- cbind(alpha.wk,alpha.wk2)
          Xalpha.wks <- Xweeks %*% t(alphaREweeks)

          #Random county effects
          eps.county <- RE.county1yr[rep(1:nrow(RE.county1yr),times=n_years),]

          #Combine components
          mu.county.log <- Xalpha + Xalpha.wks + eps.county

          #Calculate means for each county-year combination
          split.muc <- lapply(split(mu.county.log,rep(1:(n_counties*n_years),each=n_weeks)),matrix,ncol=ncol(mu.county.log))
          mu.cy <- t(sapply(split.muc,colMeans))

          #Calculate weighted mean of counties in each year
          weights.mat <- matrix(weights.open,nrow=nrow(mu.cy),ncol=ncol(mu.cy),byrow=FALSE)  
          mu.cy.weights <- mu.cy*weights.mat
          split.mucy <- lapply(split(mu.cy.weights,rep(1:n_years,each=n_counties)),matrix,ncol=ncol(mu.cy.weights))
          mu.y <- t(sapply(split.mucy,colSums))
    
          #Standardize summer index by fixed mean, SD
          mu.y.z <- (mu.y - 1.12)/0.59

        #Calculate area occupied from posterior samples  
          Wsum <- matrix(NA,nrow=n_years,ncol=n_iter)

          #Summer effect 
            for(ii in 1:n_iter){
              Wsum[,ii] <- mu.y.z[,ii]*gamma.sum[ii,1]
            }

          #Combine components and exponentiate
            mu.win <- exp(Wgamma0 + Wsum)
            
        #Add matrix of values to list
          params <- matrix(as.numeric(c(ss,mm,pp)),nrow=n_years,ncol=3,byrow=TRUE)
          params <- cbind(params,uyears)
          win.list.noes[[i]] <- cbind(params,mu.win)
          colnames(win.list.noes[[i]]) <- c('scenario','model','period','yr',paste0('i',1:n_iter))
          
        #See progress
          print(i)

      } #pp
    } #mm
  } #ss

#-----------------------------------------------------------------------------#    
#Forecasts that don't account for environmental stochasticity nor parameter uncertainty
#(only account for climate uncertainty)
  
  win.list.clim <- list()
  i <- 0
  set.seed(44)

  for(ss in 1:nrow(sc)){ #scenario

    for(mm in 1:nrow(mod)){ #GCM

      for(pp in 1:3){ #time period
        uyears <- get(paste0('years',pp))
        i <- i+1

        #Combine ecological covariates in summer submodel in dataframe
          cyw <- expand.grid(wk=uweeks,county.ind=ucounties,yr=uyears)
          cyw$yr.ind <- cyw$yr - min(cyw$yr) + 1
            #Standardize week (want same values as those used in 2004-2018 analysis)
            cyw$wk.st <- (cyw$wk-wk.m)/wk.sd
          cyw <- join(cyw,springclim[springclim$scenario==ss & springclim$model==mm & springclim$period==pp,
                                     c('yr','spGDD.st','spGDD.st2','spPCP.st','spPCP.st2')],by='yr',type='left')
          cyw <- join_all(list(cyw,cov.c[,c('county.ind','cropC.st')],gly.c[,c('county.ind','gly.st')],
                               avgGDD[avgGDD$scenario==ss & avgGDD$model==mm & avgGDD$period==pp,c('county.ind','avgGDD.st')],
                               avgPCP[avgPCP$scenario==ss & avgPCP$model==mm & avgPCP$period==pp,c('county.ind','avgPCP.st')]),
                          by='county.ind',type='left')
          cyw <- join(cyw,diffPCP[diffPCP$scenario==ss & diffPCP$model==mm & diffPCP$period==pp,c('county.ind','yr','diffPCP.st','diffPCP.st2')],
                      by=c('county.ind','yr'),type='left')
          cyw <- join(cyw,diffGDD[diffGDD$scenario==ss & diffGDD$model==mm & diffGDD$period==pp,c('county.ind','yr','wk','diffGDD.st','diffGDD.st2')],
                      by=c('county.ind','yr','wk'),type='left')
          cyw$diffavgGDD <- cyw$diffGDD.st*cyw$avgGDD.st
          cyw$diffavgPCP <- cyw$diffPCP.st*cyw$avgPCP.st
          cyw$glycrop <- cyw$gly.st*cyw$cropC.st
          cyw <- cyw[,c('county.ind','yr','yr.ind','wk','wk.st','spGDD.st','spGDD.st2','spPCP.st','spPCP.st2',
                        'avgGDD.st','diffGDD.st','diffGDD.st2','diffavgGDD',
                        'avgPCP.st','diffPCP.st','diffPCP.st2','diffavgPCP',
                        'gly.st','cropC.st','glycrop')]
          names(cyw)[6:ncol(cyw)] <- c('spGDD','spGDD2','spPCP','spPCP2',
                                       'avgGDD','diffGDD','diffGDD2','diffavgGDD',
                                       'avgPCP','diffPCP','diffPCP2','diffavgPCP',
                                       'gly','crop','glycrop')
 
        #Calculate summer population index from posterior samples
          #Fixed effects
          Xcounty <- as.matrix(cbind(rep(1,nrow(cyw)),
                                     cyw[,c('spGDD','spGDD2','spPCP','spPCP2',
                                            'avgGDD','diffGDD','diffGDD2','diffavgGDD',
                                            'avgPCP','diffPCP','diffPCP2','diffavgPCP',
                                            'gly','crop','glycrop')]))
          Xalpha <- Xcounty %*% alpha.md

          #Effect of week
          Xweeks <- matrix(c(cyw$wk.st,cyw$wk.st*cyw$wk.st),ncol=2,nrow=nrow(cyw),byrow=FALSE)  
          alphaREweeks <- matrix(c(alpha.wk.md,alpha.wk2.md),nrow=2)
          Xalpha.wks <- Xweeks %*% alphaREweeks
          
          #Random county effects
          eps.county <- RE.countyallyr.md

          #Combine components
          mu.county.log <- Xalpha + Xalpha.wks + eps.county

          #Calculate means for each county-year combination
          split.muc <- split(mu.county.log,rep(1:(n_counties*n_years),each=n_weeks))
          mu.cy <- sapply(split.muc,mean)

          #Calculate weighted mean of counties in each year
          weights.vec <- rep(weights.open,n_years)  
          mu.cy.weights <- mu.cy*weights.vec
          split.mucy <- split(mu.cy.weights,rep(1:n_years,each=n_counties))
          mu.y <- sapply(split.mucy,sum)
    
          #Standardize summer index by fixed mean, SD
          mu.y.z <- (mu.y - 1.12)/0.59

        #Calculate area occupied from posterior samples  
          #Summer effect 
            Wsum <- mu.y.z*gamma.sum.md

          #Combine components and exponentiate
            mu.win <- exp(Wgamma0.md + Wsum)
            
        #Add matrix of values to list (could add stuff so we also output summer index values on the log scale [mu.y])
          params <- matrix(as.numeric(c(ss,mm,pp)),nrow=n_years,ncol=3,byrow=TRUE)
          params <- cbind(params,uyears)
          win.list.clim[[i]] <- cbind(params,mu.win)
          colnames(win.list.clim[[i]]) <- c('scenario','model','period','yr','mu.win')
          
        #See progress
          print(i)

      } #pp
    } #mm
  } #ss
  
#-----------------------------------------------------------------------------#    
#Forecasts that don't account for environmental stochasticity and don't account for climate uncertainty
#(only account for parameter uncertainty)
  #Using CNRMESM model and moderate-high emissions scenario (SSP370)

  win.list.param33 <- list()
  i <- 0
  set.seed(44)
  
    ss <- 3

      mm <- 3

      for(pp in 1:3){ #time period
        selectp <- pp
        uyears <- get(paste0('years',selectp))
        i <- i+1
        
        #Combine ecological covariates in summer submodel in dataframe
          cyw <- expand.grid(wk=uweeks,county.ind=ucounties,yr=uyears)
          cyw$yr.ind <- cyw$yr - min(cyw$yr) + 1
            #Standardize week (want same values as those used in 2004-2018 analysis)
            cyw$wk.st <- (cyw$wk-wk.m)/wk.sd
          cyw <- join(cyw,springclim[springclim$scenario==ss & springclim$model==mm & springclim$period==pp,
                                     c('yr','spGDD.st','spGDD.st2','spPCP.st','spPCP.st2')],by='yr',type='left')
          cyw <- join_all(list(cyw,cov.c[,c('county.ind','cropC.st')],gly.c[,c('county.ind','gly.st')],
                               avgGDD[avgGDD$scenario==ss & avgGDD$model==mm & avgGDD$period==pp,c('county.ind','avgGDD.st')],
                               avgPCP[avgPCP$scenario==ss & avgPCP$model==mm & avgPCP$period==pp,c('county.ind','avgPCP.st')]),
                          by='county.ind',type='left')
          cyw <- join(cyw,diffPCP[diffPCP$scenario==ss & diffPCP$model==mm & diffPCP$period==pp,c('county.ind','yr','diffPCP.st','diffPCP.st2')],
                      by=c('county.ind','yr'),type='left')
          cyw <- join(cyw,diffGDD[diffGDD$scenario==ss & diffGDD$model==mm & diffGDD$period==pp,c('county.ind','yr','wk','diffGDD.st','diffGDD.st2')],
                      by=c('county.ind','yr','wk'),type='left')
          cyw$diffavgGDD <- cyw$diffGDD.st*cyw$avgGDD.st
          cyw$diffavgPCP <- cyw$diffPCP.st*cyw$avgPCP.st
          cyw$glycrop <- cyw$gly.st*cyw$cropC.st
          cyw <- cyw[,c('county.ind','yr','yr.ind','wk','wk.st','spGDD.st','spGDD.st2','spPCP.st','spPCP.st2',
                        'avgGDD.st','diffGDD.st','diffGDD.st2','diffavgGDD',
                        'avgPCP.st','diffPCP.st','diffPCP.st2','diffavgPCP',
                        'gly.st','cropC.st','glycrop')]
          names(cyw)[6:ncol(cyw)] <- c('spGDD','spGDD2','spPCP','spPCP2',
                                       'avgGDD','diffGDD','diffGDD2','diffavgGDD',
                                       'avgPCP','diffPCP','diffPCP2','diffavgPCP',
                                       'gly','crop','glycrop')

      #Calculate summer population index from posterior samples
        #Fixed effects
        Xcounty <- as.matrix(cbind(rep(1,nrow(cyw)),
                                   cyw[,c('spGDD','spGDD2','spPCP','spPCP2',
                                          'avgGDD','diffGDD','diffGDD2','diffavgGDD',
                                          'avgPCP','diffPCP','diffPCP2','diffavgPCP',
                                          'gly','crop','glycrop')]))
        Xalpha <- Xcounty %*% t(alpha)
        
        #Effect of week
        Xweeks <- matrix(c(cyw$wk.st,cyw$wk.st*cyw$wk.st),ncol=2,nrow=nrow(cyw),byrow=FALSE)  
        alphaREweeks <- cbind(alpha.wk,alpha.wk2)
        Xalpha.wks <- Xweeks %*% t(alphaREweeks)
        
        #Random county effects
        eps.county <- RE.county1yr[rep(1:nrow(RE.county1yr),times=n_years),]

        #Combine components
        mu.county.log <- Xalpha + Xalpha.wks + eps.county
        
        #Calculate means for each county-year combination
        split.muc <- lapply(split(mu.county.log,rep(1:(n_counties*n_years),each=n_weeks)),matrix,ncol=ncol(mu.county.log))
        mu.cy <- t(sapply(split.muc,colMeans))
        
        #Calculate weighted mean of counties in each year
        weights.mat <- matrix(weights.open,nrow=nrow(mu.cy),ncol=ncol(mu.cy),byrow=FALSE)  
        mu.cy.weights <- mu.cy*weights.mat
        split.mucy <- lapply(split(mu.cy.weights,rep(1:n_years,each=n_counties)),matrix,ncol=ncol(mu.cy.weights))
        mu.y <- t(sapply(split.mucy,colSums))
        
        #Standardize summer index by fixed mean, SD
        mu.y.z <- (mu.y - 1.12)/0.59
        
      #Calculate area occupied from posterior samples  
        Wsum <- matrix(NA,nrow=n_years,ncol=n_iter)
        
        #Summer effect 
        for(ii in 1:n_iter){
          Wsum[,ii] <- mu.y.z[,ii]*gamma.sum[ii,1]
        }
        
        #Combine components and exponentiate
        mu.win <- exp(Wgamma0 + Wsum)
        
        #Add matrix of values to list
        params <- matrix(as.numeric(c(ss,mm,pp)),nrow=n_years,ncol=3,byrow=TRUE)
        params <- cbind(params,uyears)
        win.list.param33[[i]] <- cbind(params,mu.win)
        colnames(win.list.param33[[i]]) <- c('scenario','model','period','yr',paste0('i',1:n_iter))
        
        #See progress
        print(i)
        
      } #pp
   
#-----------------------------------------------------------------------------#       
#Summarize forecasts and compare width of 90% credible intervals

  #No environmental stochasticity
    preds.byPGR.noes <- lapply(win.list.noes,colMeans)
    preds.byPGR2.noes <- do.call(rbind,preds.byPGR.noes)
    preds.byPGR.df.noes <- as.data.frame(preds.byPGR2.noes)
    index.p <- preds.byPGR.df.noes$period
    preds.p.list.noes <- split(preds.byPGR.df.noes,f=index.p)
    preds.summary.noes <- data.frame(period=unique(preds.byPGR.df.noes$period))
    preds.summary.noes$mn <- sapply(preds.p.list.noes,function(x) mean(as.matrix(x[,5:(n_iter+4)])))
    preds.summary.noes$sd <- sapply(preds.p.list.noes,function(x) sd(as.matrix(x[,5:(n_iter+4)])))
    preds.summary.noes$q0.05 <- sapply(preds.p.list.noes,function(x) quantile(as.matrix(x[,5:(n_iter+4)]),0.05))
    preds.summary.noes$q0.25 <- sapply(preds.p.list.noes,function(x) quantile(as.matrix(x[,5:(n_iter+4)]),0.25))
    preds.summary.noes$q0.5 <- sapply(preds.p.list.noes,function(x) quantile(as.matrix(x[,5:(n_iter+4)]),0.5))
    preds.summary.noes$q0.75 <- sapply(preds.p.list.noes,function(x) quantile(as.matrix(x[,5:(n_iter+4)]),0.75))
    preds.summary.noes$q0.95 <- sapply(preds.p.list.noes,function(x) quantile(as.matrix(x[,5:(n_iter+4)]),0.95))
    preds.summary.noes$CIw90 <- preds.summary.noes$q0.95 - preds.summary.noes$q0.05      
  
  #No environmental stochasticity, no parameter uncertainty (climate uncertainty only)
    preds.byPGR.clim <- lapply(win.list.clim,colMeans)
    preds.byPGR2.clim <- do.call(rbind,preds.byPGR.clim)
    preds.byPGR.df.clim <- as.data.frame(preds.byPGR2.clim)
    index.p <- preds.byPGR.df.clim$period
    preds.p.list.clim <- split(preds.byPGR.df.clim,f=index.p)
    preds.summary.clim <- data.frame(period=unique(preds.byPGR.df.clim$period))
    preds.summary.clim$mn <- sapply(preds.p.list.clim,function(x) mean(x[,5]))
    preds.summary.clim$sd <- sapply(preds.p.list.clim,function(x) sd(x[,5]))
    preds.summary.clim$q0.05 <- sapply(preds.p.list.clim,function(x) quantile(x[,5],0.05))
    preds.summary.clim$q0.25 <- sapply(preds.p.list.clim,function(x) quantile(x[,5],0.25))
    preds.summary.clim$q0.5 <- sapply(preds.p.list.clim,function(x) quantile(x[,5],0.5))
    preds.summary.clim$q0.75 <- sapply(preds.p.list.clim,function(x) quantile(x[,5],0.75))
    preds.summary.clim$q0.95 <- sapply(preds.p.list.clim,function(x) quantile(x[,5],0.95))
    preds.summary.clim$CIw90 <- preds.summary.clim$q0.95 - preds.summary.clim$q0.05      
  
  #No environmental stochasticity, no climate uncertainty (parameter uncertainty only)   
    preds.p.list.param33 <- lapply(win.list.param33,colMeans)
    preds.summary.param33 <- data.frame(period=1:3)
    preds.summary.param33$mn <- sapply(preds.p.list.param33,function(x) mean(x[5:(n_iter+4)]))
    preds.summary.param33$sd <- sapply(preds.p.list.param33,function(x) sd(x[5:(n_iter+4)]))
    preds.summary.param33$q0.05 <- sapply(preds.p.list.param33,function(x) quantile(x[5:(n_iter+4)],0.05))
    preds.summary.param33$q0.25 <- sapply(preds.p.list.param33,function(x) quantile(x[5:(n_iter+4)],0.25))
    preds.summary.param33$q0.5 <- sapply(preds.p.list.param33,function(x) quantile(x[5:(n_iter+4)],0.5))
    preds.summary.param33$q0.75 <- sapply(preds.p.list.param33,function(x) quantile(x[5:(n_iter+4)],0.75))
    preds.summary.param33$q0.95 <- sapply(preds.p.list.param33,function(x) quantile(x[5:(n_iter+4)],0.95))
    preds.summary.param33$CIw90 <- preds.summary.param33$q0.95 - preds.summary.param33$q0.05   
    
  #All uncertainties    
    preds.byPGR <- lapply(win.list,colMeans)
    preds.byPGR2 <- do.call(rbind,preds.byPGR)
    preds.byPGR.df <- as.data.frame(preds.byPGR2)
    index.p <- preds.byPGR.df$period
    preds.p.list <- split(preds.byPGR.df,f=index.p)
    preds.summary.allu <- data.frame(period=unique(preds.byPGR.df$period))
    preds.summary.allu$mn <- sapply(preds.p.list,function(x) mean(as.matrix(x[,5:(n_iter+4)])))
    preds.summary.allu$sd <- sapply(preds.p.list,function(x) sd(as.matrix(x[,5:(n_iter+4)])))
    preds.summary.allu$q0.05 <- sapply(preds.p.list,function(x) quantile(as.matrix(x[,5:(n_iter+4)]),0.05))
    preds.summary.allu$q0.25 <- sapply(preds.p.list,function(x) quantile(as.matrix(x[,5:(n_iter+4)]),0.25))
    preds.summary.allu$q0.5 <- sapply(preds.p.list,function(x) quantile(as.matrix(x[,5:(n_iter+4)]),0.5))
    preds.summary.allu$q0.75 <- sapply(preds.p.list,function(x) quantile(as.matrix(x[,5:(n_iter+4)]),0.75))
    preds.summary.allu$q0.95 <- sapply(preds.p.list,function(x) quantile(as.matrix(x[,5:(n_iter+4)]),0.95))
    preds.summary.allu$CIw90 <- preds.summary.allu$q0.95 - preds.summary.allu$q0.05
  
  #Combining all results
    preds.summary.allu$forecast <- 'all'
    preds.summary.noes$forecast <- 'noes'
    preds.summary.clim$forecast <- 'climonly'
    preds.summary.param33$forecast <- 'puonly'
    psumm <- rbind(preds.summary.allu,preds.summary.noes,preds.summary.clim,preds.summary.param33)
    
    #Calculating the mean of forecast medians in each period
    mns <- ddply(psumm,.(period),summarize,mn=mean(q0.5))
    psumm$periodmn <- mns$mn[match(psumm$period,mns$period)]
    #Relativizing to 0
    psumm$q0.05R0 <- psumm$q0.05-psumm$q0.5
    psumm$q0.95R0 <- psumm$q0.95-psumm$q0.5
    psumm$q0.25R0 <- psumm$q0.25-psumm$q0.5
    psumm$q0.75R0 <- psumm$q0.75-psumm$q0.5
 
#Figure with 4 sets of forecasts: 50 and 90% CIs, centered around zero
#(Fig. 3 in the paper)
  
    par(mfrow=c(1,1),mar=c(1.5,3,0.7,0.7),cex=0.9) 
    #climate uncertainty only
    plot(0,0,type='n',xlim=c(0.5,3.5),ylim=c(-3.8,4.1),axes=F,xaxs='i',yaxs='i')
      arrows(x0=psumm$period[psumm$forecast=='all']-0.2,x1=psumm$period[psumm$forecast=='all']-0.2,
             y0=psumm$q0.05R0[psumm$forecast=='climonly'],y1=psumm$q0.95R0[psumm$forecast=='climonly'],
             length=0,col='darkseagreen4',lwd=1)
      arrows(x0=psumm$period[psumm$forecast=='all']-0.2,x1=psumm$period[psumm$forecast=='all']-0.2,
             y0=psumm$q0.25R0[psumm$forecast=='climonly'],y1=psumm$q0.75R0[psumm$forecast=='climonly'],
             length=0,col='darkseagreen4',lwd=3)
      axis(1,at=par('usr')[1:2],labels=F,tcl=0)
      axis(1,at=1:3,labels=c('2023-2043','2050-2070','2080-2100'),tcl=-0.25,mgp=c(1.5,0.4,0),cex=0.9)
      axis(2,at=par('usr')[3:4],labels=F,tcl=0)
      axis(2,at=seq(-4,4,by=2),labels=seq(-4,4,by=2),tcl=-0.25,mgp=c(1.5,0.5,0),las=1,cex=0.9) 
    #parameter uncertainty only  
      arrows(x0=psumm$period[psumm$forecast=='all']-0.1,x1=psumm$period[psumm$forecast=='all']-0.1,
             y0=psumm$q0.05R0[psumm$forecast=='puonly'],y1=psumm$q0.95R0[psumm$forecast=='puonly'],
             length=0,col='salmon3',lwd=1)
      arrows(x0=psumm$period[psumm$forecast=='all']-0.1,x1=psumm$period[psumm$forecast=='all']-0.1,
             y0=psumm$q0.25R0[psumm$forecast=='puonly'],y1=psumm$q0.75R0[psumm$forecast=='puonly'],
             length=0,col='salmon3',lwd=3)    
    #climate + parameter uncertainty
      arrows(x0=psumm$period[psumm$forecast=='all']-0,x1=psumm$period[psumm$forecast=='all']-0,
             y0=psumm$q0.05R0[psumm$forecast=='noes'],y1=psumm$q0.95R0[psumm$forecast=='noes'],
             length=0,col='steelblue3',lwd=1)
      arrows(x0=psumm$period[psumm$forecast=='all']-0,x1=psumm$period[psumm$forecast=='all']-0,
             y0=psumm$q0.25R0[psumm$forecast=='noes'],y1=psumm$q0.75R0[psumm$forecast=='noes'],
             length=0,col='steelblue3',lwd=3)  
    #all uncertainties 
      arrows(x0=psumm$period[psumm$forecast=='all']+0.1,x1=psumm$period[psumm$forecast=='all']+0.1,
             y0=psumm$q0.05R0[psumm$forecast=='all'],y1=psumm$q0.95R0[psumm$forecast=='all'],
             length=0,col='black',lwd=1)
      arrows(x0=psumm$period[psumm$forecast=='all']+0.1,x1=psumm$period[psumm$forecast=='all']+0.1,
             y0=psumm$q0.25R0[psumm$forecast=='all'],y1=psumm$q0.75R0[psumm$forecast=='all'],
             length=0,col='black',lwd=3)  
      arrows(x0=par('usr')[1],x1=par('usr')[2],y0=0,y1=0,col='gray50',length=0)
      legend(x=0.5,y=-1.3,legend=c('Climate','Parameter','Climate + Parameter','Climate + Parameter + Environment'),
             lty=1,lwd=2,col=c('darkseagreen4','salmon3','steelblue3','black'),bty='n',cex=0.9)
      mtext('Area occupied (ha), relativized',side=2,line=1.7)

