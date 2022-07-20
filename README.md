#  Multi-season climate projections forecast declines in migratory monarch butterflies

### Erin R. Zylstra, Naresh Neupane, and Elise F. Zipkin

### Global Change Biology

### Code/Data DOI: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5597600.svg)](https://doi.org/10.5281/zenodo.5597600)
_______________________________________________________________________________________________________________________________________

## Abstract:
Climate change poses a unique threat to migratory species as it has the potential to alter environmental conditions at multiple points along a species’ migratory route. The eastern migratory population of monarch butterflies (<i>Danaus plexippus</i>) has declined markedly over the last few decades, in part due to variation in breeding-season climate. Here, we combined a retrospective, annual-cycle model for the eastern monarch population with climate projections within the spring breeding grounds in eastern Texas and across the summer breeding grounds in the midwestern U.S. and southern Ontario, Canada to evaluate how monarchs are likely to respond to climate change over the next century. Our results reveal that projected changes in breeding-season climate are likely to lead to decreases in monarch abundance, with high potential for overwintering population size to fall below the historical minimum three or more times in the next two decades. Climatic changes across the expansive summer breeding grounds will also cause shifts in the distribution of monarchs, with higher projected abundances in areas that become wetter but not appreciably hotter (e.g., northern Ohio) and declines in abundance where summer temperatures are projected to increase well above those observed in the recent past (e.g., northern Minnesota). Although climate uncertainties dominate long-term population forecasts, our analyses suggest that we can improve precision of near-term forecasts by collecting targeted data to better understand relationships between breeding-season climate variables and local monarch abundance. Overall, our results highlight the importance of accounting for the impacts of climate changes throughout the full annual cycle of migratory species. 
_______________________________________________________________________________________________________________________________________

## Code
1. [MonarchPopulationForecasts.R](Code/MonarchPopulationForecasts.R): R code to run a retrospective full-annual-cycle model for monarchs in eastern North America from 2004-2018, import and format climate projections for the spring and summer breeding grounds, and forecast monarch population sizes in summer and early winter. 
2. [FACM_2004-2018.stan](Code/FACM_2004-2018.stan): Stan model file for running the retrospective population model (called from the R script above). 

## Data for the retrospective population model
All monarch data from the overwintering grounds and covariate data for 2004-2018 are publicly available.  Monarch data from the summer breeding grounds are proprietary and are therefore not publicly available (though we provide descriptions of those data here).  
1. [YearlyData.csv](Data/Monarchs_winter.csv): Data on overwintering monarch population size and annual covariates. 
    - yr: year 
    - area.dec: area occupied (ha) by monarch butterflies in December
    - spGDD.east: annual growing degree days (GDD) accumulated between 22 Mar and 2 May (weeks 4-9) in eastern Texas
    - spPCP.east: cumulative precipitation (mm) between Feb and Apr in eastern Texas
    - NDVI: mean NDVI for the first (northern) half of autumn migration
    - denseforest: percent of area in and around overwintering colonies with dense forest cover
2. [Covariates_County.txt](Data/Covariates_County.txt): Covariate values associated with each county on the summer breeding range, in the U.S. and Canada.
    - county.ind: a unique index identifying each county (1:545)
    - state.county: combined state-county FIPS code
    - state: FIPS code for each state
    - county: FIPS code for each county in a state
    - state.name: 2-letter code for each state
    - county.name: name of each county
    - area.land.sqmi: land area of each county (sq. mi)
    - Xcentroid: longitude of county centroid
    - Ycentroid: latitude of county centroid
    - avgGDD: average of annual GDD values accumulated between 3 May and 15 Aug (weeks 10-24), 2004-2018
    - avgPCP: average of annual, cumulative precipitation (mm) between Apr and Aug, 2004-2018
    - perc.open: percent of each county that is unforested
    - perc.crop: percent of each county associated with agricultural crops
    - surveyed: indicates whether one or more monarch surveys were conducted in the county between 2004 and 2018 or not (1 or 0, respectively)
3. [Covariates_CountyYear.csv](Data/Covariates_CountyYear.csv): Annual covariate values associated with each county on the summer breeding range, in the U.S. and Canada, 2004-2018.
    - county.ind: a unique index identifying each county (1:545)
    - state.county: combined state-county FIPS code
    - yr: year
    - glyphosate: estimated proportion of corn and soy crops in each county sprayed with glyphosate
    - diffPCP: difference between annual precipitation (cumulative, Apr-Aug) and average precipitation (average of annual values, 2004-2018) in each county
4. [Covariates_CountyWeek.csv](Data/Covariates_CountyWeek.csv): Weekly covariate values associated with each county on the summer breeding range, in the U.S. and Canada, 2004-2018.
    - county.ind: a unique index identifying each county (1:545)
    - state.county: combined state-county FIPS code
    - yr: year
    - wk: week (16-24)
    - diffGDD: difference between GDD and average GDD (2004-2018) for that week and county
5. Monarchs_summer.csv: Data from monarch surveys on the summer breeding grounds, 2004-2018 (not publicly available). 
    - program: name of monitoring program
    - state.county: combined state-county FIPS code
    - lat: latitude of survey location
    - long: longitude of survey location
    - yr: year
    - wk: week (16-24)
    - monarch: total number of adult monarchs observed during survey
    - duration: number of person/party hours spent surveying
    - site.ind: a unique index identifying each survey location
    - county.ind: a unique index identifying each county
    - perc.open: percent of area immediately surrounding each survey location that is unforested 

## Climate projections
1. [SpringClimateProjections.csv](Data/SpringClimateProjections.csv): Projected GDD (accumulated between 22 Mar and 2 May) and precipitation (Feb through Apr) for the spring breeding grounds in eastern Texas for three future time periods
    - Scenario: Combinations of Shared Socioeconomic Pathways (SSPs) and Representative Concentration Pathways (RCPs) that describe different trends in atmospheric greenhouse gases (SSP126 = low emissions scenario; SSP245 = moderate emissions scenario; SSP370 = moderate-high emissions scenario; SSP585 = high emissions scenario)
    - Year: year
    - Model: GCM
    - spGDD: projected GDD
    - pcp.Feb: projected precipitation in February
    - pcp.Mar: projected precipitation in March
    - pcp.Apr: projected precipitation in April
2. [SummerClimateProjections1a.csv](Data/SummerClimateProjections1a.csv): Projected GDD and precipitation (Apr through Aug) for each county on the summer breeding grounds for the first of three future time periods (2023-2043) and 3 of the 6 GCMs.
    - Scenario: Combinations of Shared Socioeconomic Pathways (SSPs) and Representative Concentration Pathways (RCPs) that describe different trends in atmospheric greenhouse gases.
    - Year: year
    - Model: GCM
    - state.county: combined state-county FIPS code
    - gdd.wk10.21: projected GDD accumulated between 3 May and 25 Jul
    - gdd.wk10.22: projected GDD accumulated between 3 May and 1 Aug
    - gdd.wk10.23: projected GDD accumulated between 3 May and 8 Aug
    - gdd.wk10.24: projected GDD accumulated between 3 May and 15 Aug
    - pcp.Apr: projected precipitation in April
    - pcp.May: projected precipitation in May
    - pcp.Jun: projected precipitation in June
    - pcp.Jul: projected precipitation in July
    - pcp.Aug: projected precipitation in August
3. [SummerClimateProjections1b.csv](Data/SummerClimateProjections2.csv): Projected GDD and precipitation (Apr through Aug) for each county on the summer breeding grounds for the early twenty-first century (2023-2043) and the other 3 GCMs. Columns are the same as those in SummerClimateProjections1a.csv.
4. [SummerClimateProjections2a.csv](Data/SummerClimateProjections2a.csv): Projected GDD and precipitation (Apr through Aug) for each county on the summer breeding grounds for the mid twenty-first century (2050-2070) and 3 of the 6 GCMs. Columns are the same as those in SummerClimateProjections1a.csv.
5. [SummerClimateProjections2b.csv](Data/SummerClimateProjections2b.csv): Projected GDD and precipitation (Apr through Aug) for each county on the summer breeding grounds for the mid twenty-first century (2050-2070) and the other 3 GCMs. Columns are the same as those in SummerClimateProjections1a.csv.
6. [SummerClimateProjections3a.csv](Data/SummerClimateProjectionsa3.csv): Projected GDD and precipitation (Apr through Aug) for each county on the summer breeding grounds for the late twenty-first century (2080-2100) and 3 of the 6 GCMs. Columns are the same as those in SummerClimateProjections1a.csv.
7. [SummerClimateProjections3b.csv](Data/SummerClimateProjectionsa3b.csv): Projected GDD and precipitation (Apr through Aug) for each county on the summer breeding grounds for the late twenty-first century (2080-2100) and the other 3 GCMs. Columns are the same as those in SummerClimateProjections1a.csv.
