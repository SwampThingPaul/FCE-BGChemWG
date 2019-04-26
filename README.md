# FCE LTER Biogeochemistry Working Group <img src="./Hex/FCE_BGChem.png" align="right" height = "120"/>

This repo is dedicated to the sharing of information between [Biogeochemistry Working Group](http://fcelter.fiu.edu/research/working_groups/?wg=20&p=FCEIII) team members. 

# Task Lists

- [x] Develop a hex sticker for the FCE LTER Biogeochemistry working group (NERD ALERT!!). 
- [x] Develop R function to retrieve data directly from LTER main data portal, [https://portal.lternet.edu](https://portal.lternet.edu/nis/home.jsp).
- [ ] Develop R-script to retrieve, screen and plot data along Shark River and Taylor Slough transects. 
- [ ] Develop a R-shiny app to interact with BGChem data.

***

## R-function
Here is an extremely basic function to pull data from the LTER main data portal. Dataset data package, PASTA and DOI is needed. 

```
# Function to retrieve LTER data from the main data portal. 

read.lter=function(data.package,PASTA,DOI){
  prefix="http://pasta.lternet.edu/package/data/eml/"
  
  infile1=paste0(prefix,data.package,PASTA,DOI)
  dt1=read.csv(infile1,na.strings=c("-9999","-9999.00","-9999.000"))
  return(dt1)
}

###########
## Example
###########

PASTA="1075/8/"
DOI="ac7159e66cbad75abb61bd1992f8d2c0"
example=read.lter("knb-lter-fce/",PASTA,DOI)

```

## Trends
The ["../Trends"](./Trends/) sub-directory has scripts specifically analyzing trend data within the FCE monitoring network. 

Here are preliminary plots.

<img src="./Trends/Plots/TPTNDOC_WQPlots.png" align="center" width = "100%"/>

<center> Total Phosphorus, Total Nitrogen and Dissolved Organic Carbon concentration data for sites along Shark River Slough (SRS) and Taylor Slough (TS).</center>

<br>
<br>

<img src="./Trends/Plots/SRPDIN_WQPlots.png" align="center" width = "100%"/>

<center> Soluable Reactive Phosphorus and Dissolved Inorganic Nitrogen concentration data for sites along Shark River Slough (SRS) and Taylor Slough (TS).</center>