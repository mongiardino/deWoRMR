# deWoRMR: automatic web scraping of taxonomic hierarchies from WoRMS
This script takes as only input the URL of a website for a clade of interest in [WoRMS](https://www.marinespecies.org/index.php) (the World Register of Marine Species) and uses web scraping tools to download the taxonomy of the clade and save it in ways that can be easily manipulated. The final product is a table that includes every species within the clade as a row entry, and every higher-level taxon those species are assigned to as the columns. Further information such as taxon validity and synonyms is stored.

It should be noted that WoRMS has developed its own ways of accessing and downloading their data (see [here](https://www.marinespecies.org/aphia.php?p=webservice)), yet I prefer using this simple script which will find, extract, and format the entire taxonomic tree for a clade with minimal user input. As the code visits **every** website in the hierarchy downstream of the taxon of interest, running times can become considerable. As an example, scraping the taxonomy of Echinoidea (sea urchins, heart urchins, sand dollars, and close relatives), which includes 1,012 valid species as of 11/14/2023, requires visiting 6,936 websites (including 2,023 invalid species names and all higher-level taxa). This process takes about 1.5 hours, but this number will ultimately also depend on internet connectivity. More efficient options should definitely be explored for clades that are much larger than this.
