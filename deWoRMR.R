#deWoRMR: automatic web scraping of taxonomic hierarchies from WoRMS
#Written by Nicolas Mongiardino Koch 11/2023

#ONLY INPUT REQUIRED: URL to the page of the target clade (starting point)
#In this case Echinoidea, but search for your clade of interest in WoRMS and
#paste the correct link below
starting_url = 'https://www.marinespecies.org/aphia.php?p=taxdetails&id=123082'

#Install (optionally) and load required packages
packages <- c('rvest','stringr','dplyr')
new_packages <- packages[!packages %in% installed.packages()[,'Package']]
if(length(new_packages)) { install.packages(new_packages) }

library(rvest)
library(stringr)
library(dplyr)

################################  PART 1  ######################################
###    Traverse the taxonomic hierarchy of WoRMS and finds URLs of species   ###
################################################################################

read_url = read_html(starting_url)

#this is the URL address of all WoRMS pages
worms_base = 'https://www.marinespecies.org/aphia.php?p=taxdetails&id='

#links_to_visit contains links to child taxa that need to be searched
links_to_visit = html_attr(html_nodes(read_url, '#ChildTaxa a'), 'href')
links_to_visit = sapply(strsplit(links_to_visit, 'id='), '[', 2)

#vector to save the the information on species
species = c()
#place to save aphia IDs of places we've already been to
already_visited = sapply(strsplit(starting_url, 'id='), '[', 2)

#loop that keeps going as long as there are websites that still need visiting
i = 1
while(i <= length(links_to_visit)) {
  #if we have not already visited this website, go in, otherwise simply delete
  if(!links_to_visit[i] %in% already_visited) {
    url = paste0(worms_base, links_to_visit[i])
    read_url = read_html(url)
    
    #what kind of taxon are we visiting?
    taxon_level = str_squish(gsub("[\r\n]", "", 
                                  html_text(html_nodes(read_url, 
                                                       '#Rank .leave_image_space'))))
    
    #if it is not a fossil (sometimes these show up inconsistently)
    if(length(taxon_level) > 0) {
      #if it is not a species, check to see if it has direct children and add
      #these to the list of links we need to visit
      if(taxon_level != 'Species') {
        if(length(html_attr(html_nodes(read_url, '#ChildTaxa a'), 'href')) > 0) {
          to_add = html_attr(html_nodes(read_url, '#ChildTaxa a'), 'href')
          links_to_visit = c(links_to_visit, 
                             sapply(strsplit(to_add, 'id='), '[', 2))
        }
        
        #if it is a species, download its information
      } else {
        status = html_text(html_nodes(read_url, '#Status span'))
        
        if(status %in% c('accepted', 'temporary name')) {
          correct_species = html_text(html_nodes(read_url, 
                                                 '#AcceptedName .leave_image_space'))
        } else {
          correct_species = 'NA'
        }
        
        #get the taxonomy and clean
        this_taxonomy = html_text(html_nodes(read_url, 
                                             '#Classification .leave_image_space'))
        this_taxonomy = str_squish(gsub("[\r\n]", "", this_taxonomy))
        taxon_info = paste(links_to_visit[i], this_taxonomy, 
                           status, correct_species, sep = ', ')
        
        #save if not an ichnotaxon
        if(!grepl('ichnotaxa', taxon_info)) {
          species = c(species, taxon_info)
        }
      }
    }
  }

  #no matter what you did, send current link from those waiting to be visited
  #to those already visited
  already_visited = c(already_visited, links_to_visit[i])
  links_to_visit = links_to_visit[-i]
}

## Optional (but this took a while so probably safer)
save(species, file = 'taxonomy_raw.Rda')

################################  PART 2  ######################################
###                   Give the data a more useful format                     ###
################################################################################

taxonomy = sapply(strsplit(species, ', '), '[', 2)

#find a list of subgenera
subgenera = c()
for(i in 1:length(taxonomy)) {
  this_taxonomy_splitted = unlist(strsplit(taxonomy[i], ' '))
  if(any(this_taxonomy_splitted == '(Subgenus)')) {
    if(any(this_taxonomy_splitted == '†')) {
      pos = which(this_taxonomy_splitted == '(Subgenus)') - 2
    } else {
      pos = which(this_taxonomy_splitted == '(Subgenus)') - 1
    }
    subgenera = c(subgenera, this_taxonomy_splitted[pos])
  }
}
subgenera = unique(subgenera)

#get all text strings included within brackets, and remove the subgenera to
#retain only taxonomic ranks. Also remove very infrequently used ranks (present
#in < 5% of taxa)
hierarchy = unlist(strsplit(taxonomy, ' '))[grep('\\(', unlist(strsplit(taxonomy, ' ')))]
hierarchy = sort(table(hierarchy), decreasing = T)
hierarchy = hierarchy[!names(hierarchy) %in% subgenera]
hierarchy = names(hierarchy[hierarchy > max(hierarchy)*0.05])

df = as_tibble(matrix(NA, ncol = length(hierarchy), 
                      nrow = length(species)), .name_repair = 'unique')
df = df %>% rename_at(vars(colnames(df)), ~ hierarchy) %>% 
  mutate_if(is.logical, as.character)

#fill in the taxonomic df based on the rank hierarchy
for(i in 1:length(taxonomy)) {
  a = unlist(strsplit(taxonomy[i], ' '))
  matches = which(a %in% colnames(df))
  for(j in length(matches):1) {
    if(j != 1) {
      to_use = paste(a[(matches[j-1]+1):(matches[j]-1)], collapse = ' ')
    } else {
      to_use = paste(a[1:(matches[j]-1)], collapse = ' ')
    }
    df[i, which(colnames(df) == a[matches[j]])] = to_use
  }
}

#add other informatio gathered besides taxonomy
aphiaID = sapply(strsplit(species, ', '), '[', 1)
status = sapply(strsplit(species, ', '), '[', 3)
valid_species = sapply(strsplit(species, ', '), '[', 4)
valid_species = str_squish(gsub("[\r\n]", "", valid_species))
valid_species[is.na(valid_species)] = 'NA'

#get rid of authorities for valid species
for(i in 1:length(valid_species)) {
  if(valid_species[i] != 'NA') {
    clean_sp = valid_species[i]
    while(!clean_sp %in% df$`(Species)`) {
      clean_sp = unlist(strsplit(clean_sp, ' '))
      
      #couldn't be cleaned of authorities, keep full name
      if(length(clean_sp) == 1) {
        clean_sp = valid_species[i]
        break
      }
      
      clean_sp = paste(clean_sp[1:(length(clean_sp)-1)], collapse = ' ')
    }
    valid_species[i] = clean_sp
  }
}

#combine all elements
final_taxonomy = df %>% mutate(aphiaID = aphiaID, Status = status, 
                               Valid_species = valid_species) %>%
  dplyr::select(aphiaID, everything())

#make sure we didn't include taxa from other clades (trust me, it happens)
read_url = read_html(starting_url)
this_taxonomy = html_text(html_nodes(read_url, 
                                     '#Classification .leave_image_space'))
this_taxonomy = str_squish(gsub("[\r\n]", "", this_taxonomy))
this_taxonomy = unlist(strsplit(this_taxonomy, ' '))

#filter off target hits and rename columns to remove parentheses
target_rank = tail(this_taxonomy, n = 1)
target_clade = tail(which(this_taxonomy %in% colnames(df)), n = 2)
target_clade = paste0(this_taxonomy[which(1:length(this_taxonomy) > target_clade[1] & 
                                            1:length(this_taxonomy) < target_clade[2])], 
                      collapse = ' ')

final_taxonomy = final_taxonomy %>% 
  filter(!!as.symbol(target_rank) == target_clade) %>% 
  rename_at(vars(colnames(df)), ~ gsub('^.|.$', '', colnames(df)))

#remove columns that are constant (i.e., taxonomic categories above that of
#the clade of interest)
constant = names(which(apply(final_taxonomy, 2, 
                             function(x) length(unique(x))) == 1))
final_taxonomy = final_taxonomy %>% select(-one_of(constant))

#attempt at ordering columns following their position in the taxonomic hierarchy
to_sort = colnames(final_taxonomy)
to_sort = to_sort[2:(length(to_sort)-2)]

#remove from original data the ones that turned out to be off target and split
#to find the relative position of taxonomic ranks
taxonomy = taxonomy[which(aphiaID %in% unlist(final_taxonomy$aphiaID))]
split_taxonomy = strsplit(taxonomy, ' ')
for(i in 1:length(split_taxonomy)) {
  split_taxonomy[[i]] = split_taxonomy[[i]][grep('\\(', split_taxonomy[[i]])]
  split_taxonomy[[i]] = split_taxonomy[[i]][split_taxonomy[[i]] %in% 
                                              paste0('(', 
                                                     colnames(final_taxonomy), 
                                                     ')')]
}

#start with a species that had a long hierarchy and try to place any missing
#ranks relative to those already present within it
longest = which.max(sapply(split_taxonomy, length))
starting_order = split_taxonomy[[longest]]
still_missing = to_sort[which(!to_sort %in% gsub('^.|.$', '', starting_order))]
while(length(starting_order) != length(to_sort)) {
  for(i in 1:length(still_missing)) {
    smaller_split_taxonomy = split_taxonomy[which(sapply(split_taxonomy, function(x) 
      any(grepl(still_missing[i], x))))]
    
    pos = sapply(smaller_split_taxonomy, function(x) grep(still_missing[i], x))
    before = unlist(Map('[', smaller_split_taxonomy, pos-1))
    
    #if the missing rank is always before or always after one who's position we
    #already know, the it is added to its correct placement
    if(length(unique(before)) == 1) {
      goes_after = gsub('^.|.$', '', unique(before))
      starting_order = c(starting_order[1:grep(goes_after, starting_order)], 
                         paste0('(', still_missing[i], ')'),
                         starting_order[(grep(goes_after, starting_order)+1):
                                          length(starting_order)])
    } else {
      after = unlist(Map('[', smaller_split_taxonomy, pos+1))
      if(length(unique(after)) == 1) {
        goes_before = gsub('^.|.$', '', unique(after))
        starting_order = c(starting_order[1:(grep(goes_before, starting_order)-1)], 
                           paste0('(', still_missing[i], ')'),
                           starting_order[grep(goes_after, starting_order):
                                            length(starting_order)])
      } else {
        #inferring the position failed, will be left at the end
        starting_order = c(starting_order, 
                           paste0('(', still_missing[i], ')'))
      }
    }
  }
}

#sort the final taxonomy
final_positions = as.numeric(na.omit(match(gsub('^.|.$', '', starting_order), 
                                           colnames(final_taxonomy))))
final_positions = c(1, final_positions, (ncol(final_taxonomy)-1):ncol(final_taxonomy))
ordered_names = colnames(final_taxonomy)[final_positions]

final_taxonomy = final_taxonomy %>% relocate(any_of(ordered_names)) %>% 
  mutate(across(where(is.character), ~na_if(., 'NA'))) %>% arrange(Species)

#save results as csv and Rda files
name = paste0(target_clade, '_WoRMS_taxonomy_', Sys.Date())
write.csv(final_taxonomy, file = paste0(name, '.csv'), quote = F, row.names = F)
save(final_taxonomy, file = paste0(name, '.Rda'))