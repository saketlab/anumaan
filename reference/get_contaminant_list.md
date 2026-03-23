# Get Contaminant List from Reference File

Loads contaminant organisms from common_commensals.csv and filters by
syndrome and/or specimen type.

## Usage

``` r
get_contaminant_list(syndrome = NULL, specimen_type = NULL, return_all = FALSE)
```

## Arguments

- syndrome:

  Character. Optional syndrome name to filter by (e.g., "Bloodstream
  infections", "Urinary tract infections").

- specimen_type:

  Character. Optional specimen type to filter by (e.g., "Blood culture",
  "Urine culture").

- return_all:

  Logical. If TRUE, returns all contaminants from all syndromes. Default
  FALSE.

## Value

Character vector of contaminant organism names (lowercase)

## Examples

``` r
# Get all blood culture contaminants
get_contaminant_list(syndrome = "Bloodstream infections")
#> $names
#>  [1] "Coagulase-negative Staphylococcus" "Staphylococcus epidermidis"       
#>  [3] "Staphylococcus hominis"            "Staphylococcus haemolyticus"      
#>  [5] "Staphylococcus capitis"            "Staphylococcus warneri"           
#>  [7] "Corynebacterium species"           "Corynebacterium striatum"         
#>  [9] "Corynebacterium jeikeium"          "Cutibacterium"                    
#> [11] "Cutibacterium acnes"               "Micrococcus species"              
#> [13] "Micrococcus luteus"                "Bacillus (non-anthracis)"         
#> [15] "Bacillus subtilis"                 "Bacillus cereus"                  
#> [17] "Viridans streptococci"             "Aerococcus spp."                  
#> [19] "Kocuria spp."                      "Dermacoccus spp."                 
#> [21] "Rothia spp."                      
#> 
#> $patterns
#> $patterns[[1]]
#> $patterns[[1]]$original
#> [1] "Coagulase-negative Staphylococcus"
#> 
#> $patterns[[1]]$patterns
#> [1] "coagulase-negative staphylococcus" "^c\\.\\s"                         
#> [3] "^c\\.\\s*staphylococcus"          
#> 
#> 
#> $patterns[[2]]
#> $patterns[[2]]$original
#> [1] "Staphylococcus epidermidis"
#> 
#> $patterns[[2]]$patterns
#> [1] "staphylococcus epidermidis" "^s\\.\\s"                  
#> [3] "epidermidis"                "^s\\.\\s*epidermidis"      
#> [5] "staph.*epidermidis"         "coag.*neg"                 
#> [7] "coagulase.*negative"       
#> 
#> 
#> $patterns[[3]]
#> $patterns[[3]]$original
#> [1] "Staphylococcus hominis"
#> 
#> $patterns[[3]]$patterns
#> [1] "staphylococcus hominis" "^s\\.\\s"               "hominis"               
#> [4] "^s\\.\\s*hominis"       "staph.*hominis"         "coag.*neg"             
#> [7] "coagulase.*negative"   
#> 
#> 
#> $patterns[[4]]
#> $patterns[[4]]$original
#> [1] "Staphylococcus haemolyticus"
#> 
#> $patterns[[4]]$patterns
#> [1] "staphylococcus haemolyticus" "^s\\.\\s"                   
#> [3] "haemolyticus"                "^s\\.\\s*haemolyticus"      
#> [5] "staph.*haemolyticus"         "coag.*neg"                  
#> [7] "coagulase.*negative"        
#> 
#> 
#> $patterns[[5]]
#> $patterns[[5]]$original
#> [1] "Staphylococcus capitis"
#> 
#> $patterns[[5]]$patterns
#> [1] "staphylococcus capitis" "^s\\.\\s"               "capitis"               
#> [4] "^s\\.\\s*capitis"       "staph.*capitis"         "coag.*neg"             
#> [7] "coagulase.*negative"   
#> 
#> 
#> $patterns[[6]]
#> $patterns[[6]]$original
#> [1] "Staphylococcus warneri"
#> 
#> $patterns[[6]]$patterns
#> [1] "staphylococcus warneri" "^s\\.\\s"               "warneri"               
#> [4] "^s\\.\\s*warneri"       "staph.*warneri"         "coag.*neg"             
#> [7] "coagulase.*negative"   
#> 
#> 
#> $patterns[[7]]
#> $patterns[[7]]$original
#> [1] "Corynebacterium species"
#> 
#> $patterns[[7]]$patterns
#> [1] "corynebacterium species" "corynebacterium"        
#> [3] "^c\\.\\s"                "^c\\.\\s*species"       
#> [5] "coryno"                  "diphtheroids?"          
#> 
#> 
#> $patterns[[8]]
#> $patterns[[8]]$original
#> [1] "Corynebacterium striatum"
#> 
#> $patterns[[8]]$patterns
#> [1] "corynebacterium striatum" "^c\\.\\s"                
#> [3] "striatum"                 "^c\\.\\s*striatum"       
#> [5] "coryno"                   "diphtheroids?"           
#> 
#> 
#> $patterns[[9]]
#> $patterns[[9]]$original
#> [1] "Corynebacterium jeikeium"
#> 
#> $patterns[[9]]$patterns
#> [1] "corynebacterium jeikeium" "^c\\.\\s"                
#> [3] "jeikeium"                 "^c\\.\\s*jeikeium"       
#> [5] "coryno"                   "diphtheroids?"           
#> 
#> 
#> $patterns[[10]]
#> $patterns[[10]]$original
#> [1] "Cutibacterium"
#> 
#> $patterns[[10]]$patterns
#> [1] "cutibacterium"  "^c\\.\\s"       "propioni"       "cuti"          
#> [5] "p\\.?\\s*acnes"
#> 
#> 
#> $patterns[[11]]
#> $patterns[[11]]$original
#> [1] "Cutibacterium acnes"
#> 
#> $patterns[[11]]$patterns
#> [1] "cutibacterium acnes" "^c\\.\\s"            "acnes"              
#> [4] "^c\\.\\s*acnes"      "propioni"            "cuti"               
#> [7] "p\\.?\\s*acnes"     
#> 
#> 
#> $patterns[[12]]
#> $patterns[[12]]$original
#> [1] "Micrococcus species"
#> 
#> $patterns[[12]]$patterns
#> [1] "micrococcus species" "micrococcus"         "^m\\.\\s"           
#> [4] "^m\\.\\s*species"    "micro"              
#> 
#> 
#> $patterns[[13]]
#> $patterns[[13]]$original
#> [1] "Micrococcus luteus"
#> 
#> $patterns[[13]]$patterns
#> [1] "micrococcus luteus" "^m\\.\\s"           "luteus"            
#> [4] "^m\\.\\s*luteus"    "micro"             
#> 
#> 
#> $patterns[[14]]
#> $patterns[[14]]$original
#> [1] "Bacillus (non-anthracis)"
#> 
#> $patterns[[14]]$patterns
#> [1] "bacillus (non-anthracis)" "^b\\.\\s"                
#> [3] "(non-anthracis)"          "^b\\.\\s*(non-anthracis)"
#> [5] "bacil"                   
#> 
#> 
#> $patterns[[15]]
#> $patterns[[15]]$original
#> [1] "Bacillus subtilis"
#> 
#> $patterns[[15]]$patterns
#> [1] "bacillus subtilis" "^b\\.\\s"          "subtilis"         
#> [4] "^b\\.\\s*subtilis" "bacil"            
#> 
#> 
#> $patterns[[16]]
#> $patterns[[16]]$original
#> [1] "Bacillus cereus"
#> 
#> $patterns[[16]]$patterns
#> [1] "bacillus cereus" "^b\\.\\s"        "cereus"          "^b\\.\\s*cereus"
#> [5] "bacil"          
#> 
#> 
#> $patterns[[17]]
#> $patterns[[17]]$original
#> [1] "Viridans streptococci"
#> 
#> $patterns[[17]]$patterns
#> [1] "viridans streptococci" "^v\\.\\s"              "streptococci"         
#> [4] "^v\\.\\s*streptococci"
#> 
#> 
#> $patterns[[18]]
#> $patterns[[18]]$original
#> [1] "Aerococcus spp."
#> 
#> $patterns[[18]]$patterns
#> [1] "aerococcus spp." "aerococcus"      "^a\\.\\s"        "^a\\.\\s*spp."  
#> 
#> 
#> $patterns[[19]]
#> $patterns[[19]]$original
#> [1] "Kocuria spp."
#> 
#> $patterns[[19]]$patterns
#> [1] "kocuria spp."  "kocuria"       "^k\\.\\s"      "^k\\.\\s*spp."
#> 
#> 
#> $patterns[[20]]
#> $patterns[[20]]$original
#> [1] "Dermacoccus spp."
#> 
#> $patterns[[20]]$patterns
#> [1] "dermacoccus spp." "dermacoccus"      "^d\\.\\s"         "^d\\.\\s*spp."   
#> 
#> 
#> $patterns[[21]]
#> $patterns[[21]]$original
#> [1] "Rothia spp."
#> 
#> $patterns[[21]]$patterns
#> [1] "rothia spp."   "rothia"        "^r\\.\\s"      "^r\\.\\s*spp."
#> 
#> 
#> 

# Get all contaminants for a specific specimen type
get_contaminant_list(specimen_type = "Blood culture")
#> $names
#>  [1] "Coagulase-negative Staphylococcus" "Staphylococcus epidermidis"       
#>  [3] "Staphylococcus hominis"            "Staphylococcus haemolyticus"      
#>  [5] "Staphylococcus capitis"            "Staphylococcus warneri"           
#>  [7] "Corynebacterium species"           "Corynebacterium striatum"         
#>  [9] "Corynebacterium jeikeium"          "Cutibacterium"                    
#> [11] "Cutibacterium acnes"               "Micrococcus species"              
#> [13] "Micrococcus luteus"                "Bacillus (non-anthracis)"         
#> [15] "Bacillus subtilis"                 "Bacillus cereus"                  
#> [17] "Viridans streptococci"             "Aerococcus spp."                  
#> [19] "Kocuria spp."                      "Dermacoccus spp."                 
#> [21] "Rothia spp."                      
#> 
#> $patterns
#> $patterns[[1]]
#> $patterns[[1]]$original
#> [1] "Coagulase-negative Staphylococcus"
#> 
#> $patterns[[1]]$patterns
#> [1] "coagulase-negative staphylococcus" "^c\\.\\s"                         
#> [3] "^c\\.\\s*staphylococcus"          
#> 
#> 
#> $patterns[[2]]
#> $patterns[[2]]$original
#> [1] "Staphylococcus epidermidis"
#> 
#> $patterns[[2]]$patterns
#> [1] "staphylococcus epidermidis" "^s\\.\\s"                  
#> [3] "epidermidis"                "^s\\.\\s*epidermidis"      
#> [5] "staph.*epidermidis"         "coag.*neg"                 
#> [7] "coagulase.*negative"       
#> 
#> 
#> $patterns[[3]]
#> $patterns[[3]]$original
#> [1] "Staphylococcus hominis"
#> 
#> $patterns[[3]]$patterns
#> [1] "staphylococcus hominis" "^s\\.\\s"               "hominis"               
#> [4] "^s\\.\\s*hominis"       "staph.*hominis"         "coag.*neg"             
#> [7] "coagulase.*negative"   
#> 
#> 
#> $patterns[[4]]
#> $patterns[[4]]$original
#> [1] "Staphylococcus haemolyticus"
#> 
#> $patterns[[4]]$patterns
#> [1] "staphylococcus haemolyticus" "^s\\.\\s"                   
#> [3] "haemolyticus"                "^s\\.\\s*haemolyticus"      
#> [5] "staph.*haemolyticus"         "coag.*neg"                  
#> [7] "coagulase.*negative"        
#> 
#> 
#> $patterns[[5]]
#> $patterns[[5]]$original
#> [1] "Staphylococcus capitis"
#> 
#> $patterns[[5]]$patterns
#> [1] "staphylococcus capitis" "^s\\.\\s"               "capitis"               
#> [4] "^s\\.\\s*capitis"       "staph.*capitis"         "coag.*neg"             
#> [7] "coagulase.*negative"   
#> 
#> 
#> $patterns[[6]]
#> $patterns[[6]]$original
#> [1] "Staphylococcus warneri"
#> 
#> $patterns[[6]]$patterns
#> [1] "staphylococcus warneri" "^s\\.\\s"               "warneri"               
#> [4] "^s\\.\\s*warneri"       "staph.*warneri"         "coag.*neg"             
#> [7] "coagulase.*negative"   
#> 
#> 
#> $patterns[[7]]
#> $patterns[[7]]$original
#> [1] "Corynebacterium species"
#> 
#> $patterns[[7]]$patterns
#> [1] "corynebacterium species" "corynebacterium"        
#> [3] "^c\\.\\s"                "^c\\.\\s*species"       
#> [5] "coryno"                  "diphtheroids?"          
#> 
#> 
#> $patterns[[8]]
#> $patterns[[8]]$original
#> [1] "Corynebacterium striatum"
#> 
#> $patterns[[8]]$patterns
#> [1] "corynebacterium striatum" "^c\\.\\s"                
#> [3] "striatum"                 "^c\\.\\s*striatum"       
#> [5] "coryno"                   "diphtheroids?"           
#> 
#> 
#> $patterns[[9]]
#> $patterns[[9]]$original
#> [1] "Corynebacterium jeikeium"
#> 
#> $patterns[[9]]$patterns
#> [1] "corynebacterium jeikeium" "^c\\.\\s"                
#> [3] "jeikeium"                 "^c\\.\\s*jeikeium"       
#> [5] "coryno"                   "diphtheroids?"           
#> 
#> 
#> $patterns[[10]]
#> $patterns[[10]]$original
#> [1] "Cutibacterium"
#> 
#> $patterns[[10]]$patterns
#> [1] "cutibacterium"  "^c\\.\\s"       "propioni"       "cuti"          
#> [5] "p\\.?\\s*acnes"
#> 
#> 
#> $patterns[[11]]
#> $patterns[[11]]$original
#> [1] "Cutibacterium acnes"
#> 
#> $patterns[[11]]$patterns
#> [1] "cutibacterium acnes" "^c\\.\\s"            "acnes"              
#> [4] "^c\\.\\s*acnes"      "propioni"            "cuti"               
#> [7] "p\\.?\\s*acnes"     
#> 
#> 
#> $patterns[[12]]
#> $patterns[[12]]$original
#> [1] "Micrococcus species"
#> 
#> $patterns[[12]]$patterns
#> [1] "micrococcus species" "micrococcus"         "^m\\.\\s"           
#> [4] "^m\\.\\s*species"    "micro"              
#> 
#> 
#> $patterns[[13]]
#> $patterns[[13]]$original
#> [1] "Micrococcus luteus"
#> 
#> $patterns[[13]]$patterns
#> [1] "micrococcus luteus" "^m\\.\\s"           "luteus"            
#> [4] "^m\\.\\s*luteus"    "micro"             
#> 
#> 
#> $patterns[[14]]
#> $patterns[[14]]$original
#> [1] "Bacillus (non-anthracis)"
#> 
#> $patterns[[14]]$patterns
#> [1] "bacillus (non-anthracis)" "^b\\.\\s"                
#> [3] "(non-anthracis)"          "^b\\.\\s*(non-anthracis)"
#> [5] "bacil"                   
#> 
#> 
#> $patterns[[15]]
#> $patterns[[15]]$original
#> [1] "Bacillus subtilis"
#> 
#> $patterns[[15]]$patterns
#> [1] "bacillus subtilis" "^b\\.\\s"          "subtilis"         
#> [4] "^b\\.\\s*subtilis" "bacil"            
#> 
#> 
#> $patterns[[16]]
#> $patterns[[16]]$original
#> [1] "Bacillus cereus"
#> 
#> $patterns[[16]]$patterns
#> [1] "bacillus cereus" "^b\\.\\s"        "cereus"          "^b\\.\\s*cereus"
#> [5] "bacil"          
#> 
#> 
#> $patterns[[17]]
#> $patterns[[17]]$original
#> [1] "Viridans streptococci"
#> 
#> $patterns[[17]]$patterns
#> [1] "viridans streptococci" "^v\\.\\s"              "streptococci"         
#> [4] "^v\\.\\s*streptococci"
#> 
#> 
#> $patterns[[18]]
#> $patterns[[18]]$original
#> [1] "Aerococcus spp."
#> 
#> $patterns[[18]]$patterns
#> [1] "aerococcus spp." "aerococcus"      "^a\\.\\s"        "^a\\.\\s*spp."  
#> 
#> 
#> $patterns[[19]]
#> $patterns[[19]]$original
#> [1] "Kocuria spp."
#> 
#> $patterns[[19]]$patterns
#> [1] "kocuria spp."  "kocuria"       "^k\\.\\s"      "^k\\.\\s*spp."
#> 
#> 
#> $patterns[[20]]
#> $patterns[[20]]$original
#> [1] "Dermacoccus spp."
#> 
#> $patterns[[20]]$patterns
#> [1] "dermacoccus spp." "dermacoccus"      "^d\\.\\s"         "^d\\.\\s*spp."   
#> 
#> 
#> $patterns[[21]]
#> $patterns[[21]]$original
#> [1] "Rothia spp."
#> 
#> $patterns[[21]]$patterns
#> [1] "rothia spp."   "rothia"        "^r\\.\\s"      "^r\\.\\s*spp."
#> 
#> 
#> 

# Get all contaminants across all syndromes
get_contaminant_list(return_all = TRUE)
#> $names
#>  [1] "Coagulase-negative Staphylococcus"                                                                                                          
#>  [2] "Staphylococcus epidermidis"                                                                                                                 
#>  [3] "Staphylococcus hominis"                                                                                                                     
#>  [4] "Staphylococcus haemolyticus"                                                                                                                
#>  [5] "Staphylococcus capitis"                                                                                                                     
#>  [6] "Staphylococcus warneri"                                                                                                                     
#>  [7] "Corynebacterium species"                                                                                                                    
#>  [8] "Corynebacterium striatum"                                                                                                                   
#>  [9] "Corynebacterium jeikeium"                                                                                                                   
#> [10] "Cutibacterium"                                                                                                                              
#> [11] "Cutibacterium acnes"                                                                                                                        
#> [12] "Micrococcus species"                                                                                                                        
#> [13] "Micrococcus luteus"                                                                                                                         
#> [14] "Bacillus (non-anthracis)"                                                                                                                   
#> [15] "Bacillus subtilis"                                                                                                                          
#> [16] "Bacillus cereus"                                                                                                                            
#> [17] "Viridans streptococci"                                                                                                                      
#> [18] "Aerococcus spp."                                                                                                                            
#> [19] "Kocuria spp."                                                                                                                               
#> [20] "Dermacoccus spp."                                                                                                                           
#> [21] "Rothia spp."                                                                                                                                
#> [22] "Same as blood (skin flora contamination during needle insertion)"                                                                           
#> [23] "Viridans group streptococci"                                                                                                                
#> [24] "Non-pathogenic Neisseria spp."                                                                                                              
#> [25] "Haemophilus parainfluenzae"                                                                                                                 
#> [26] "Anaerobic oral flora (Prevotella, Porphyromonas, Fusobacterium)"                                                                            
#> [27] "Actinomyces spp."                                                                                                                           
#> [28] "Same as blood"                                                                                                                              
#> [29] "Corynebacterium species (diphtheroids)"                                                                                                     
#> [30] "Cutibacterium (Propionibacterium) acnes"                                                                                                    
#> [31] "Bacillus species"                                                                                                                           
#> [32] "Bacteroides-Prevotella group"                                                                                                               
#> [33] "Eubacterium"                                                                                                                                
#> [34] "Lactobacillus"                                                                                                                              
#> [35] "Faecalibacterium"                                                                                                                           
#> [36] "Lactobacillus species"                                                                                                                      
#> [37] "Coagulase-negative Staphylococcus (Staphylococcus epidermidis, Staphylococcus saprophyticus)"                                               
#> [38] "Coagulase-negative Staphylococcus (Staphylococcus epidermidis, Staphylococcus hominis, Staphylococcus haemolyticus, Staphylococcus capitis)"
#> [39] "Neisseria spp. (non-pathogenic)"                                                                                                            
#> [40] "Haemophilus spp."                                                                                                                           
#> [41] "Moraxella catarrhalis"                                                                                                                      
#> [42] "No typical commensals (growth usually clinically significant)"                                                                              
#> 
#> $patterns
#> $patterns[[1]]
#> $patterns[[1]]$original
#> [1] "Coagulase-negative Staphylococcus"
#> 
#> $patterns[[1]]$patterns
#> [1] "coagulase-negative staphylococcus" "^c\\.\\s"                         
#> [3] "^c\\.\\s*staphylococcus"          
#> 
#> 
#> $patterns[[2]]
#> $patterns[[2]]$original
#> [1] "Staphylococcus epidermidis"
#> 
#> $patterns[[2]]$patterns
#> [1] "staphylococcus epidermidis" "^s\\.\\s"                  
#> [3] "epidermidis"                "^s\\.\\s*epidermidis"      
#> [5] "staph.*epidermidis"         "coag.*neg"                 
#> [7] "coagulase.*negative"       
#> 
#> 
#> $patterns[[3]]
#> $patterns[[3]]$original
#> [1] "Staphylococcus hominis"
#> 
#> $patterns[[3]]$patterns
#> [1] "staphylococcus hominis" "^s\\.\\s"               "hominis"               
#> [4] "^s\\.\\s*hominis"       "staph.*hominis"         "coag.*neg"             
#> [7] "coagulase.*negative"   
#> 
#> 
#> $patterns[[4]]
#> $patterns[[4]]$original
#> [1] "Staphylococcus haemolyticus"
#> 
#> $patterns[[4]]$patterns
#> [1] "staphylococcus haemolyticus" "^s\\.\\s"                   
#> [3] "haemolyticus"                "^s\\.\\s*haemolyticus"      
#> [5] "staph.*haemolyticus"         "coag.*neg"                  
#> [7] "coagulase.*negative"        
#> 
#> 
#> $patterns[[5]]
#> $patterns[[5]]$original
#> [1] "Staphylococcus capitis"
#> 
#> $patterns[[5]]$patterns
#> [1] "staphylococcus capitis" "^s\\.\\s"               "capitis"               
#> [4] "^s\\.\\s*capitis"       "staph.*capitis"         "coag.*neg"             
#> [7] "coagulase.*negative"   
#> 
#> 
#> $patterns[[6]]
#> $patterns[[6]]$original
#> [1] "Staphylococcus warneri"
#> 
#> $patterns[[6]]$patterns
#> [1] "staphylococcus warneri" "^s\\.\\s"               "warneri"               
#> [4] "^s\\.\\s*warneri"       "staph.*warneri"         "coag.*neg"             
#> [7] "coagulase.*negative"   
#> 
#> 
#> $patterns[[7]]
#> $patterns[[7]]$original
#> [1] "Corynebacterium species"
#> 
#> $patterns[[7]]$patterns
#> [1] "corynebacterium species" "corynebacterium"        
#> [3] "^c\\.\\s"                "^c\\.\\s*species"       
#> [5] "coryno"                  "diphtheroids?"          
#> 
#> 
#> $patterns[[8]]
#> $patterns[[8]]$original
#> [1] "Corynebacterium striatum"
#> 
#> $patterns[[8]]$patterns
#> [1] "corynebacterium striatum" "^c\\.\\s"                
#> [3] "striatum"                 "^c\\.\\s*striatum"       
#> [5] "coryno"                   "diphtheroids?"           
#> 
#> 
#> $patterns[[9]]
#> $patterns[[9]]$original
#> [1] "Corynebacterium jeikeium"
#> 
#> $patterns[[9]]$patterns
#> [1] "corynebacterium jeikeium" "^c\\.\\s"                
#> [3] "jeikeium"                 "^c\\.\\s*jeikeium"       
#> [5] "coryno"                   "diphtheroids?"           
#> 
#> 
#> $patterns[[10]]
#> $patterns[[10]]$original
#> [1] "Cutibacterium"
#> 
#> $patterns[[10]]$patterns
#> [1] "cutibacterium"  "^c\\.\\s"       "propioni"       "cuti"          
#> [5] "p\\.?\\s*acnes"
#> 
#> 
#> $patterns[[11]]
#> $patterns[[11]]$original
#> [1] "Cutibacterium acnes"
#> 
#> $patterns[[11]]$patterns
#> [1] "cutibacterium acnes" "^c\\.\\s"            "acnes"              
#> [4] "^c\\.\\s*acnes"      "propioni"            "cuti"               
#> [7] "p\\.?\\s*acnes"     
#> 
#> 
#> $patterns[[12]]
#> $patterns[[12]]$original
#> [1] "Micrococcus species"
#> 
#> $patterns[[12]]$patterns
#> [1] "micrococcus species" "micrococcus"         "^m\\.\\s"           
#> [4] "^m\\.\\s*species"    "micro"              
#> 
#> 
#> $patterns[[13]]
#> $patterns[[13]]$original
#> [1] "Micrococcus luteus"
#> 
#> $patterns[[13]]$patterns
#> [1] "micrococcus luteus" "^m\\.\\s"           "luteus"            
#> [4] "^m\\.\\s*luteus"    "micro"             
#> 
#> 
#> $patterns[[14]]
#> $patterns[[14]]$original
#> [1] "Bacillus (non-anthracis)"
#> 
#> $patterns[[14]]$patterns
#> [1] "bacillus (non-anthracis)" "^b\\.\\s"                
#> [3] "(non-anthracis)"          "^b\\.\\s*(non-anthracis)"
#> [5] "bacil"                   
#> 
#> 
#> $patterns[[15]]
#> $patterns[[15]]$original
#> [1] "Bacillus subtilis"
#> 
#> $patterns[[15]]$patterns
#> [1] "bacillus subtilis" "^b\\.\\s"          "subtilis"         
#> [4] "^b\\.\\s*subtilis" "bacil"            
#> 
#> 
#> $patterns[[16]]
#> $patterns[[16]]$original
#> [1] "Bacillus cereus"
#> 
#> $patterns[[16]]$patterns
#> [1] "bacillus cereus" "^b\\.\\s"        "cereus"          "^b\\.\\s*cereus"
#> [5] "bacil"          
#> 
#> 
#> $patterns[[17]]
#> $patterns[[17]]$original
#> [1] "Viridans streptococci"
#> 
#> $patterns[[17]]$patterns
#> [1] "viridans streptococci" "^v\\.\\s"              "streptococci"         
#> [4] "^v\\.\\s*streptococci"
#> 
#> 
#> $patterns[[18]]
#> $patterns[[18]]$original
#> [1] "Aerococcus spp."
#> 
#> $patterns[[18]]$patterns
#> [1] "aerococcus spp." "aerococcus"      "^a\\.\\s"        "^a\\.\\s*spp."  
#> 
#> 
#> $patterns[[19]]
#> $patterns[[19]]$original
#> [1] "Kocuria spp."
#> 
#> $patterns[[19]]$patterns
#> [1] "kocuria spp."  "kocuria"       "^k\\.\\s"      "^k\\.\\s*spp."
#> 
#> 
#> $patterns[[20]]
#> $patterns[[20]]$original
#> [1] "Dermacoccus spp."
#> 
#> $patterns[[20]]$patterns
#> [1] "dermacoccus spp." "dermacoccus"      "^d\\.\\s"         "^d\\.\\s*spp."   
#> 
#> 
#> $patterns[[21]]
#> $patterns[[21]]$original
#> [1] "Rothia spp."
#> 
#> $patterns[[21]]$patterns
#> [1] "rothia spp."   "rothia"        "^r\\.\\s"      "^r\\.\\s*spp."
#> 
#> 
#> $patterns[[22]]
#> $patterns[[22]]$original
#> [1] "Same as blood (skin flora contamination during needle insertion)"
#> 
#> $patterns[[22]]$patterns
#> [1] "same as blood (skin flora contamination during needle insertion)"
#> [2] "^s\\.\\s"                                                        
#> [3] "as"                                                              
#> [4] "^s\\.\\s*as"                                                     
#> [5] "same as"                                                         
#> 
#> 
#> $patterns[[23]]
#> $patterns[[23]]$original
#> [1] "Viridans group streptococci"
#> 
#> $patterns[[23]]$patterns
#> [1] "viridans group streptococci" "^v\\.\\s"                   
#> [3] "group"                       "^v\\.\\s*group"             
#> [5] "viridans group"             
#> 
#> 
#> $patterns[[24]]
#> $patterns[[24]]$original
#> [1] "Non-pathogenic Neisseria spp."
#> 
#> $patterns[[24]]$patterns
#> [1] "non-pathogenic neisseria spp." "^n\\.\\s"                     
#> [3] "^n\\.\\s*neisseria"            "non-pathogenic neisseria"     
#> 
#> 
#> $patterns[[25]]
#> $patterns[[25]]$original
#> [1] "Haemophilus parainfluenzae"
#> 
#> $patterns[[25]]$patterns
#> [1] "haemophilus parainfluenzae" "^h\\.\\s"                  
#> [3] "parainfluenzae"             "^h\\.\\s*parainfluenzae"   
#> 
#> 
#> $patterns[[26]]
#> $patterns[[26]]$original
#> [1] "Anaerobic oral flora (Prevotella, Porphyromonas, Fusobacterium)"
#> 
#> $patterns[[26]]$patterns
#> [1] "anaerobic oral flora (prevotella, porphyromonas, fusobacterium)"
#> [2] "^a\\.\\s"                                                       
#> [3] "oral"                                                           
#> [4] "^a\\.\\s*oral"                                                  
#> [5] "anaerobic oral"                                                 
#> 
#> 
#> $patterns[[27]]
#> $patterns[[27]]$original
#> [1] "Actinomyces spp."
#> 
#> $patterns[[27]]$patterns
#> [1] "actinomyces spp." "actinomyces"      "^a\\.\\s"         "^a\\.\\s*spp."   
#> 
#> 
#> $patterns[[28]]
#> $patterns[[28]]$original
#> [1] "Same as blood"
#> 
#> $patterns[[28]]$patterns
#> [1] "same as blood" "^s\\.\\s"      "as"            "^s\\.\\s*as"  
#> [5] "same as"      
#> 
#> 
#> $patterns[[29]]
#> $patterns[[29]]$original
#> [1] "Corynebacterium species (diphtheroids)"
#> 
#> $patterns[[29]]$patterns
#> [1] "corynebacterium species (diphtheroids)"
#> [2] "corynebacterium"                       
#> [3] "^c\\.\\s"                              
#> [4] "^c\\.\\s*species"                      
#> [5] "corynebacterium species"               
#> [6] "coryno"                                
#> [7] "diphtheroids?"                         
#> 
#> 
#> $patterns[[30]]
#> $patterns[[30]]$original
#> [1] "Cutibacterium (Propionibacterium) acnes"
#> 
#> $patterns[[30]]$patterns
#> [1] "cutibacterium (propionibacterium) acnes"
#> [2] "^c\\.\\s"                               
#> [3] "(propionibacterium)"                    
#> [4] "^c\\.\\s*(propionibacterium)"           
#> [5] "cutibacterium (propionibacterium)"      
#> [6] "propioni"                               
#> [7] "cuti"                                   
#> [8] "p\\.?\\s*acnes"                         
#> 
#> 
#> $patterns[[31]]
#> $patterns[[31]]$original
#> [1] "Bacillus species"
#> 
#> $patterns[[31]]$patterns
#> [1] "bacillus species" "bacillus"         "^b\\.\\s"         "^b\\.\\s*species"
#> [5] "bacil"           
#> 
#> 
#> $patterns[[32]]
#> $patterns[[32]]$original
#> [1] "Bacteroides-Prevotella group"
#> 
#> $patterns[[32]]$patterns
#> [1] "bacteroides-prevotella group" "^b\\.\\s"                    
#> [3] "^b\\.\\s*group"              
#> 
#> 
#> $patterns[[33]]
#> $patterns[[33]]$original
#> [1] "Eubacterium"
#> 
#> $patterns[[33]]$patterns
#> [1] "eubacterium" "^e\\.\\s"   
#> 
#> 
#> $patterns[[34]]
#> $patterns[[34]]$original
#> [1] "Lactobacillus"
#> 
#> $patterns[[34]]$patterns
#> [1] "lactobacillus" "^l\\.\\s"      "bacil"         "lacto"        
#> 
#> 
#> $patterns[[35]]
#> $patterns[[35]]$original
#> [1] "Faecalibacterium"
#> 
#> $patterns[[35]]$patterns
#> [1] "faecalibacterium" "^f\\.\\s"        
#> 
#> 
#> $patterns[[36]]
#> $patterns[[36]]$original
#> [1] "Lactobacillus species"
#> 
#> $patterns[[36]]$patterns
#> [1] "lactobacillus species" "lactobacillus"         "^l\\.\\s"             
#> [4] "^l\\.\\s*species"      "bacil"                 "lacto"                
#> 
#> 
#> $patterns[[37]]
#> $patterns[[37]]$original
#> [1] "Coagulase-negative Staphylococcus (Staphylococcus epidermidis, Staphylococcus saprophyticus)"
#> 
#> $patterns[[37]]$patterns
#> [1] "coagulase-negative staphylococcus (staphylococcus epidermidis, staphylococcus saprophyticus)"
#> [2] "^c\\.\\s"                                                                                    
#> [3] "^c\\.\\s*staphylococcus"                                                                     
#> [4] "coagulase-negative staphylococcus"                                                           
#> 
#> 
#> $patterns[[38]]
#> $patterns[[38]]$original
#> [1] "Coagulase-negative Staphylococcus (Staphylococcus epidermidis, Staphylococcus hominis, Staphylococcus haemolyticus, Staphylococcus capitis)"
#> 
#> $patterns[[38]]$patterns
#> [1] "coagulase-negative staphylococcus (staphylococcus epidermidis, staphylococcus hominis, staphylococcus haemolyticus, staphylococcus capitis)"
#> [2] "^c\\.\\s"                                                                                                                                   
#> [3] "^c\\.\\s*staphylococcus"                                                                                                                    
#> [4] "coagulase-negative staphylococcus"                                                                                                          
#> 
#> 
#> $patterns[[39]]
#> $patterns[[39]]$original
#> [1] "Neisseria spp. (non-pathogenic)"
#> 
#> $patterns[[39]]$patterns
#> [1] "neisseria spp. (non-pathogenic)" "neisseria"                      
#> [3] "^n\\.\\s"                        "^n\\.\\s*spp."                  
#> [5] "neisseria spp."                 
#> 
#> 
#> $patterns[[40]]
#> $patterns[[40]]$original
#> [1] "Haemophilus spp."
#> 
#> $patterns[[40]]$patterns
#> [1] "haemophilus spp." "haemophilus"      "^h\\.\\s"         "^h\\.\\s*spp."   
#> 
#> 
#> $patterns[[41]]
#> $patterns[[41]]$original
#> [1] "Moraxella catarrhalis"
#> 
#> $patterns[[41]]$patterns
#> [1] "moraxella catarrhalis" "^m\\.\\s"              "catarrhalis"          
#> [4] "^m\\.\\s*catarrhalis" 
#> 
#> 
#> $patterns[[42]]
#> $patterns[[42]]$original
#> [1] "No typical commensals (growth usually clinically significant)"
#> 
#> $patterns[[42]]$patterns
#> [1] "no typical commensals (growth usually clinically significant)"
#> [2] "^n\\.\\s"                                                     
#> [3] "typical"                                                      
#> [4] "^n\\.\\s*typical"                                             
#> [5] "no typical"                                                   
#> 
#> 
#> 
```
