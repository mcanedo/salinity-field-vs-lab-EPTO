# ==============================================================================
# GLOBAL EPTO SALINITY CONGRUENCE PROJECT (HIGH-PRECISION SPATIAL PIPELINE)
# Framework: Lab Tolerance (Standartox API) vs. Field Realization (FRED + GlobSalt CSV)
# Spatial Constraints: Strict 2 km maximum distance threshold for stream reaches
# Standardized to Practical Salinity Scale (PSS) using PSS-78
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. DEPENDENCIES & PACKAGE LOADING
# ------------------------------------------------------------------------------
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
if (!requireNamespace("standartox", quietly = TRUE)) remotes::install_github("andschar/standartox")

library(standartox)  # Direct API access to lab ecotox data
library(tidyverse)   # Tidy data wrangling and ggplot2 visualization
library(sf)          # Simple Features for vector spatial operations (Point-to-Point Joins)
library(data.table)  # High-performance processing for large CSVs (GlobSalt)
library(lme4)        # Linear mixed-effects modeling
library(wql)         # Aquatic monitoring tools (contains ec2pss function)

# ------------------------------------------------------------------------------
# 2. HELPER FUNCTION: UNIT STANDARDIZATION
# ------------------------------------------------------------------------------
#' Convert Electrical Conductivity to Practical Salinity Units (PSS)
convert_ec_to_pss <- function(ec_val, input_unit = "uS/cm", temp = 25) {
  ec_ms <- ifelse(input_unit == "uS/cm", ec_val / 1000, ec_val)
  pss_value <- wql::ec2pss(ec = ec_ms, t = temp, p = 0)
  return(pss_value)
}

# ------------------------------------------------------------------------------
# 3. STATIC DATA IMPORT & LOCAL STANDARDIZATION (EXPLICIT SCHEMA MATCHING)
# ------------------------------------------------------------------------------
message("Bypassing decommissioned API server...")
message("Fetching standalone Standartox database components from Zenodo...")

if (!requireNamespace("fst", quietly = TRUE)) install.packages("fst")
library(fst)

# Static source links for the official Standartox v0.3 data components on Zenodo
test_url  <- "https://zenodo.org/records/15642312/files/test_fin.fst?download=1"
taxa_url  <- "https://zenodo.org/records/15642312/files/taxa.fst?download=1"
phch_url  <- "https://zenodo.org/records/15642312/files/phch.fst?download=1" 

temp_test <- tempfile(fileext = ".fst")
temp_taxa <- tempfile(fileext = ".fst")
temp_phch <- tempfile(fileext = ".fst")

message("Downloading toxicity arrays, taxonomy tables, and chemical properties...")
download.file(test_url, destfile = temp_test, mode = "wb", quiet = TRUE)
download.file(taxa_url, destfile = temp_taxa, mode = "wb", quiet = TRUE)
download.file(phch_url, destfile = temp_phch, mode = "wb", quiet = TRUE)

# Read directly into highly optimized data.tables
stx_tests <- read_fst(temp_test) %>% as.data.table()
stx_taxa  <- read_fst(temp_taxa)  %>% as.data.table()
stx_phch  <- read_fst(temp_phch)  %>% as.data.table()

message("Executing explicit relational table joins...")

# Step 1: Bind biological taxonomic hierarchies via Taxon Link ID (tl_id)
merged_meta <- merge(stx_tests, stx_taxa, by = "tl_id", all.x = TRUE)

# Step 2: Bind chemical names and descriptors via Chemical Link ID (cl_id)
standartox_local <- merge(merged_meta, stx_phch, by = "cl_id", all.x = TRUE)

message("Executing targeted filtering and salinity normalization...")
lab_tolerance <- standartox_local %>%
  # Filter using verified schema categories and lowercase strings
  filter(
    endpoint_group == "XX50",
    tolower(effect) == "mortality",
    duration >= 0 & duration <= 96
  ) %>%
  # Explicitly map schema designations to match your downstream script
  rename(
    order   = tax_order,
    family  = tax_family,
    genus   = tax_genus,
    species = tax_taxon,
    unit    = concentration_unit
  ) %>%
  # Isolate diagnostic target parameters
  filter(order %in% c("Ephemeroptera", "Plecoptera", "Trichoptera", "Odonata")) %>%
  filter(str_detect(tolower(cname), "sodium chloride|salinity|sea water|salt")) %>%
  # Standardize units using a case-insensitive match configuration
  mutate(unit_lower = tolower(unit)) %>%
  filter(unit_lower %in% c("us/cm", "ms/cm", "g/l", "mg/l")) %>%
  mutate(
    lab_salinity_pss = case_when(
      unit_lower == "us/cm" ~ convert_ec_to_pss(concentration, input_unit = "uS/cm", temp = 25),
      unit_lower == "ms/cm" ~ convert_ec_to_pss(concentration, input_unit = "mS/cm", temp = 25),
      unit_lower == "g/l"   ~ concentration,
      unit_lower == "mg/l"  ~ concentration / 1000,
      TRUE                  ~ NA_real_
    )
  ) %>%
  filter(!is.na(lab_salinity_pss)) %>%
  # Aggregate raw toxicity values down to a single species-level median threshold
  group_by(order, family, genus, species) %>%
  summarise(
    median_lab_LC50_pss = median(lab_salinity_pss, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(species) & species != "")

# Clean up local system disk storage footprints
unlink(c(temp_test, temp_taxa, temp_phch))
message(paste("Success! Locally synchronized and filtered", nrow(lab_tolerance), "lab-tested EPTO species."))

# ------------------------------------------------------------------------------
# 4. ENVIRONMENTAL LAYER PROCESSING (GLOBALSALT .CSV)
# ------------------------------------------------------------------------------
message("Loading and aggregating GlobSalt CSV file...")
globsalt_raw <- fread("GlobSalt_v2.0.csv")

globsalt_aggregated <- globsalt_raw %>%
  # 1. Drop rows missing vital coordinates or conductivity metrics
  filter(!is.na(Conductivity) & !is.na(x) & !is.na(y)) %>%
  
  # 2. Group by the actual schema names for station and coordinates
  group_by(Station_ID, y, x) %>% 
  
  # 3. Calculate the 95th percentile max EC per monitoring location
  summarise(
    station_max_ec_us = quantile(Conductivity, 0.95, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  
  # 4. Rename columns to match your downstream script requirements
  rename(
    station_id = Station_ID,
    latitude   = y,
    longitude  = x
  ) %>%
  
  # 5. Convert the aggregated EC values into practical salinity units (PSS)
  mutate(
    station_salinity_pss = convert_ec_to_pss(station_max_ec_us, "uS/cm")
  )

message("Converting GlobSalt stations into spatial vector data...")
globsalt_sf <- st_as_sf(globsalt_aggregated, coords = c("longitude", "latitude"), crs = 4326)

# ------------------------------------------------------------------------------
# 5. BIOLOGICAL DISTRIBUTIONS (FRED) & STRICT NEAREST NEIGHBOR SPATIAL JOIN
# ------------------------------------------------------------------------------
message("Loading FRED field occurrence points...")
fred_raw <- fread("Global_EPTO_Database.csv")

fred_epto <- fred_raw %>%
  # 1. Filter using the exact lowercase columns from the database structure
  filter(
    order %in% c("Ephemeroptera", "Plecoptera", "Trichoptera", "Odonata"),
    !is.na(longitude) & !is.na(latitude)
  ) %>%
  
  # 2. Select columns and map them back to your script's preferred capitalization
  select(
    ID        = occurrence_ID,
    Order     = order,
    Family    = family,
    Genus     = genus,
    Species   = species,
    Longitude = longitude,
    Latitude  = latitude
  )

fred_sf <- st_as_sf(fred_epto, coords = c("Longitude", "Latitude"), crs = 4326)

message("Performing Spatial Nearest-Neighbor Join (Robinson metric projection)...")
fred_metric <- st_transform(fred_sf, crs = "ESRI:54030")
globsalt_metric <- st_transform(globsalt_sf, crs = "ESRI:54030")

# Find index of the physically closest water station for every biological observation
nearest_station_indices <- st_nearest_feature(fred_metric, globsalt_metric)

# Pull matching records out, calculate true distances, and filter strictly
fred_joined <- fred_epto %>%
  mutate(
    matched_station_id  = globsalt_aggregated$station_id[nearest_station_indices],
    field_salinity_pss  = globsalt_aggregated$station_salinity_pss[nearest_station_indices],
    spatial_distance_m  = as.numeric(st_distance(fred_metric, globsalt_metric[nearest_station_indices, ], by_element = TRUE))
  ) %>%
  # CRITICAL UPDATE: Enforce maximum 2 km (2000 meters) distance threshold
  filter(spatial_distance_m <= 2000)

# Calculate the realized environmental niche threshold per species based on high-precision data
field_tolerance <- fred_joined %>%
  group_by(Order, Family, Genus, Species) %>%
  summarise(
    n_field_obs = n(),
    max_field_sal_pss = quantile(field_salinity_pss, 0.95, na.rm = TRUE),
    mean_lat = mean(Latitude),
    mean_lon = mean(Longitude),
    .groups = "drop"
  ) %>%
  filter(n_field_obs >= 5) # Retain minimum occurrence threshold to keep niche estimates robust

# ------------------------------------------------------------------------------
# 6. DATA HARMONIZATION & CONGRUENCE MAPPING
# ------------------------------------------------------------------------------
library(dplyr)
library(stringr)

message("Standardizing FRED nomenclature for deep structural matching...")

fred_clean_matches <- fred_joined %>%
  mutate(
    # 1. Convert underscores to standard spacing
    species_clean = str_replace_all(Species, "_", " "),
    # 2. Strip out trailing author names/years (isolate the first two words)
    species_clean = str_extract(species_clean, "^\\S+\\s+\\S+"),
    # 3. Ensure case alignment across all keys
    species_clean = tolower(species_clean),
    Genus_clean   = tolower(Genus),
    Family_clean  = tolower(Family)
  )

# Prepare lab dataset keys for smooth merging
lab_lookup <- lab_tolerance %>%
  mutate(
    species_clean = tolower(species),
    genus_clean   = tolower(genus),
    family_clean  = tolower(family)
  )

# ------------------------------------------------------------------------------
# HIERARCHICAL JOIN SEQUENCE (Preserves all ~1M rows with best available lab data)
# ------------------------------------------------------------------------------
message("Executing multi-tier taxonomic join to preserve spatial replicates...")

# Tier 1: Direct Species-to-Species matches
match_species <- lab_lookup %>% 
  filter(!is.na(species_clean) & species_clean != "") %>%
  select(species_clean, lab_LC50_species = median_lab_LC50_pss)

# Tier 2: Genus-level pool averages (for field observations missing exact species lab tests)
match_genus <- lab_lookup %>%
  group_by(genus_clean) %>%
  summarise(lab_LC50_genus = median(median_lab_LC50_pss, na.rm = TRUE), .groups = "drop")

# Tier 3: Family-level pool averages (the ultimate safety net)
match_family <- lab_lookup %>%
  group_by(family_clean) %>%
  summarise(lab_LC50_family = median(median_lab_LC50_pss, na.rm = TRUE), .groups = "drop")


# Map benchmarks back to the unaggregated 1-million-row database
congruence_expanded <- fred_clean_matches %>%
  left_join(match_species, by = "species_clean") %>%
  left_join(match_genus,   by = c("Genus_clean" = "genus_clean")) %>%
  left_join(match_family,  by = c("Family_clean" = "family_clean")) %>%
  
  # Coalesce: Use Species lab data if available; if not, use Genus; if not, use Family
  mutate(
    assigned_lab_LC50_pss = coalesce(lab_LC50_species, lab_LC50_genus, lab_LC50_family),
    match_resolution = case_when(
      !is.na(lab_LC50_species) ~ "Species Level",
      !is.na(lab_LC50_genus)   ~ "Genus Level",
      !is.na(lab_LC50_family)  ~ "Family Level",
      TRUE                     ~ "No Match"
    )
  ) %>%
  # Drop anything that couldn't be matched at any taxonomic level
  filter(match_resolution != "No Match") %>%
  # Remove temporary lowercase working variables
  select(-species_clean, -Genus_clean, -Family_clean, -lab_LC50_species, -lab_LC50_genus, -lab_LC50_family)

# ------------------------------------------------------------------------------
# DIAGNOSTIC PROFILE
# ------------------------------------------------------------------------------
message("--- Final Congruence Matrix Profile ---")
print(table(congruence_expanded$match_resolution))
message(paste("Total observations retained for spatial analysis:", nrow(congruence_expanded)))

library(dplyr)
library(ggplot2)
library(maps)

# ------------------------------------------------------------------------------
# SETUP: Calculate Mismatch Metrics
# ------------------------------------------------------------------------------
mismatch_db <- congruence_expanded %>%
  mutate(
    mismatch_val = field_salinity_pss - assigned_lab_LC50_pss,
    is_mismatch  = field_salinity_pss > assigned_lab_LC50_pss
  )

# ------------------------------------------------------------------------------
# Q1: Overall Mismatch Count & Rate
# ------------------------------------------------------------------------------
total_cases   <- nrow(mismatch_db)
mismatch_rows <- sum(mismatch_db$is_mismatch)
mismatch_rate <- mean(mismatch_db$is_mismatch) * 100

cat("--- Q1: OVERALL MISMATCH METRICS ---\n")
cat(paste("Total Observations Analyzed:", total_cases, "\n"))
cat(paste("Number of Mismatch Cases (Field > Lab):", mismatch_rows, "\n"))
cat(paste("Overall Mismatch Frequency:", round(mismatch_rate, 2), "%\n\n"))

# ------------------------------------------------------------------------------
# Q2 & Q3: Genus-Level Frequency and Internal Variability
# ------------------------------------------------------------------------------
genus_analysis <- mismatch_db %>%
  group_by(Order, Genus) %>%
  summarise(
    sample_size = n(),
    mismatch_count = sum(is_mismatch),
    mismatch_percent = mean(is_mismatch) * 100,
    # Standard deviation tracks intra-taxa variance/internal variability
    mismatch_sd = sd(mismatch_val, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Filter out genera with fewer than 30 observations to ensure statistical reliability
  filter(sample_size >= 30)

# Top 5 Genera where mismatch happens most frequently (Q2)
top_mismatch_genera <- genus_analysis %>% arrange(desc(mismatch_percent)) %>% head(5)

# Top 5 Genera with the highest internal variability / spatial divergence (Q3)
top_variable_genera <- genus_analysis %>% arrange(desc(mismatch_sd)) %>% head(5)

cat("--- Q2: TOP 5 GENERA WITH HIGHEST MISMATCH FREQUENCY ---\n")
print(top_mismatch_genera)

cat("\n--- Q3: TOP 5 GENERA WITH HIGHEST INTERNAL VARIABILITY (SD) ---\n")
print(top_variable_genera)

# ------------------------------------------------------------------------------
# VISUALIZATION 1: Genus Variability Boxplot
# ------------------------------------------------------------------------------
# Pick top 15 most well-sampled genera to visualize intra-taxa variability cleanly
top_15_genera <- genus_analysis %>% arrange(desc(sample_size)) %>% head(15) %>% pull(Genus)

p1 <- ggplot(mismatch_db %>% filter(Genus %in% top_15_genera), 
             aes(x = reorder(Genus, mismatch_val, FUN = median), y = mismatch_val, fill = Order)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 0.5, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 1) +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Intra-Taxa Variability in Field vs. Lab Mismatches",
    subtitle = "Values above red line indicate field survival exceeding lab lethal limits",
    x = "Genus",
    y = "Mismatch Magnitude (Field Salinity - Lab LC50 PSS)"
  )
print(p1)

# ------------------------------------------------------------------------------
# Q4: Spatial Mapping of Mismatches
# ------------------------------------------------------------------------------
# Retrieve a base world map outline
world_map <- map_data("world")

p2 <- ggplot() +
  # Draw global map backdrop
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
               fill = "grey90", color = "white", linewidth = 0.1) +
  # Plot spatial data points colored by mismatch magnitude
  geom_point(data = mismatch_db, aes(x = Longitude, y = Latitude, color = mismatch_val), 
             alpha = 0.6, size = 0.8) +
  scale_color_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    name = "Mismatch\nMagnitude\n(ΔS)"
  ) +
  coord_quickmap() +
  theme_minimal() +
  labs(
    title = "Global Spatial Distribution of Eco-Toxicological Mismatches",
    subtitle = "Red dots expose areas where field salinity outpaces known lab lethal tolerances",
    x = "Longitude",
    y = "Latitude"
  )
print(p2)
