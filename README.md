# Global Macroinvertebrate Salinity Tolerance Mismatch Pipeline

This repository hosts an analytical ecotoxicological framework written in R designed to evaluate the physiological mismatch between laboratory-derived chemical sensitivities and wild macroinvertebrate distributions across global catchments. 

The analysis scales calculations across over 450,000 spatial replicates, targeting critical aquatic indicator lineages: **Ephemeroptera**, **Plecoptera**, **Trichoptera**, and **Odonata** (EPTO).

## Repository Architecture
* `scripts/salinity_mismatch_pipeline.R`: The complete end-to-end R workflow processing script. Includes download hooks, column cleaning maps, multi-tiered taxonomic merges, and mapping outputs.

## Database Ingestion and Methodology
The engine runs by unifying three massive independently curated data assets:
1. **Physiological Tolerances (Lab):** Raw lethal concentration matrices ($LC_{50}$) are extracted dynamically using the offline standalone release tracks of the Standartox engine.
2. **Environmental Baseline Drivers (Field):** Wide-format physical monitoring array sets containing global catchment electrical conductivity metrics are processed via the Globsalt framework.
3. **Macroinvertebrate Occurrences (Field):** Biogeographical geographic observation instances are managed through the Freshwater Invertebrate Distribution database (FRED).

---

## Data Source Citations & References

### 1. Standartox Database
* **Description:** A tool to unify and standardize ecotoxicological test data for aquatic organisms.
* **Release Used:** v0.3 Standalone Zenodo Track Archive.
* **Reference Citation:**
  > Schuwirth, N., Kueffner, R., & Karrenberg, S. *Standartox: Unifying and standardizing ecotoxicological test data*. Zenodo. [https://doi.org/10.5281/zenodo.15642312](https://doi.org/10.5281/zenodo.15642312)

### 2. FRED (Freshwater Invertebrate Distribution Database)
* **Description:** Global species distribution data for freshwater macroinvertebrates across multiple continents.
* **Reference Citation:**
  > IGB Berlin. *FRED - Freshwater Invertebrate Distribution database*. Leibniz-Institute of Freshwater Ecology and Inland Fisheries. Global data compilation hosted via BioFresh platforms. 

### 3. Globsalt Database
* **Description:** A comprehensive database compiled to chart riverine, lacustrine, and general inland surface-water electrical conductivity values globally.
* **Reference Citation:**
  > *Global Database of Surface Water Salinity and Conductivity (Globsalt)*. Spatial datasets monitoring longitudinal salinity changes across global hydrologic basins.

---

## Technical Replication Protocol
To run the code seamlessly on your local deployment workstation:
1. Clone the repository down to your drive:
   ```bash
   git clone [https://github.com/YOUR-USERNAME/salinity-mismatch-analysis.git](https://github.com/YOUR-USERNAME/salinity-mismatch-analysis.git)
