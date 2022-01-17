# era5land-utilities
A set of python scripts allowing to download, reggrid, aggregate and concatenate ERA5-land data from the CDS (https://cds.climate.copernicus.eu/cdsapp#!/home).
The four steps prepare files for the main meteorological inputs of the LISVAP evapostranspiration model and the LISFLOOD hydrological model.


# Content:

## 1. main files

* "01_ERA5land_downloader.py": CDS API script to use CDS service to retrieve hourly ERA5* variables and iterate over all months in the specified years. 
* "02_ERA5land_yearly_files.py": creates yearly files of daily ERA5* variables and from aggregated hourly values and creates yealry files from monthly files in the specified years.
* "03_ERA5land_interpolate.py": allows to reggrid yearly files of daily ERA5* variables to another resolution
    * Original grid: 0.1 degrees
    * Objective grid: 1 arc minute (approx 0.016667 deg)
* "04_ERA5land_concatenate.py": concatenates yearly files of daily ERA5* variables into a single file for the whole specified
 
## 2. additional files

* "11_ERA5land_compression_par.py": Uses calculated add offset and scale factor from "02_ERA5land_yearly_files.py" to estimate the optimal parameters for the whole period (1981-2020)
* "12_ERA5land_adapt_time_HPC.py": script that shifts the variable dates by one. Useful for inputs of LISVAP and LISFLOOD.
