## modERN data download utils

From the parent directory, run:
util/ChIP_Manifest.downloader.sh

It reads the Excel File in the parent directory (DATA) https://github.com/meekrob/ELT-2-ChIP-revision/blob/388f8f68a8c81f9cd0d3ad5d2169abe6c70126f3/DATA/ChIP_Source_Data_Manifest.xlsx
using exportXLsheet.RScript, and downloads the data files using wget.
See the ChIP_Source_Data_Manifest.xlsx for the downloaded files and their source location.
