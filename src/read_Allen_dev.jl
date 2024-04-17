using CSV
using DataFrames
using GeoDataFrames
using GeoJSON

const DATA_DIR = "C:/Users/rlippman/Documents/development/ADRIA_data/spatial_datasets/GBR-Bathy10m"
const ALLEN_DIR = joinpath(DATA_DIR, "..", "AllenAtlas_GBR-20231118074407")

const criteria =  ["Benthic-Map", 
                    "Geomorphic-Map", "Reef-Extent"]          #readdir(DATA_DIR)    # Criteria = all files/folders in the directory

# Function to read and print GeoJSON file
function read_and_print_geojson(filepath, criteria)
    for c in criteria        
    directory_path = joinpath(filepath, c)                                                  # GeoJSON file containing directory
    geojson_files = filter(file -> endswith(file, ".geojson"), readdir(directory_path))     # Filter GeoJSON files
    geojson_files_path = [joinpath(directory_path, geojson_file) for geojson_file in geojson_files]
    if !isempty(geojson_files)                                                                  # if there are GeoJSON files
            for i in 1:length(geojson_files_path)                                                    # Loop through GeoJSON files
                df = GeoDataFrames.read(geojson_files_path[i])
                # df = DataFrame(fc)
                # display(df)
                println("$(c)\tSize: ", size(df), "\tNames: ", names(df))
                println(first(df, 5))

            # Write data to a CSV file
            try 
                CSV.write(joinpath(ALLEN_DIR, "outputs", "df_data_$(c).csv"), df)#, delim=',', writeheader=false)
                println("Data written to CSV file.")
            catch err
                println("An error occurred while writing data to CSV: ", err)
            end 
        end

        else
            println("No GeoJSON files found in the directory.")
        end
    end
end

read_and_print_geojson(ALLEN_DIR, criteria)

