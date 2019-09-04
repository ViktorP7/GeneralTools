# 19-08-19 Ported over R script to plot map of sequenced isolates

using Plots, DataFrames, Luxor
gr()

# Path to polygons
pathCoords = "C:/Users/UCD/Documents/Lab/Cork MAP/PolygonCoords/"

# Create an array of counties
counties = ["Antrim","Armagh","Carlow","Cavan","Clare",
              "Cork","Donegal","Down",
              "Dublin","Fermanagh","Galway","Kerry","Kildare",
              "Kilkenny","Laois","Leitrim",
              "Limerick","Derry","Longford","Louth","Mayo",
              "Meath","Monaghan","Offaly",
              "Roscommon","Sligo","Tipperary","Tyrone",
              "Waterford","Westmeath","Wexford","Wicklow"]

# Define function to get polygons
function getPolyCoords(counties, path)

  # Define dictionary to store coords
  coords = Dict()

  # Loop thru each county and pull out relevant coords into list
  for index in 1:length(counties)

    fileName = string(path, "PolygonCoords_", counties[index], ".txt")

    current = readtable(fileName, header = true, separator = '\t')

    coords[[counties[index]]] = readtable(fileName, header = true, separator = '\t')
  end

  return coords

end

# Get the polygons
polygonCoords = getPolyCoords(counties, pathCoords)

# Define function to get map limits
function mapLimits(polygonCoords, counties)

  # Initialise vectors to store the mins and maxes of each county
  minX = Array{Float64}(undef, length(counties))
  maxX = Array{Float64}(undef, length(counties))
  minY = Array{Float64}(undef, length(counties))
  maxY = Array{Float64}(undef, length(counties))

  for index in 1:length(counties)

    # Describe returns a df of summary stats, including min max
    minX[index] = describe(polygonCoords[[counties[index]]])[1,3]
    maxX[index] = describe(polygonCoords[[counties[index]]])[1,5]
    minY[index] = describe(polygonCoords[[counties[index]]])[2,3]
    maxY[index] = describe(polygonCoords[[counties[index]]])[2,5]
  end

  # Create object ranges with overall mins and maxes
  # use minimum and maximum cos min/max suck!
  ranges = [minimum(minX), maximum(maxX), minimum(minY), maximum(maxY)]

  return(ranges)

end

# Get the map limits
limits = mapLimits(polygonCoords, counties)

# Plot a blank window using limits
plot(x=NaN, y=NaN, xlims=(limits[1], limits[2]),
    ylims=(limits[3], limits[4]), grid = false,
    framestyle = :none)

# Define function to plot polygons on the blank plot
function plotPolygons(polygonCoords, counties)

  for index in 1:length(counties)

    #Plot without initialising a new window
    poly(polygonCoords[[counties[index]]])

  end
end

plotPolygons(polygonCoords, counties)

savefig(ireland, "ireland.pdf")
