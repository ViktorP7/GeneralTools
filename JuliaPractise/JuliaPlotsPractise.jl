# Julia Plotly tutorial
using Plots
gr()

x=1:20; y=rand(20,2)

p = plot(x,y,title="A Random Plot");

# Note, plotly plots seem to only save as .htmls, refuse pdf or png
savefig(p, "ARandoPlot.pdf")

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

using DataFrames
# Define function to get polygons
function getPolyCoords(counties, path)

  # Define dictionary to store coords
  coords = Dict()

  # Loop thru each county and pull out relevant coords into list
  for index in 1:length(counties)

    fileName = string(path, "PolygonCoords_", counties[index], ".txt")

    coords[[counties[index]]] = readtable(fileName, header = true, separator = "\t")
  end

  return coords

end

polygonCoords = getPolyCoords(counties, pathCoords)
