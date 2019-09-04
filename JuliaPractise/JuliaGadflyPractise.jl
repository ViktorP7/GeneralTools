
# Plotting with Gadfly tutorial
using Gadfly
plot(y=[1,2,3])

# Import iris dataset from RDatasets
using Gadfly, RDatasets
iris = dataset("datasets", "iris")

# Plot sepal length vs sepal width
p = plot(iris, x=:SepalLength, y=:SepalWidth, Geom.point);

# Type name of object to plot it
p

# Can save as a certain filetype - create file and draw to it
img = SVG("iris_plot.svg", 14cm, 8cm)
draw(img, p)

# Can plot directly without saving plot in variable
plot(iris, x=:SepalLength, y=:SepalWidth)

# Can add lines to plot
plot(iris, x=:SepalLength, y=:SepalWidth, Geom.point, Geom.line)

# Can plot arrays/vectors too
SepLen = iris[:SepalLength]
SepWid = iris[:SepalWidth]
plot(x=SepLen, y=SepWid, Geom.point,
    Guide.xlabel("Sepal Length"), Guide.ylabel("Sepal Width"))

# Can add color to the plot aswell
plot(iris, x=:SepalLength, y=:SepalWidth, Geom.point, color=:Species)

# Can transform scales to get better visual, for example log
set_default_plot_size(21cm ,8cm)
mammals = dataset("MASS", "mammals")
p1 = plot(mammals, x=:Body, y=:Brain, label=:Mammal, Geom.point, Geom.label)
p2 = plot(mammals, x=:Body, y=:Brain, label=:Mammal, Geom.point, Geom.label,
     Scale.x_log10, Scale.y_log10)
hstack(p1, p2)

# Can also transform with _sqrt, _log, _log2, _log10 and _asinh
using Printf
Diamonds = dataset("ggplot2","diamonds")
p3= plot(Diamonds, x=:Price, y=:Carat, Geom.histogram2d(xbincount=25, ybincount=25),
    Scale.x_continuous(format=:engineering) )
p4= plot(Diamonds, x=:Price, y=:Carat, Geom.histogram2d(xbincount=25, ybincount=25),
    Scale.x_continuous(format=:plain),
    Scale.y_sqrt(labels=y->@sprintf("%i", y^2)),
    Scale.color_log10(minvalue=1.0, maxvalue=10^4),
    Guide.yticks(ticks=sqrt.(0:5)) )
hstack(p3, p4)

# Can have discrete scales
mtcars = dataset("datasets","mtcars")
 labeldict = Dict(4=>"four", 6=>"six", 8=>"eight")
p5 = plot(mtcars, x=:Cyl, color=:Cyl, Geom.histogram,
    Scale.x_discrete(levels=[4,6,8]), Scale.color_discrete(levels=[4,6,8]) )
p6 = plot(mtcars, x=:Cyl, color=:Cyl, Geom.histogram,
    Scale.x_discrete(labels=i->labeldict[i], levels=[8,6,4]),
    Scale.color_discrete(levels=[8,6,4]) )
hstack(p5, p6)

# Gadly can guess scales and guides when not supplied
gasoline = dataset("Ecdat", "Gasoline")
plot(gasoline, x=:Year, y=:LGasPCar, color=:Country, Geom.point, Geom.line)

# Gadfly plots long form data (ie, same type in one column)
# But some frames have wide data, as below with male/female in separate cols
births = RDatasets.dataset("HistData", "Arbuthnot")[[:Year, :Males, :Females]]

# Need to stack the male and female columns
stack(births, [:Males, :Females])

# Plot stacked data
births = RDatasets.dataset("HistData", "Arbuthnot")[[:Year, :Males, :Females]] # hide
plot(stack(births, [:Males, :Females]), x=:Year, y=:value, color=:variable,
     Geom.line)

# Can actually avoid stacking and refer as if long form
plot(births, x=:Year, y=Col.value(:Males, :Females),
    color=Col.index(:Males, :Females), Geom.line)
nothing # hide

# Can make multiple plots with same dataset and axis
plot(iris, xgroup="Species", x="SepalLength", y="SepalWidth",
     Geom.subplot_grid(Geom.point))

# Can do multiple layers on the same plot
xdata = sort(iris[:SepalWidth])
ydata = cumsum(xdata)
line = layer(x=xdata, y=ydata, Geom.line, Theme(default_color="red"))
bars = layer(iris, x=:SepalWidth, Geom.bar)
plot(line, bars)

plot(iris,
     layer(x=:SepalLength, y=:SepalWidth),
     layer(x=:PetalLength, y=:PetalWidth, Theme(default_color="red")),
     Guide.xlabel("length"), Guide.ylabel("width"), Guide.title("Iris data"),
     Guide.manual_color_key("",["Sepal","Petal"],
                            [Gadfly.current_theme().default_color,"red"]))
