using DrWatson
@quickactivate "cmip5"

println(
"""
Currently active project is: $(projectname())

Path of active project: $(projectdir())

Have fun with your new project!
"""
)
