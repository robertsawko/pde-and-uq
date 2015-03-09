# Suggestion on how to do mesh refinement with a save to every directory.
# To be used with parun -b 
nList[0].nn=11
so.name = p.name+"_"+`nList[0].nn`+"_"
start
nList[0].nn=21
so.name = p.name+"_"+`nList[0].nn`+"_"
start
nList[0].nn=41
so.name = p.name+"_"+`nList[0].nn`+"_"
start
nList[0].nn=81
so.name = p.name+"_"+`nList[0].nn`+"_"
start
quit
