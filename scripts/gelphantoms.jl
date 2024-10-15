const gelphantoms = (
    T₁ = Dict(
        #              1   2   3   4   5   6   7   8   9   10  11   12   13   14   15   16   17   18
        292 => Float64[200,299,296,463,450,444,604,596,754,745, 903,1448, 966,1034,1160,1276,1262,1415],
        296 => Float64[223,334,331,516,502,496,674,666,841,831,1007,1615,1078,1153,1293,1422,1407,1576],
        300 => Float64[249,372,368,574,559,552,750,741,935,924,1120,1796,1199,1282,1437,1581,1563,1750]
    ),
    T₂ = Dict(
        292 => Float64[52 ,73 ,113,53 ,94 ,154,95 ,136,116,157,137 ,390 ,224 ,167 ,214 ,204 ,184 ,174 ],
        296 => Float64[50 ,70 ,111,49 ,89 ,151,89 ,129,109,149,128 ,373 ,212 ,156 ,201 ,190 ,171 ,161 ],
        300 => Float64[48 ,67 ,110,46 ,85 ,148,84 ,124,103,142,121 ,359 ,203 ,148 ,191 ,180 ,161 ,151 ]
        )
)

# interpolate to 294K
gelphantoms.T₁[294] = (gelphantoms.T₁[292] + gelphantoms.T₁[296])/2
gelphantoms.T₂[294] = (gelphantoms.T₂[292] + gelphantoms.T₂[296])/2
