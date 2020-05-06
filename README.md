# How to compile the project

`make`

The Makefile allows us to choose the MaxSAT solver to our liking.  To change the default solver ([TT-Open-WBO-INC]) change the variables: SUPERSOLVERNAME and SUPERSOLVERNAMEID.  

# How to run the project

`./timetabler data/PESP/set-01/R1L1.xml -opt-time=2  [solver options]`


The solver option depend on the solver used. Please read the solver documentation.

The following options are available for all solvers:

## Timetabler Options
### Entry time variables: 0 - For all section and all time; 1 - For all time; 2 - Smart time
```-opt-time= <int32>  [   0 ..    2]      (default: 2)```

# Dependencies
   
 [TT-Open-WBO-INC](https://drive.google.com/file/d/140d8jDHZHo5d7WuoNpLqZXmHasgYkH38/view) solver
 [Loandra](https://maxsat-evaluations.github.io/2019/mse19-solver-src/incomplete/Loandra.zip)
 [Open-WBO-Inc](https://github.com/sbjoshi/Open-WBO-Inc)
 [LinSBPS](https://maxsat-evaluations.github.io/2019/mse19-solver-src/incomplete/LinSBPS2018.zip)
 [SATLike](https://maxsat-evaluations.github.io/2019/mse19-solver-src/incomplete/SATLike3.0-c-weighted.zip)
 c++ compiler.
 [Rapid JSON Parser](https://rapidjson.org/)
  
# Data Sets
   
[PESP benckmark](http://num.math.uni-goettingen.de/~m.goerigk/pesplib/)
[SBB benchmark](https://github.com/potassco/train-scheduling-with-hybrid-asp)
