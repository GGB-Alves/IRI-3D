# IRI-3D
Trace a line between two chosen atoms and calculates the sign(λ2)ρ (related to IRI, look for it in Multiwfn guide).

When running a Interactive Region Indicator (IRI) with mMltiwfn you will find yourself with 2 cube files, one for sign(λ2)ρ and the other for the IRI (a surface).
You will use VMD to load the result (on VMD Tk console: source IRIfill.vmd), and controlling your isosurface you can have multiple different results (surfaces).

Well, with this Fortran program you can choose 2 atoms and calulate the values of sign(λ2)ρ between them (in a line) and display "how repulsive it is".
  sign(λ2)ρ is responsible for the color map, thus defining what kind of interaction you have (again, refer to the guide).

The program workd as intended, HOWEVER, the linear approximation approach (LAA) is too simplistic.
I intend to work on two gradual changes:
  i) Add a circular path approximation;
  ii) Make it 3D, using a cylinder (for the line) or torus (for the circular path);

Achieving i) is relatively easy, only require coding the building of circle segment, it still is 1D so the interpolation should work just fine! Of course I'd not abandon LAA, rather I'll just add a option to choose between two approaches, because LAA still has some value (not much though);

Achieving ii) will be quite challenging, I'll need to go from 1D to 3D, instead of values on a line will be values on a volume. I still have to sit down and figure out how to approach it (study math).

I'm no professional coder, nor I have great coding skills, I'm sure my code is not optional (I can't even tell) and it's probably slow (it's pretty slow reading and storing the sign(λ2)ρ values, but it's almost two million points);

About the code: 
  JMol and Povray interactions are completely avoidable if you are sure your XYZ structures are the same in the func1.cub file and your reference file.
  If your Jmol, for some miracle, always loads and display your molecule properly, Povray can be tossed away;
  If your Jmol always loads poorly when called from another routine, you can choose to run it headless, I did not have the time to fix this, I wrote on a huge time constrain;
  Multiwfn always outputs as func1.cub and func2.cub, I just followed this logic and didn't botter adding a name reading part, you just run it directly with func1.cub in the folder and you are good to go =D
  Lastly, I tried to leave as many explanation comments as possible so it's clear what I'm doing (for my own sake in case I come back to it after some months).
