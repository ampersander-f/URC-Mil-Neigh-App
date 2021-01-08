# URC-Mil-Neigh-App

## A brief overview of the project

My lab is interested in nanoparticles, specifically collections of nanoparticles spread in a single layer. This is often called a "2D material", and the ones I studied are known to be unusually strong--imagine a single sheet of tissue paper being able to lift a bowling ball without tearing. These materials would make excellent candidates for technologies such as ultra-thin microelectronics for solar cells, or perhaps embedding chemical catalysts (e.g. the metals inside your car's catalytic converter) among the nanoparticles then having the 2D material fold itself into a useful shape (e.g. high surface-area forms, which increase catalytic efficiency).

Unfortunately, we still don't know *how* these 2D materials work. The 2D film is made by depositing the nanoparticles on the surface of water, and then using barriers to confine the monolayer to a rectangular shape. The nanoparticles are naturally attracted to one another, but we found that when we use these barriers to slowly compress the monolayer, the strength of the material increases exponentially. We measure this strength using Wilhelmy plates, which are small metal rectangles hung upright via string so that their bottom edge are in the monolayer. As the strength of the monolayer increases, it pulls the Wilhelmy plates down into the water; the tugging on the string is recorded and when integrated over monolayer area, is the strength. This strength eventually reaches a critical point, and the monolayer finally collapses--imagine a piece of paper wrinkling when you push two opposite edges together. The strength decreases to low values nearly instantaneously after the collapse, so knowing when collapse will happen will be vital in developing technologies down the road.

I worked on implementing these experiments as simulations. The first goal was to demonstrate that the simulations replicated experimental values, so the strength curves between my simulations and the experiments had to match (and they do). The ongoing goal is to understand the mechanism of collapse in the monolayer using only the thermodynamic variables (pressure, free energy, entropy, temperature, etc.), and to use these relationships to control/direct collapse in useful ways. 

## How does it run?

I designed this code to be as streamlined as possible, so that it is possible to run very different simulations but only having to change one file to do so. The 'simulation' is a *class* rather than a function, allowing for easy parallel computing as each simulation is an instance and does not depend on shared resources. This idea was extended to the analysis and graphical methods, and the simulation instance will automatically pass on the correct files and constants to the analysis and ~prevents me from using the wrong set of values for two weeks~ reduces human error between runs.

The simulation is specified in the 'runner' file, on here as `runner_CAURS.py`. All of the necessary python libraries are imported, including those necessary for the classes instantiated later (i.e. each library has global access only so there aren't 6 versions of numpy open at the same time). `DASH` is the molecular dynamics (MD) engine and we specify the build + location. The constants and parameters (e.g. potential function, box size, compression rate, etc) are specified, and last are the command codes (they designate number of time steps, whether to compress or not, and the 'fixes' for DASH, which are C++ commands in a python wrapper). A simulation instance is created, and the the run() command runs the simulation. 

The simulation class (`simulator.py`) mostly functions to automate populating the DASH functions for the specified parameters, which makes it much easier to run these simulations--no need to rewrite the simulations or have a million copy-pasted versions with minor parameter changes. Since DASH requires CUDA programs that specify each potential function (included in `/DASH_RL_potentials`), I ended up making simulation subclasses for each potential as we needed it (`SimulatorRL.py`, to allow for different numbers of instance variables across simulations.

Once the simulation is done, an .xyz file is spit out, and it specifies the xyz location of each particle at each timestep of the simulation. This filename is chosen ahead of time, and so the runner automatically finds the new .xyz file. 

We next specify the `Grapher` instance we would like to use (allows us to graph the same data multiple ways without re-running the analysis or saving the otherwise prohibitively large data set). 

The `PressureSensor` instance is then created (and requires a `Grapher` object), which defines a region in the simulation where we pretend there is a Wilhelmy plate (as with the `simulation` instance, we have subclasses, `PressureSensorRL.py` and `PressureSensor2Body.py`, which define different potential functions and their integrated forms). 

At last the `analyzer.py` instance is created, and the `PressureSensor` and `Grapher` passed to it. The .xyz file is opened and read, and the main `analyzer` function (`perform_analysis()`) is called.

Finally, when done, the .xyz file is closed to preserve the data. The .xyz file and all graphs, charts, and data files created by the `Grapher`s are tidily stored in the same directory, and the full simulation/analysis/visualization is complete.


## Some Background Information and Technical Details

Two-dimensional materials, often called single-layer materials, can be formed using nanoparticles. These materials may be tailored to suit a wide range of needs, with applications that span photovoltaics, microelectronics, and surface chemistry/catalysis (e.g. water purification). My lab was interested in 2D materials with unusually high tensile strengths (resistance to breaking under mechanical stresses like shear force), specifically thiol-ligated gold nanoparticle monolayers. These nanoparticles have a solid gold core, made of a few hundred gold atoms, which were then covered in hydrocarbon 'hairs' attached via sulfur atoms ("thiol-ligated"). The tensile strength of these layers is attributed to a velcro-like effect of the hairs with neighboring particles.  

Our goal for this project was to develop an empirical theory of the properties of these 2D films, so that it would be possible to precisely manipulate the films (potential direct applications included dispersing a chemical catalyst throughout the film, then transforming it into useful shapes like cylinders or high-surface-area forms, e.g.). One group on this project did wet lab work, where the nanoparticles are deposited on the surface of water, then compressed in different ways to force the films to self-assemble. Key to this was the use of Wilhelmy plates, which are very thin metal rectangles suspended by fishing wire so that the bottom of the plates are submerged. By orienting one perpendicular and the other parallel to the direction of compression, we may observe the change in surface pressure as the 'velcro effect' draws the nanoparticles closer together. The surface pressures are then integrated over the surface area of the film to find two useful measures of tensile strength: the Strain Modulus (resistance to stretching) and the Shear Modulus (resistance to twisting). Eventually these forces grow to the point of collapse, and the monolayer begins 'folding' into a bilayer.  

My project aimed to replicate these experimental results using simulations. A Pritzker Institute for Molecular Engineering graduate student built a molecular dynamics simulation engine that we adapted for use with our supercomputing cluster ([link to github page for the engine](https://github.com/dreid1991/md_engine)). By doing many first-principles MD simulations of two, three, and seven nanoparticles drawing closer together (i.e. simulating the nanoparticles by monitoring each carbon atom on the 'hairs' and measuring their attractive energies), we were able to develop empirical equations for the potential function--in other words, we had an equation that if we knew how far apart two nanoparticles were, we could calculate the strength of the velcro effect between them. The engine and the potential functions were then used for Monte Carlo methods to model how the monolayer changes under compression.  

The strain and stress moduli curves from the simulations were eventually found to match the experimental ones (using 7 nanoparticles to find the potential function yielded results on the same order as the results from experiment, whereas using only 2 or 3 only matched the general shape of the curves but overestimated their strength). From this, we were confident that the simulations could provide useful mechanistic insight to the collapse of the monolayer at the moments of high stress (and it's much cheaper to run simulations than to buy gold nanoparticles for each tiny change you want to make). Additionally, the simulations verified that the equilibrium assumptions of the experimental conditions held, which were important in many of its applications; the most important was that the collapse patterns that we observed, which formed diagonal striations across the layer, were a result of only interparticle interactions and totally independent of the compressions we were doing (i.e. the compression of the monolayer was slow enough that the system was at equilibrium at all times).   

In its current form, the project is focused on evaluating the effect of each of the key variables on monolayer collapse: 
* height distributions (do random fluctuations in particle height determine where collapse occurs?)
* packing and lattice effects (transition state theory, or in essence, do misalignments in packing density determine collapse location?)
* local pressure fluctuations (what does the surface pressure look like at points of interest, not just wherever a Wilhelmy plate is?)
* local free energy fluctuations (how much energy from attractions between nanoparticles is released when the bilayer collapses? is there an energy threshhold that determines when collapse happens? are the striations we observe a consequence of minimizing the energy of the bilayer?)
