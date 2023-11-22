# BinaryCelestialSim
Generates Simulations of a binary star system by numerically solving twelve first order equations of motion for two stars.

## 3-D Simulations
Here we will show some examples of different cases and types of orbits that were simulated using the binary star system script. 

### Case 1: Sinusoidal Orbits
#### Animated Gif of Orbits
The mass of the stars is the same and the initial velocity vector of star one is V<sub>0</sub> = (0 0.5 0). This results in a very circular orbit where both stars follow the same path.

![](https://github.com/GrahamDrive/BinaryCelestialSim/blob/main/Example%20Simulations/M1=1%20M2=1%20v01%20=%20%5B0%200.5%200%5D/StarOrbit%20Graph%20for%20M1=1%20M2=1.gif?raw=true)

#### Y-Axis Velocity
When we look at y component of the velocity we can see that this type of orbit creates a sinusoidal shape which we like to call a sinusoidal orbit.

![](https://github.com/GrahamDrive/BinaryCelestialSim/blob/main/Example%20Simulations/M1=1%20M2=1%20v01%20=%20%5B0%200.5%200%5D/Y-Velocity%20Graph%20for%20M1=1%20M2=1.png?raw=true)

#### Velocity Magnitudes
Now if we look at the magnitude of the two stars velocity we can see that they are identical they remain relatively stable. The sinusoid that appears in the graphs we suspect are more due to the inaccuracies in ode45 differentiation rather than and would not match experimental data.

| Star One| Second Star|
|---------|------------|
|![](https://github.com/GrahamDrive/BinaryCelestialSim/blob/main/Example%20Simulations/M1=1%20M2=1%20v01%20=%20%5B0%200.5%200%5D/Star%201%20Velocity%20Magnitude%201%20M2=1.png?raw=true)|![](https://github.com/GrahamDrive/BinaryCelestialSim/blob/main/Example%20Simulations/M1=1%20M2=1%20v01%20=%20%5B0%200.5%200%5D/Star%202%20Velocity%20Magnitude%201%20M2=1.png?raw=true)|


### Case 2: Unbalanced System
#### Animated Gif of Orbits 
This second system demonstrates the effect of unbalanced masses on a binary system. The gif is a little fast but we can see how the smaller star M1 has a much larger orbit as it is "tossed around" by its bigger relative. One interesting characteristic of this binary specifically is that the orbits of the two stars actually cross each other far from the center of mass.

![](https://github.com/GrahamDrive/BinaryCelestialSim/blob/main/Example%20Simulations/M1=2%20M2=6%20v01%20=%20%5B1.2%20-0.6%201%5D/StarOrbit%20Graph%20for%20M1=2%20M2=6.gif?raw=true)

#### Y-Axis Velocity
The y axis of the velocity vector now makes this interesting sawtooth shape. We can see the effects of the mass differences here quantitatively. It appears that the larger star moves at a fraction of the speed of the smaller star.

![](https://github.com/GrahamDrive/BinaryCelestialSim/blob/main/Example%20Simulations/M1=2%20M2=6%20v01%20=%20%5B1.2%20-0.6%201%5D/Y-Velocity%20Graph%20for%20M1=2%20M2=6.png?raw=true)

#### Velocity Magnitudes
The velocity magnitude graphs for these two stars are interesting as they are identical in shape but are simply scaled due to the size difference of that stars.

| Star One| Second Star|
|---------|------------|
|![](https://github.com/GrahamDrive/BinaryCelestialSim/blob/main/Example%20Simulations/M1=2%20M2=6%20v01%20=%20%5B1.2%20-0.6%201%5D/Star%201%20Velocity%20Magnitude%202%20M2=6.png?raw=true)|![](https://github.com/GrahamDrive/BinaryCelestialSim/blob/main/Example%20Simulations/M1=2%20M2=6%20v01%20=%20%5B1.2%20-0.6%201%5D/Star%202%20Velocity%20Magnitude%202%20M2=6.png?raw=true)|