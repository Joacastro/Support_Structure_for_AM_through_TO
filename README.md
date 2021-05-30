# Generating Supports Structures for Additive Manufacturing through Topology optimization 
## This is a MATLAB code to generate support using topology optimization with length scale control

This code use the density approach to generate efficiente and easy to remove supports structures for additive manufacturing.
The problem is formulated as is shown in the follow figure and considering the a vertical build direction.

![image](https://user-images.githubusercontent.com/52135295/120120254-02a51700-c16a-11eb-8d17-a5aa17a85a28.png)

We focus in generate self-supporting support structures. For that, we gather the geometric guidelines of the Ultimaker S2+ 3D Printer 
by printing the different pieces to get the following specifications.

- Minimum wall thickness
- Minimum hole diameter 
- Minimum separation between walls
- Optimum aspect ratio
- LPBL (Longest Printable Bridge Length)

We included this specifications as a restrictions in our topology optimization framework.

## Base Code
The Matlab codes are based in the [88 lines Matlab code for TO](https://www.topopt.mek.dtu.dk/Apps-and-software/Efficient-topology-optimization-in-MATLAB). The codes are adapted to inclued the additive manufacturing contraints developed by Fernandez et. al. in [An aggregation strategy of maximum size constraints in density-based topology optimization](https://link.springer.com/article/10.1007/s00158-019-02313-8) and [Imposing minimum and maximum member size, minimum cavity size, and minimum separation distance between solid members in topology optimization](https://www.sciencedirect.com/science/article/abs/pii/S004578252030342X).
