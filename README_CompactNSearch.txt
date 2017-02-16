

CompactNSearch.h/cpp contains a class to efficiently compute neighborhood information of particles in a given radius r around each particle.

The class is initialized using the constructor via

CompactNSearch(Real r, bool erase_empty_cells = false);,

where "r" represents the search radius. The second parameter's value can remain "false".

The implementation supports the usage of multiple particle systems in order to be able to distinguish between fluid and boundary systems (or even different fluid phases). Usually, no neighborhood information for boundary particles has to be computed which reduces the computational effort.

Using the method

unsigned int add_discretization(Vector3d* x_, std::size_t n_, bool is_moving_, bool search_neighbors_)

a particle system can be added to the object. The particle coordinates are passed to the function using a pointer "x" to a contiguous array of three-dimensional vectors. The number of particles is defined by the second parameter "n". "is_moving" defines if the particles of the added set are dynamic or static. The last parameter should only be set to "false" for boundary particle systems. Finally, the method is returning an id that should be stored by the user in order to associate the passed particle set with the returned id.

The actual neighborhood search can then be invoked using the method

void neighborhood_search();.

Note that the invocation for fluid simulations is usually performed only once per timestep.

In order to finally determine the neighboring particles of particle "i" of a particle system with id "j" the discretizations array has to be retrieved via

std::vector<SPHDiscretization> const& discretizations() const;.

Using the "j"th object in the vector, the number of particles neighboring particle "i" can then be determined using

unsigned int n_neighbors(unsigned int i) const.

By calling

ParticleID const& neighbor(unsigned int i, unsigned int k) const
{
    return m_neighbors[i][k];
},

the index of a neighboring particle can be retrieved, where k € [0,  n_neighbors(i)). The method returns a ParticleID object which stores a pair of integers representing the particle system id and the particle index of the neighboring particle.




The usage of the following information is optional and therefore not required in order to correctly determine particle neighborhood information!

The class is additionally able to sort the particles according to a space-filling z-curve. Therefore, spatially neighboring particles are also spatially close in memory which drastically reduces cache misses. In order to apply the sort algorithm the method

void z_sort();

has to be called. During the simulation the particles move and should therefore be repeatedly sorted after a user-defined number of timesteps (~1000). An even more frequent usage of the sort algorithm is not advantageous as the algorithm is computationally rather expensive.
As the particles carry additional information, e.g. mass, density, velocity etc., the according arrays also have to be sorted in order to remain consistent with the coordinate array. Therefore, the user has to call the method

template <typename T>
void SPHDiscretization::sort_field(T* lst) const;

after invocation of "z_sort()" for each property array (including the coordinate array). The parameter passed to the method is a pointer to the first element of the contiguous property array.
