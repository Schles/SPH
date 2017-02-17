#include "BoundaryPsi.h"

BoundaryPsi::BoundaryPsi(CompactNSearch::NeighborhoodSearch* nhSearch, ParticleManager* manager, Parameters *params, Kernel *kernel) :
  m_neighborhoodsearch(nhSearch),
  _parameters(params),
  particleManager(manager),
  m_kernel(kernel)
{
  
}

void BoundaryPsi::calc_boundary_psi(){

  // disable all 
  for(int pointSet = 0; pointSet < particleManager->getParticleGroupSize(); pointSet++){
      m_neighborhoodsearch->point_set(pointSet).enable_neighborsearch(false);
  }

  // search only static
  for(int pointSet = 0; pointSet <  particleManager->getParticleGroupSize(); pointSet++){
    if(! ( particleManager->getObject(pointSet)->isDynamic || particleManager->getObject(pointSet)->isFluid))
      m_neighborhoodsearch->point_set(pointSet).enable_neighborsearch(true);
  }

  m_neighborhoodsearch->find_neighbors();

  for(int pointSet = 0; pointSet <  particleManager->getParticleGroupSize(); pointSet++){
    if(! ( particleManager->getObject(pointSet)->isDynamic || particleManager->getObject(pointSet)->isFluid))
      calc_boundary_psi_i(pointSet);
  }

  
  
  for(int pointSet = 0; pointSet < particleManager->getParticleGroupSize(); pointSet++){
    for(int j = 0; j < particleManager->getParticleGroupSize(); j++){
      m_neighborhoodsearch->point_set(j).enable_neighborsearch(false);
    }

    if(particleManager->getObject(pointSet)->isDynamic && !particleManager->getObject(pointSet)->isFluid){
      m_neighborhoodsearch->point_set(pointSet).enable_neighborsearch(true);
      m_neighborhoodsearch->find_neighbors();
      calc_boundary_psi_i(pointSet);
    }

  }

  // Activate Fluids
  for(int pointSet = 0; pointSet < particleManager->getParticleGroupSize(); pointSet++){
    if(particleManager->getObject(pointSet)->isFluid)
      m_neighborhoodsearch->point_set(pointSet).enable_neighborsearch(true);
    else
      m_neighborhoodsearch->point_set(pointSet).enable_neighborsearch(false);

  }


}


void BoundaryPsi::calc_boundary_psi_i(unsigned int pointSet){
    CompactNSearch::PointSet const& ps = m_neighborhoodsearch->point_set(pointSet);

    Boundary* boundaryObject = particleManager->getBoundaryObject(pointSet);
    
    for(unsigned int i = 0; i < particleManager->getObjectSize(pointSet); i++){

      double delta = 0;
      
      for (int k = 0; k < ps.n_neighbors(i); k++) {

	CompactNSearch::PointID nid = ps.neighbor(i, k);
	
	Eigen::Vector3d q = particleManager->getPosition(pointSet, i) - particleManager->getPosition(nid.point_set_id, nid.point_id);


	if(!particleManager->isFluid(nid.point_set_id)){
	  delta += m_kernel->w(q, _parameters->h);
	}
      }
      
      double volume = 0.0;

      if(delta > 0)
	volume = 1.0 / delta;
      
      boundaryObject->m_boundaryPsi[i] = _parameters->restDensity * volume;
    }
    
}
