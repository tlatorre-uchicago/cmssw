#include "RecoVertex/KalmanVertexFit/interface/KalmanSmoothedVertexChi2Estimator.h"
// #include "RecoVertex/KalmanVertexFit/interface/KalmanVertexTrackCompatibilityEstimator.h"


template <unsigned int N>
typename KalmanSmoothedVertexChi2Estimator<N>::BDpair
KalmanSmoothedVertexChi2Estimator<N>::estimate(const CachingVertex<N> & vertex) const
{
//initial vertex part
  float v_part = 0.;
  float returnChi = 0.;

  if (vertex.hasPrior()) {
    v_part = helper.vertexChi2(vertex.priorVertexState(), vertex.vertexState());
  }

//vector of tracks part
  typedef typename CachingVertex<N>::RefCountedVertexTrack RefCountedVertexTrack;
  std::vector< RefCountedVertexTrack > tracks = vertex.tracks();
  float sum = 0.;
  bool success = true;
  for(typename std::vector<RefCountedVertexTrack>::iterator i = tracks.begin(); i != tracks.end(); i++)
  {
    BDpair result = helper.trackParameterChi2((*i)->linearizedTrack(), (*i)->refittedState());
    success = success && result.first;
    sum += (*i)->weight() * result.second;
  }
 returnChi = v_part + sum;
 if (returnChi < 0) {
    fprintf(stderr, "v_part = %f\n", v_part);
    fprintf(stderr, "sum = %f\n", sum);
    fprintf(stderr, "KalmanSmoothedVertexChi2Estimator: returnChi < 0!\n");
      for(typename std::vector<RefCountedVertexTrack>::iterator i = tracks.begin(); i != tracks.end(); i++)
      {
        BDpair result = helper.trackParameterChi2((*i)->linearizedTrack(), (*i)->refittedState());
        success = success && result.first;
        fprintf(stderr, "(*i)->weight = %f\n", (*i)->weight());
        fprintf(stderr, "second = %f\n", result.second);
      }
 }
 return BDpair(success, returnChi);
}

template class KalmanSmoothedVertexChi2Estimator<5>;
template class KalmanSmoothedVertexChi2Estimator<6>;
