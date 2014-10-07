#include <cmath>
#include "cblas.h"
#include <glm/glm.hpp>

//Computes the distance between the two line segments define by p0, p1, q0, q1.
//The output variables gamma0 and gamma1 represent the closest points on the segments: (1-alpha)*p0 + alpha*p1 and (1-beta)*q0 + beta*q1
inline float DistanceLinesegLineseg(const glm::vec3 p0, const glm::vec3 p1, const glm::vec3 q0, const glm::vec3 q1, float& alpha, float& beta)
{
  float u[3], v[3], w[3];
  float a, b, c, d, e, D, sc, sN, sD, tc, tN, tD;
  float SMALL_NUM = 1e-10f;

  for (int k = 0; k < 3; k++)
  {
    u[k] = p1[k] - p0[k];
    v[k] = q1[k] - q0[k];
    w[k] = p0[k] - q0[k];
  }

  a = u[0]*u[0]+u[1]*u[1]+u[2]*u[2];
  b = u[0]*v[0]+u[1]*v[1]+u[2]*v[2];
  c = v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
  d = u[0]*w[0]+u[1]*w[1]+u[2]*w[2];
  e = v[0]*w[0]+v[1]*w[1]+v[2]*w[2];
  D = a*c - b*b;        // always >= 0
  sD = D;
  tD = D;

  // compute the line parameters of the two closest points
  if (D < SMALL_NUM) {  // the lines are almost parallel
      sN = 0.0;         // force using point P0 on segment S1
      sD = 1.0;         // to prevent possible division by 0.0 later
      tN = e;
      tD = c;
  }
  else {                 // get the closest points on the infinite lines
      sN = (b*e - c*d);
      tN = (a*e - b*d);
      if (sN < 0.0) {        // sc < 0 => the s=0 edge is visible
          sN = 0.0;
          tN = e;
          tD = c;
      }
      else if (sN > sD) {  // sc > 1  => the s=1 edge is visible
          sN = sD;
          tN = e + b;
          tD = c;
      }
  }

  if (tN < 0.0) {            // tc < 0 => the t=0 edge is visible
      tN = 0.0;
      // recompute sc for this edge
      if (-d < 0.0)
          sN = 0.0;
      else if (-d > a)
          sN = sD;
      else {
          sN = -d;
          sD = a;
      }
  }
  else if (tN > tD) {      // tc > 1  => the t=1 edge is visible
      tN = tD;
      // recompute sc for this edge
      if ((-d + b) < 0.0)
          sN = 0;
      else if ((-d + b) > a)
          sN = sD;
      else {
          sN = (-d +  b);
          sD = a;
      }
  }
  // finally do the division to get sc and tc
  sc = (abs(sN) < SMALL_NUM ? 0.0f : sN / sD);
  tc = (abs(tN) < SMALL_NUM ? 0.0f : tN / tD);

  // get the difference of the two closest points
  float dP[3];
  for (int k = 0; k < 3; k++)
  {
    dP[k] = w[k] + sc*u[k] - tc*v[k];
  }
  alpha = sc;
  beta  = tc;

  return cblas_snrm2(3,dP,1);   // return the closest distance
}
