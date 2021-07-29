#include <cmath>

float calculate_acute_angle_to_beamline(float xline, float yline, float zline) {
  /*
    Simple trigonometric relationship.
    Calculates the angle between: 
    - a point lying on a line which passes through the origin
    - the beam line
    This is correct only if:
    1) the origin lies at (0,0,0) and the beamline aligns with 
    the z direction.
    2) the lines passes thorugh the origin
    A more generic function is required for other cases.
   */
  float oppSide = std::sqrt( xline*xline + yline*yline );
  float adjSide = zline;
  return std::atan( oppSide / adjSide );
}
