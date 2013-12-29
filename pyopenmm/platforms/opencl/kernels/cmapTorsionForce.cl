const real PI = 3.14159265358979323846f;

// Compute the first angle.

real4 v0a = (real4) (pos1.xyz-pos2.xyz, 0.0f);
real4 v1a = (real4) (pos3.xyz-pos2.xyz, 0.0f);
real4 v2a = (real4) (pos3.xyz-pos4.xyz, 0.0f);
real4 cp0a = cross(v0a, v1a);
real4 cp1a = cross(v1a, v2a);
real cosangle = dot(normalize(cp0a), normalize(cp1a));
real angleA;
if (cosangle > 0.99f || cosangle < -0.99f) {
    // We're close to the singularity in acos(), so take the cross product and use asin() instead.

    real4 cross_prod = cross(cp0a, cp1a);
    real scale = dot(cp0a, cp0a)*dot(cp1a, cp1a);
    angleA = asin(SQRT(dot(cross_prod, cross_prod)/scale));
    if (cosangle < 0.0f)
        angleA = PI-angleA;
}
else
   angleA = acos(cosangle);
angleA = (dot(v0a, cp1a) >= 0 ? angleA : -angleA);
angleA = fmod(angleA+2.0f*PI, 2.0f*PI);

// Compute the second angle.

real4 v0b = (real4) (pos5.xyz-pos6.xyz, 0.0f);
real4 v1b = (real4) (pos7.xyz-pos6.xyz, 0.0f);
real4 v2b = (real4) (pos7.xyz-pos8.xyz, 0.0f);
real4 cp0b = cross(v0b, v1b);
real4 cp1b = cross(v1b, v2b);
cosangle = dot(normalize(cp0b), normalize(cp1b));
real angleB;
if (cosangle > 0.99f || cosangle < -0.99f) {
    // We're close to the singularity in acos(), so take the cross product and use asin() instead.

    real4 cross_prod = cross(cp0b, cp1b);
    real scale = dot(cp0b, cp0b)*dot(cp1b, cp1b);
    angleB = asin(SQRT(dot(cross_prod, cross_prod)/scale));
    if (cosangle < 0.0f)
        angleB = PI-angleB;
}
else
   angleB = acos(cosangle);
angleB = (dot(v0b, cp1b) >= 0 ? angleB : -angleB);
angleB = fmod(angleB+2.0f*PI, 2.0f*PI);

// Identify which patch this is in.

int2 pos = MAP_POS[MAPS[index]];
int size = pos.y;
real delta = 2*PI/size;
int s = (int) (angleA/delta);
int t = (int) (angleB/delta);
float4 c[4];
int coeffIndex = pos.x+4*(s+size*t);
c[0] = COEFF[coeffIndex];
c[1] = COEFF[coeffIndex+1];
c[2] = COEFF[coeffIndex+2];
c[3] = COEFF[coeffIndex+3];
real da = angleA/delta-s;
real db = angleB/delta-t;

// Evaluate the spline to determine the energy and gradients.

real torsionEnergy = 0.0f;
real dEdA = 0.0f;
real dEdB = 0.0f;
torsionEnergy = da*torsionEnergy + ((c[3].w*db + c[3].z)*db + c[3].y)*db + c[3].x;
dEdA = db*dEdA + (3.0f*c[3].w*da + 2.0f*c[2].w)*da + c[1].w;
dEdB = da*dEdB + (3.0f*c[3].w*db + 2.0f*c[3].z)*db + c[3].y;
torsionEnergy = da*torsionEnergy + ((c[2].w*db + c[2].z)*db + c[2].y)*db + c[2].x;
dEdA = db*dEdA + (3.0f*c[3].z*da + 2.0f*c[2].z)*da + c[1].z;
dEdB = da*dEdB + (3.0f*c[2].w*db + 2.0f*c[2].z)*db + c[2].y;
torsionEnergy = da*torsionEnergy + ((c[1].w*db + c[1].z)*db + c[1].y)*db + c[1].x;
dEdA = db*dEdA + (3.0f*c[3].y*da + 2.0f*c[2].y)*da + c[1].y;
dEdB = da*dEdB + (3.0f*c[1].w*db + 2.0f*c[1].z)*db + c[1].y;
torsionEnergy = da*torsionEnergy + ((c[0].w*db + c[0].z)*db + c[0].y)*db + c[0].x;
dEdA = db*dEdA + (3.0f*c[3].x*da + 2.0f*c[2].x)*da + c[1].x;
dEdB = da*dEdB + (3.0f*c[0].w*db + 2.0f*c[0].z)*db + c[0].y;
dEdA /= delta;
dEdB /= delta;
energy += torsionEnergy;

// Apply the force to the first torsion.

real normCross1 = dot(cp0a, cp0a);
real normSqrBC = dot(v1a, v1a);
real normBC = SQRT(normSqrBC);
real normCross2 = dot(cp1a, cp1a);
real dp = 1.0f/normSqrBC;
real4 ff = (real4) ((-dEdA*normBC)/normCross1, dot(v0a, v1a)*dp, dot(v2a, v1a)*dp, (dEdA*normBC)/normCross2);
real4 force1 = ff.x*cp0a;
real4 force4 = ff.w*cp1a;
real4 d = ff.y*force1 - ff.z*force4;
real4 force2 = d-force1;
real4 force3 = -d-force4;

// Apply the force to the second torsion.

normCross1 = dot(cp0b, cp0b);
normSqrBC = dot(v1b, v1b);
normBC = SQRT(normSqrBC);
normCross2 = dot(cp1b, cp1b);
dp = 1.0f/normSqrBC;
ff = (real4) ((-dEdB*normBC)/normCross1, dot(v0b, v1b)*dp, dot(v2b, v1b)*dp, (dEdB*normBC)/normCross2);
real4 force5 = ff.x*cp0b;
real4 force8 = ff.w*cp1b;
d = ff.y*force5 - ff.z*force8;
real4 force6 = d-force5;
real4 force7 = -d-force8;
