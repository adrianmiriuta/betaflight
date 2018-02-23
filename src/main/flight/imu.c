/*
 * This file is part of Cleanflight.
 *
 * Cleanflight is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Cleanflight is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Cleanflight.  If not, see <http://www.gnu.org/licenses/>.
 */

// Inertial Measurement Unit (IMU)

#include <stdbool.h>
#include <stdint.h>
#include <math.h>

#include "platform.h"

#include "build/build_config.h"
#include "build/debug.h"

#include "common/axis.h"

#include "pg/pg.h"
#include "pg/pg_ids.h"

#include "drivers/time.h"

#include "fc/runtime_config.h"

#include "flight/imu.h"
#include "flight/mixer.h"
#include "flight/pid.h"

#include "io/gps.h"

#include "sensors/acceleration.h"
#include "sensors/barometer.h"
#include "sensors/compass.h"
#include "sensors/gyro.h"
#include "sensors/sensors.h"

#if defined(SIMULATOR_BUILD) && defined(SIMULATOR_MULTITHREAD)
#include <stdio.h>
#include <pthread.h>

static pthread_mutex_t imuUpdateLock;

#if defined(SIMULATOR_IMU_SYNC)
static uint32_t imuDeltaT = 0;
static bool imuUpdated = false;
#endif

#define IMU_LOCK pthread_mutex_lock(&imuUpdateLock)
#define IMU_UNLOCK pthread_mutex_unlock(&imuUpdateLock)

#else

#define IMU_LOCK
#define IMU_UNLOCK

#endif

// the limit (in degrees/second) beyond which we stop integrating
// omega_I. At larger spin rates the DCM PI controller can get 'dizzy'
// which results in false gyro drift. See
// http://gentlenav.googlecode.com/files/fastRotations.pdf

#define SPIN_RATE_LIMIT 20

int32_t accSum[XYZ_AXIS_COUNT];

uint32_t accTimeSum = 0;        // keep track for integration of acc
int accSumCount = 0;
float accVelScale;

static float throttleAngleScale;
static float fc_acc;
static float smallAngleCosZ = 0;

static imuRuntimeConfig_t imuRuntimeConfig;


// quaternion of sensor frame relative to earth frame
STATIC_UNIT_TESTED quaternion qAttitude = QUATERNION_INITIALIZE;
STATIC_UNIT_TESTED quaternion qGyro = QUATERNION_INITIALIZE;
STATIC_UNIT_TESTED quaternion qGyroB = QUATERNION_INITIALIZE;
STATIC_UNIT_TESTED quaternion qGyroBinverse = QUATERNION_INITIALIZE;
STATIC_UNIT_TESTED quaternionProducts qpAttitude = QUATERNION_PRODUCTS_INITIALIZE;
STATIC_UNIT_TESTED quaternionProducts qpGyro = QUATERNION_PRODUCTS_INITIALIZE;
// headfree quaternions
quaternion qHeadfree = QUATERNION_INITIALIZE;
quaternion qOffset = QUATERNION_INITIALIZE;

// absolute angle inclination in multiple of 0.1 degree    180 deg = 1800
attitudeEulerAngles_t attitude = EULER_INITIALIZE;

PG_REGISTER_WITH_RESET_TEMPLATE(imuConfig_t, imuConfig, PG_IMU_CONFIG, 0);

PG_RESET_TEMPLATE(imuConfig_t, imuConfig,
    .dcm_kp = 2500,                // 1.0 * 10000
    .dcm_ki = 0,                   // 0.003 * 10000
    .small_angle = 25,
    .accDeadband = {.xy = 40, .z= 40},
    .acc_unarmedcal = 1
);

/*
* Calculate RC time constant used in the accZ lpf.
*/
static float calculateAccZLowPassFilterRCTimeConstant(float accz_lpf_cutoff)
{
    return 0.5f / (M_PIf * accz_lpf_cutoff);
}

static float calculateThrottleAngleScale(uint16_t throttle_correction_angle)
{
    return (1800.0f / M_PIf) * (900.0f / throttle_correction_angle);
}

void imuConfigure(uint16_t throttle_correction_angle)
{
    imuRuntimeConfig.dcm_kp = imuConfig()->dcm_kp / 10000.0f;
    imuRuntimeConfig.dcm_ki = imuConfig()->dcm_ki / 10000.0f;
    imuRuntimeConfig.acc_unarmedcal = imuConfig()->acc_unarmedcal;
    imuRuntimeConfig.small_angle = imuConfig()->small_angle;

    fc_acc = calculateAccZLowPassFilterRCTimeConstant(5.0f); // Set to fix value
    throttleAngleScale = calculateThrottleAngleScale(throttle_correction_angle);
}

void imuInit(void)
{
    smallAngleCosZ = cos_approx(degreesToRadians(imuRuntimeConfig.small_angle));
    accVelScale = 9.80665f / acc.dev.acc_1G / 10000.0f;


#if defined(SIMULATOR_BUILD) && defined(SIMULATOR_MULTITHREAD)
    if (pthread_mutex_init(&imuUpdateLock, NULL) != 0) {
        printf("Create imuUpdateLock error!\n");
    }
#endif
}

void imuResetAccelerationSum(void)
{
    accSum[0] = 0;
    accSum[1] = 0;
    accSum[2] = 0;
    accSumCount = 0;
    accTimeSum = 0;
}

#if defined(USE_ALT_HOLD)
// rotate acc into Earth frame and calculate acceleration in it
static void imuCalculateAcceleration(timeDelta_t deltaT)
{
    static float accZoffset = 0;
    static float accz_smooth = 0;

    // deltaT is measured in us ticks
    const float dT = (float)deltaT * 1e-6f;

    quaternion accel_ned;
    accel_ned.x = acc.accADC[X];
    accel_ned.y = acc.accADC[Y];
    accel_ned.z = acc.accADC[Z];
    quaternionTransformVectorBodyToEarth(&accel_ned, &qAttitude);

    if (imuRuntimeConfig.acc_unarmedcal == 1) {
        if (!ARMING_FLAG(ARMED)) {
            accZoffset -= accZoffset / 64;
            accZoffset += accel_ned.z;
        }
        accel_ned.z -= accZoffset / 64;  // compensate for gravitation on z-axis
    } else {
        accel_ned.z -= acc.dev.acc_1G;
    }

    accz_smooth = accz_smooth + (dT / (fc_acc + dT)) * (accel_ned.z - accz_smooth); // low pass filter

    // apply Deadband to reduce integration drift and vibration influence
    accSum[X] += applyDeadband(lrintf(accel_ned.x), imuRuntimeConfig.accDeadband.xy);
    accSum[Y] += applyDeadband(lrintf(accel_ned.y), imuRuntimeConfig.accDeadband.xy);
    accSum[Z] += applyDeadband(lrintf(accz_smooth), imuRuntimeConfig.accDeadband.z);

    // sum up Values for later integration to get velocity and distance
    accTimeSum += deltaT;
    accSumCount++;
}
#endif // USE_ALT_HOLD

static float invSqrt(float x)
{
    return 1.0f / sqrtf(x);
}

static bool imuUseFastGains(void)
{
    return !ARMING_FLAG(ARMED) && millis() < 20000;
}

static float imuGetPGainScaleFactor(void)
{
    if (imuUseFastGains()) {
        return 10.0f;
    }
    else {
        return 1.0f;
    }
}

static void imuMahonyAHRSupdate(float dt, float gx, float gy, float gz,
                                bool useAcc, float ax, float ay, float az,
                                bool useMag, float mx, float my, float mz,
                                bool useYaw, float yawError)
{
    /*
    static float integralFBx = 0.0f,  integralFBy = 0.0f, integralFBz = 0.0f;    // integral error terms scaled by Ki

    // Calculate general spin rate (rad/s)
    const float spin_rate = sqrtf(sq(gx) + sq(gy) + sq(gz));

    // Use raw heading error (from GPS or whatever else)
    float ex = 0, ey = 0, ez = 0;
    if (useYaw) {
        while (yawError >  M_PIf) yawError -= (2.0f * M_PIf);
        while (yawError < -M_PIf) yawError += (2.0f * M_PIf);
        ez += sin_approx(yawError / 2.0f);
    }

#ifdef USE_MAG
    // Use measured magnetic field vector
    float recipMagNorm = sq(mx) + sq(my) + sq(mz);
    if (useMag && recipMagNorm > 0.01f) {
        // Normalise magnetometer measurement
        recipMagNorm = invSqrt(recipMagNorm);
        mx *= recipMagNorm;
        my *= recipMagNorm;
        mz *= recipMagNorm;

        // For magnetometer correction we make an assumption that magnetic field is perpendicular to gravity (ignore Z-component in EF).
        // This way magnetic field will only affect heading and wont mess roll/pitch angles

        // (hx; hy; 0) - measured mag field vector in EF (assuming Z-component is zero)
        // (bx; 0; 0) - reference mag field vector heading due North in EF (assuming Z-component is zero)
        const float hx = (1.0f - 2.0f * qpAttitude.yy - 2.0f * qpAttitude.zz) * mx + (2.0f * (qpAttitude.xy + -qpAttitude.wz)) * my + (2.0f * (qpAttitude.xz - -qpAttitude.wy)) * mz;
        const float hy = (2.0f * (qpAttitude.xy - -qpAttitude.wz)) * mx + (1.0f - 2.0f * qpAttitude.xx - 2.0f * qpAttitude.zz) * my + (2.0f * (qpAttitude.yz + -qpAttitude.wx)) * mz;
        const float bx = sqrtf(hx * hx + hy * hy);

        // magnetometer error is cross product between estimated magnetic north and measured magnetic north (calculated in EF)
        const float ez_ef = -(hy * bx);

        // Rotate mag error vector back to BF and accumulate
        ex += (2.0f * (qpAttitude.xz + -qpAttitude.wy)) * ez_ef;
        ey += (2.0f * (qpAttitude.yz - -qpAttitude.wx)) * ez_ef;
        ez += (1.0f - 2.0f * qpAttitude.xx - 2.0f * qpAttitude.yy) * ez_ef;
    }
#else
    UNUSED(useMag);
    UNUSED(mx);
    UNUSED(my);
    UNUSED(mz);
#endif

    // Use measured acceleration vector
    float recipAccNorm = sq(ax) + sq(ay) + sq(az);
    if (useAcc && recipAccNorm > 0.01f) {
        // Normalise accelerometer measurement
        recipAccNorm = invSqrt(recipAccNorm);
        ax *= recipAccNorm;
        ay *= recipAccNorm;
        az *= recipAccNorm;

        // Error is sum of cross product between estimated direction and measured direction of gravity
        ex += (ay * (1.0f - 2.0f * qpAttitude.xx - 2.0f * qpAttitude.yy) - az * (2.0f * (qpAttitude.yz - -qpAttitude.wx)));
        ey += (az * (2.0f * (qpAttitude.xz + -qpAttitude.wy)) - ax * (1.0f - 2.0f * qpAttitude.xx - 2.0f * qpAttitude.yy));
        ez += (ax * (2.0f * (qpAttitude.yz - -qpAttitude.wx)) - ay * (2.0f * (qpAttitude.xz + -qpAttitude.wy)));
    }

    // Compute and apply integral feedback if enabled
    if (imuRuntimeConfig.dcm_ki > 0.0f) {
        // Stop integrating if spinning beyond the certain limit
        if (spin_rate < DEGREES_TO_RADIANS(SPIN_RATE_LIMIT)) {
            const float dcmKiGain = imuRuntimeConfig.dcm_ki;
            integralFBx += dcmKiGain * ex * dt;    // integral error scaled by Ki
            integralFBy += dcmKiGain * ey * dt;
            integralFBz += dcmKiGain * ez * dt;
        }
    } else {
        integralFBx = 0.0f;    // prevent integral windup
        integralFBy = 0.0f;
        integralFBz = 0.0f;
    }

    // Calculate kP gain. If we are acquiring initial attitude (not armed and within 20 sec from powerup) scale the kP to converge faster
    const float dcmKpGain = imuRuntimeConfig.dcm_kp * imuGetPGainScaleFactor();

    // Apply proportional and integral feedback
    gx += dcmKpGain * ex + integralFBx;
    gy += dcmKpGain * ey + integralFBy;
    gz += dcmKpGain * ez + integralFBz;
    */


    // old bf method
    // has positions of high drift +-90° 45-45°
    // identical with adapted version (no diff with rot inverted)
    // Integrate rate of change of quaternion
    /*
    gx *= (0.5f * dt);
    gy *= (0.5f * dt);
    gz *= (0.5f * dt);
    quaternion buffer;
    quaternionCopy(&qAttitude, &buffer);
    // construct new quaternion from old quaternion and rate of change gyro data
    qAttitude.w += (-buffer.x * gx - buffer.y * gy - buffer.z * gz);
    qAttitude.x += (+buffer.w * gx + buffer.y * gz - buffer.z * gy);
    qAttitude.y += (+buffer.w * gy - buffer.x * gz + buffer.z * gx);
    qAttitude.z += (+buffer.w * gz + buffer.x * gy - buffer.y * gx);
    quaternionNormalize(&qAttitude);
    */

    // old bf method adapted
    // has positions of high drift +-90° 45-45°
    // problem high drift around +-90° drift pitch roll 1°/s
    // old BF vs old BF adapted no diff
    quaternion qBuff, qDiff;
    qDiff.w = 0;
    qDiff.x = gx * 0.5f * dt;
    qDiff.y = gy * 0.5f * dt;
    qDiff.z = gz * 0.5f * dt;
    quaternionMultiply(&qGyro, &qDiff, &qBuff);
    quaternionAdd(&qGyro, &qBuff, &qGyro);
    quaternionNormalize(&qGyro);
    //quaternionInverse(&qGyroB, &qGyroBinverse);


    // test method b
    // some sort of aproximation!!!
    // https://github.com/malloch/Arduino_IMU
    // problem singularities around +-90° without normalization
    // problem singularities around +-90° not sin_approx related
    // probelm high drift 45°45°
    // subjective lower static drift 0.1°/s
    // variation qDiffNorm Ko
    // variation sin cos Ko
    // no drift 0.1° differeces roll pitch +-90° to BF calculus
    /*
    //quaternion qDiff;
    //const float qDiffNorm = sqrt(gx*gx + gy*gy + gz*gz);
    qDiff.w = cos_approx((gx + gy + gz) * 0.5f * dt);
    //qDiff.w = cos(qDiffNorm * 0.5f * dt);
    qDiff.x = sin(gx * 0.5f * dt);
    qDiff.y = sin(gy * 0.5f * dt);
    qDiff.z = sin(gz * 0.5f * dt);
    //quaternionMultiply(&qGyro, &qDiff, &qGyro);
    //quaternionNormalize(&qGyro);
    quaternionMultiply(&qGyroB, &qDiff, &qGyroB);
    quaternionInverse(&qGyroB, &qGyroBinverse);
    */


    // test method c
    // https://math.stackexchange.com/questions/1693067/differences-between-quaternion-integration-methods
    // more diff than malloch vs BF calculus (on high speed movements) no drift
    /*
    //quaternion qDiff;
    const float qDiffNorm = sqrt(gx*gx + gy*gy + gz*gz);
    if (qDiffNorm > 0.0000001f) {
      qDiff.w = cos(qDiffNorm * 0.5f * dt);
      qDiff.x = (gx * sin(qDiffNorm * 0.5f * dt)) / qDiffNorm;
      qDiff.y = (gy * sin(qDiffNorm * 0.5f * dt)) / qDiffNorm;
      qDiff.z = (gz * sin(qDiffNorm * 0.5f * dt)) / qDiffNorm;
      //quaternionMultiply(&qGyro, &qDiff, &qGyro);
      quaternionMultiply(&qGyroB, &qDiff, &qGyroB);
      quaternionInverse(&qGyroB, &qGyroBinverse);
    }
    */


    // test method d
    // test incremental rotation
    // singularities circle around +-90° sin_approx cos_approx related
    // worse than BF when shaiking gently (more diff from should be position)
    // incremental vs BF large diff vs BF calculus (on higher speed movements)
    /*
    quaternion qDiff;
    qDiff.w = cos(gx * dt * 0.5f);
    qDiff.x = sin(gx * dt * 0.5f);
    qDiff.y = 0;
    qDiff.z = 0;
    quaternionMultiply(&qGyro, &qDiff, &qGyro);
    //quaternionMultiply(&qGyroB, &qDiff, &qGyroB);
    qDiff.w = cos(gy * dt * 0.5f);
    qDiff.x = 0;
    qDiff.y = sin(gy * dt * 0.5f);
    qDiff.z = 0;
    quaternionMultiply(&qGyro, &qDiff, &qGyro);
    //quaternionMultiply(&qGyroB, &qDiff, &qGyroB);
    qDiff.w = cos(gz * dt * 0.5f);
    qDiff.x = 0;
    qDiff.y = 0;
    qDiff.z = sin(gz * dt * 0.5f);
    quaternionMultiply(&qGyro, &qDiff, &qGyro);
    //quaternionMultiply(&qGyroB, &qDiff, &qGyroB);
    //quaternionInverse(&qGyroB, &qGyroBinverse);
    */


    // test method e
    // single rotation quaternion
    // singularity +-90° when using sin_approx
    // single rot vs bf same as incremental vs BF
    // incremental rot vs single rot big difference (rotation order???)
    /*
    //quaternion qDiff;
    const float cy = cos(gz * dt * 0.5);
    const float sy = sin(gz * dt * 0.5);
    const float cr = cos(gx * dt * 0.5);
    const float sr = sin(gx * dt * 0.5);
    const float cp = cos(gy * dt * 0.5);
    const float sp = sin(gy * dt * 0.5);

    qDiff.w = cy * cr * cp + sy * sr * sp;
    qDiff.x = cy * sr * cp - sy * cr * sp;
    qDiff.y = cy * cr * sp + sy * sr * cp;
    qDiff.z = sy * cr * cp - cy * sr * sp;
    //quaternionMultiply(&qGyro, &qDiff, &qGyro);

    quaternionMultiply(&qGyroB, &qDiff, &qGyroB);
    quaternionInverse(&qGyroB, &qGyroBinverse);
    */


    // acc
    quaternion vAcc, qAcc;
    quaternionProducts qpAcc;
    vAcc.w = 0;
    vAcc.x = ax;
    vAcc.y = ay;
    vAcc.z = az;
    quaternionNormalize(&vAcc);

    /*
    quaternion qAccRoll;
    float u = 0.1f;
    //float rollHalf = atan2_approx(vAcc.y,vAcc.z)/2; // mmax v1 ROT xyz
    //float rollHalf = atan(vAcc.y/vAcc.z)/2; // mmax v1 ROT xyz
    //float rollHalf = atan2_approx(vAcc.y,sqrtf(vAcc.z*vAcc.z + vAcc.x*vAcc.x))/2; // mmax v2 ROT xyz only z>0 no inverted position
    float rollHalf = atan2(vAcc.y,(float)copysign(1.0f,vAcc.z) * sqrtf(vAcc.z*vAcc.z + u*vAcc.x*vAcc.x))/2.0f; // AN3461 xyz
    //float rollHalf = atan(vAcc.y/copysign(1.0,vAcc.z)*sqrtf(vAcc.z*vAcc.z + u*vAcc.x*vAcc.x))/2; // AN3461 xyz

    qAccRoll.w = cos(rollHalf);
    qAccRoll.x = sin(rollHalf);
    qAccRoll.y = 0;
    qAccRoll.z = 0;

    quaternion qAccPitch;
    float pitchHalf= atan2(-vAcc.x, sqrtf(vAcc.y * vAcc.y + vAcc.z * vAcc.z) )/2.0f; // AN3461 xyz
    //float pitchHalf= atan(-vAcc.x/sqrtf(vAcc.y * vAcc.y + vAcc.z * vAcc.z) )/2; // AN3461 xyz
    //float pitchHalf= atan(-vAcc.x/vAcc.z)/2; // AN3461 xyz
    //pitchHalf= constrainf(pitch2, -M_PIf4, M_PIf4);
    qAccPitch.w = cos(pitchHalf);
    qAccPitch.x = 0;
    qAccPitch.y = sin(pitchHalf);
    qAccPitch.z = 0;

    quaternionMultiply(&qAccRoll, &qAccPitch, &qAcc); //xyz
    */


    //https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4570372/

    if (imuIsAccelerometerHealthy()) {
      //if (vAcc.z >= -0.9) {
      qAcc.w = +sqrtf((vAcc.z + 1) / 2.0f);
      qAcc.x = +vAcc.y/(2 * qAcc.w);
      qAcc.y = -vAcc.x/(2 * qAcc.w);
      qAcc.z = 0;
      /*} else {
      // Ko
      qAcc.w = +vAcc.y/sqrtf(2.0f * (1 - vAcc.z));
      qAcc.x = +sqrtf((1 - vAcc.z)/2.0f);
      qAcc.y = 0;
      qAcc.z = +vAcc.x/sqrtf(2.0f * (1 - vAcc.z));
      }*/

    }



    quaternionNormalize(&qAcc);

    quaternionComputeProducts(&qAcc, &qpAcc);
    quaternion qAccYaw;
    const float AccYawHalf = atan2((+2.0f * (qpAcc.wz + qpAcc.xy)), (+1.0f - 2.0f * (qpAcc.yy + qpAcc.zz))) / 2.0f;
    qAccYaw.w = cos(AccYawHalf);
    qAccYaw.x = 0;
    qAccYaw.y = 0;
    qAccYaw.z = sin(AccYawHalf);
    quaternionInverse(&qAccYaw,&qAccYaw);
    // remove yaw rotation
    //quaternionMultiply(&qAccYaw, &qAcc, &qAcc);

    //introduces roll pitch drift !
    /*
    quaternionInverse(&qAcc, &qAccInverse);
    quaternion qGyroYaw;
    quaternionMultiply(&qGyro, &qAccInverse, &qGyroYaw);
    quaternionMultiply(&qAcc, &qGyroYaw, &qAcc);*/


    quaternionComputeProducts(&qGyro, &qpGyro);
    quaternion qGyroYaw;
    const float yaw = atan2_approx((+2.0f * (qpGyro.wz + qpGyro.xy)), (+1.0f - 2.0f * (qpGyro.yy + qpGyro.zz)));
    qGyroYaw.w = cos_approx(yaw/2);
    qGyroYaw.x = 0;
    qGyroYaw.y = 0;
    qGyroYaw.z = sin_approx(yaw/2);
    //quaternionMultiply(&qGyroYaw, &qAcc, &qAcc);


    //quaternionMinimumDistance(&qAcc, &qGyro);
    //quaternionSlerp(&qAcc, &qGyro,  &qAttitude, 0.995);
    //ko
    //quaternionSlerp(&qAcc, &qGyro,  &qAttitude, constrainf(quaternionDotProduct(&qAcc, &qAttitude),0.5f,0.999f));
    //quaternionCopy(&qAttitude, &qGyro);


    quaternionCopy(&qAcc, &qAttitude);
    //quaternionNormalize(&qAttitude);
    //quaternionCopy(&qGyro, &qAttitude);
    //quaternionMultiply(&qGyroBinverse, &qGyro, &qAttitude);
    quaternionComputeProducts(&qAttitude, &qpAttitude);
    quaternionCopy(&qAttitude, &qGyro);
}

STATIC_UNIT_TESTED void imuUpdateEulerAngles(void) {
    quaternionProducts buffer;

    if (FLIGHT_MODE(HEADFREE_MODE)) {
        quaternionMultiply(&qOffset, &qAttitude, &qHeadfree);
        quaternionComputeProducts(&qHeadfree, &buffer);
    } else {
        quaternionComputeProducts(&qAttitude, &buffer);
    }

    attitude.values.roll = lrintf(atan2_approx((+2.0f * (buffer.wx + buffer.yz)), (+1.0f - 2.0f * (buffer.xx + buffer.yy))) * (1800.0f / M_PIf));
    attitude.values.pitch = lrintf(((0.5f * M_PIf) - acos_approx(+2.0f * (buffer.wy - buffer.xz))) * (1800.0f / M_PIf));
    attitude.values.yaw = lrintf((-atan2_approx((+2.0f * (buffer.wz + buffer.xy)), (+1.0f - 2.0f * (buffer.yy + buffer.zz))) * (1800.0f / M_PIf)));

    if (attitude.values.yaw < 0)
        attitude.values.yaw += 3600;

    if (getCosTiltAngle() > smallAngleCosZ) {
        ENABLE_STATE(SMALL_ANGLE);
    } else {
        DISABLE_STATE(SMALL_ANGLE);
    }
}

static bool imuIsAccelerometerHealthy(void)
{
    float accMagnitude = 0;
    for (int axis = 0; axis < 3; axis++) {
        const float a = acc.accADC[axis];
        accMagnitude += a * a;
    }

    accMagnitude = accMagnitude * 100 / (sq((int32_t)acc.dev.acc_1G));

    // Accept accel readings only in range 0.90g - 1.10g
    return (81 < accMagnitude) && (accMagnitude < 121);
}

static void imuCalculateEstimatedAttitude(timeUs_t currentTimeUs)
{
    static timeUs_t previousIMUUpdateTime;
    float rawYawError = 0;
    bool useAcc = false;
    bool useMag = false;
    bool useYaw = false;

    const timeDelta_t deltaT = currentTimeUs - previousIMUUpdateTime;
    previousIMUUpdateTime = currentTimeUs;

    if (imuIsAccelerometerHealthy()) {
        useAcc = true;
    }

#ifdef USE_MAG
    if (sensors(SENSOR_MAG) && compassIsHealthy()) {
        useMag = true;
    }
#endif
#if defined(USE_GPS)
    else if (STATE(FIXED_WING) && sensors(SENSOR_GPS) && STATE(GPS_FIX) && gpsSol.numSat >= 5 && gpsSol.groundSpeed >= 300) {
        // In case of a fixed-wing aircraft we can use GPS course over ground to correct heading
        rawYawError = DECIDEGREES_TO_RADIANS(attitude.values.yaw - gpsSol.groundCourse);
        useYaw = true;
    }
#endif

#if defined(SIMULATOR_BUILD) && defined(SKIP_IMU_CALC)
    UNUSED(imuMahonyAHRSupdate);
    UNUSED(useAcc);
    UNUSED(useMag);
    UNUSED(useYaw);
    UNUSED(rawYawError);
#else

#if defined(SIMULATOR_BUILD) && defined(SIMULATOR_IMU_SYNC)
//  printf("[imu]deltaT = %u, imuDeltaT = %u, currentTimeUs = %u, micros64_real = %lu\n", deltaT, imuDeltaT, currentTimeUs, micros64_real());
    deltaT = imuDeltaT;
#endif
    float gyroAverage[XYZ_AXIS_COUNT];
    gyroGetAccumulationAverage(gyroAverage);
    float accAverage[XYZ_AXIS_COUNT];
    if (!accGetAccumulationAverage(accAverage)) {
        useAcc = false;
    }
    imuMahonyAHRSupdate(deltaT * 1e-6f,
                        DEGREES_TO_RADIANS(gyroAverage[X]), DEGREES_TO_RADIANS(gyroAverage[Y]), DEGREES_TO_RADIANS(gyroAverage[Z]),
                        useAcc, accAverage[X], accAverage[Y], accAverage[Z],
                        useMag, mag.magADC[X], mag.magADC[Y], mag.magADC[Z],
                        useYaw, rawYawError);

    imuUpdateEulerAngles();
#endif
#if defined(USE_ALT_HOLD)
    imuCalculateAcceleration(deltaT); // rotate acc vector into earth frame
#endif
}

void imuUpdateAttitude(timeUs_t currentTimeUs)
{
    if (sensors(SENSOR_ACC) && acc.isAccelUpdatedAtLeastOnce) {
        IMU_LOCK;
#if defined(SIMULATOR_BUILD) && defined(SIMULATOR_IMU_SYNC)
        if (imuUpdated == false) {
            IMU_UNLOCK;
            return;
        }
        imuUpdated = false;
#endif
        imuCalculateEstimatedAttitude(currentTimeUs);
        IMU_UNLOCK;
    } else {
        acc.accADC[X] = 0;
        acc.accADC[Y] = 0;
        acc.accADC[Z] = 0;
    }
}

float getCosTiltAngle(void) {
    return (1.0f - 2.0f * (qpAttitude.xx + qpAttitude.yy));
}

int16_t calculateThrottleAngleCorrection(uint8_t throttle_correction_value)
{
    /*
    * Use 0 as the throttle angle correction if we are inverted, vertical or with a
    * small angle < 0.86 deg
    * TODO: Define this small angle in config.
    */
    if (getCosTiltAngle() <= 0.015f) {
        return 0;
    }
    int angle = lrintf(acos_approx(getCosTiltAngle()) * throttleAngleScale);
    if (angle > 900)
        angle = 900;
    return lrintf(throttle_correction_value * sin_approx(angle / (900.0f * M_PIf / 2.0f)));
}

#ifdef SIMULATOR_BUILD
void imuSetAttitudeRPY(float roll, float pitch, float yaw)
{
    IMU_LOCK;

    attitude.values.roll = roll * 10;
    attitude.values.pitch = pitch * 10;
    attitude.values.yaw = yaw * 10;

    IMU_UNLOCK;
}
void imuSetAttitudeQuat(float w, float x, float y, float z)
{
    IMU_LOCK;

    qAttitude.w = w;
    qAttitude.x = x;
    qAttitude.y = y;
    qAttitude.z = z;

    imuUpdateEulerAngles();

    IMU_UNLOCK;
}
#endif
#if defined(SIMULATOR_BUILD) && defined(SIMULATOR_IMU_SYNC)
void imuSetHasNewData(uint32_t dt)
{
    IMU_LOCK;

    imuUpdated = true;
    imuDeltaT = dt;

    IMU_UNLOCK;
}
#endif

bool quaternionHeadfreeOffsetSet(void) {
      if ((ABS(attitude.values.roll) < 450)  && (ABS(attitude.values.pitch) < 450)) {
        const float yaw = atan2_approx((+2.0f * (qpAttitude.wz + qpAttitude.xy)), (+1.0f - 2.0f * (qpAttitude.yy + qpAttitude.zz)));
        qOffset.w = cos_approx(yaw/2);
        qOffset.x = 0;
        qOffset.y = 0;
        qOffset.z = sin_approx(yaw/2);
        quaternionInverse(&qOffset, &qOffset);
        return(true);
    } else {
        return false;
    }
}
