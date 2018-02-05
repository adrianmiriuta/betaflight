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

#define IMU_LOCK pthread_mutex_unlock(&imuUpdateLock)
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
STATIC_UNIT_TESTED quaternion qMahonyAHRS = QUATERNION_INITIALIZE;
STATIC_UNIT_TESTED quaternion qGyroAHRS = QUATERNION_INITIALIZE;
STATIC_UNIT_TESTED quaternionProducts qP = QUATERNION_PRODUCTS_INITIALIZE;
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
static void quaternionTransformVectorBodyToEarth(t_fp_vector * v) {
  const float x = (1.0f - 2.0f * qP.yy - 2.0f * qP.zz) * v->V.X + (2.0f * (qP.xy + -qP.wz)) * v->V.Y + (2.0f * (qP.xz - -qP.wy)) * v->V.Z;
  const float y = (2.0f * (qP.xy - -qP.wz)) * v->V.X + (1.0f - 2.0f * qP.xx - 2.0f * qP.zz) * v->V.Y + (2.0f * (qP.yz + -qP.wx)) * v->V.Z;
  const float z = (2.0f * (qP.xz + -qP.wy)) * v->V.X + (2.0f * (qP.yz - -qP.wx)) * v->V.Y + (1.0f - 2.0f * qP.xx - 2.0f * qP.yy) * v->V.Z;

  v->V.X = x;
  v->V.Y = -y;
  v->V.Z = z;
}

// rotate acc into Earth frame and calculate acceleration in it
static void imuCalculateAcceleration(uint32_t deltaT)
{
    static float accZoffset = 0;
    static float accz_smooth = 0;

    // deltaT is measured in us ticks
    const float dT = (float)deltaT * 1e-6f;

    t_fp_vector accel_ned;
    accel_ned.V.X = acc.accADC[X];
    accel_ned.V.Y = acc.accADC[Y];
    accel_ned.V.Z = acc.accADC[Z];

    quaternionTransformVectorBodyToEarth(&accel_ned);

    if (imuRuntimeConfig.acc_unarmedcal == 1) {
        if (!ARMING_FLAG(ARMED)) {
            accZoffset -= accZoffset / 64;
            accZoffset += accel_ned.V.Z;
        }
        accel_ned.V.Z -= accZoffset / 64;  // compensate for gravitation on z-axis
    } else {
        accel_ned.V.Z -= acc.dev.acc_1G;
    }

    accz_smooth = accz_smooth + (dT / (fc_acc + dT)) * (accel_ned.V.Z - accz_smooth); // low pass filter

    // apply Deadband to reduce integration drift and vibration influence
    accSum[X] += applyDeadband(lrintf(accel_ned.V.X), imuRuntimeConfig.accDeadband.xy);
    accSum[Y] += applyDeadband(lrintf(accel_ned.V.Y), imuRuntimeConfig.accDeadband.xy);
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
    static float integralFBx = 0.0f,  integralFBy = 0.0f, integralFBz = 0.0f;    // integral error terms scaled by Ki

    // Calculate general spin rate (rad/s)
    const float spin_rate = sqrtf(sq(gx) + sq(gy) + sq(gz));

    // Use raw heading error (from GPS or whatever else)
    float ez = 0;
    if (useYaw) {
        while (yawError >  M_PIf) yawError -= (2.0f * M_PIf);
        while (yawError < -M_PIf) yawError += (2.0f * M_PIf);

        ez += sin_approx(yawError / 2.0f);
    }

    // Use measured magnetic field vector
    float ex = 0, ey = 0;
    float recipNorm = sq(mx) + sq(my) + sq(mz);
    if (useMag && recipNorm > 0.01f) {
        // Normalise magnetometer measurement
        recipNorm = invSqrt(recipNorm);
        mx *= recipNorm;
        my *= recipNorm;
        mz *= recipNorm;

        // For magnetometer correction we make an assumption that magnetic field is perpendicular to gravity (ignore Z-component in EF).
        // This way magnetic field will only affect heading and wont mess roll/pitch angles

        // (hx; hy; 0) - measured mag field vector in EF (assuming Z-component is zero)
        // (bx; 0; 0) - reference mag field vector heading due North in EF (assuming Z-component is zero)
        const float hx = (1.0f - 2.0f * qP.yy - 2.0f * qP.zz) * mx + (2.0f * (qP.xy + -qP.wz)) * my + (2.0f * (qP.xz - -qP.wy)) * mz;
        const float hy = (2.0f * (qP.xy - -qP.wz)) * mx + (1.0f - 2.0f * qP.xx - 2.0f * qP.zz) * my + (2.0f * (qP.yz + -qP.wx)) * mz;
        const float bx = sqrtf(hx * hx + hy * hy);

        // magnetometer error is cross product between estimated magnetic north and measured magnetic north (calculated in EF)
        const float ez_ef = -(hy * bx);

        // Rotate mag error vector back to BF and accumulate
        ex += (2.0f * (qP.xz + -qP.wy)) * ez_ef;
        ey += (2.0f * (qP.yz - -qP.wx)) * ez_ef;
        ez += (1.0f - 2.0f * qP.xx - 2.0f * qP.yy) * ez_ef;
    }

    // Use measured acceleration vector
    recipNorm = sq(ax) + sq(ay) + sq(az);
    if (useAcc && recipNorm > 0.01f) {
        // Normalise accelerometer measurement
        recipNorm = invSqrt(recipNorm);
        ax *= recipNorm;
        ay *= recipNorm;
        az *= recipNorm;

        // Error is sum of cross product between estimated direction and measured direction of gravity
        ex += (ay * (1.0f - 2.0f * qP.xx - 2.0f * qP.yy) - az * (2.0f * (qP.yz - -qP.wx)));
        ey += (az * (2.0f * (qP.xz + -qP.wy)) - ax * (1.0f - 2.0f * qP.xx - 2.0f * qP.yy));
        ez += (ax * (2.0f * (qP.yz - -qP.wx)) - ay * (2.0f * (qP.xz + -qP.wy)));
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
    }
    else {
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

    // test new method
    quaternion qGyro;
    qGyro.w = cos_approx((gx + gy + gz) * 0.5f * dt);
    qGyro.x = sin_approx(gx * 0.5f * dt);
    qGyro.y = sin_approx(gy * 0.5f * dt);
    qGyro.z = sin_approx(gz * 0.5f * dt);
    quaternionMultiply(&qMahonyAHRS, &qGyro, &qMahonyAHRS);
    quaternionMultiply(&qGyroAHRS, &qGyro, &qGyroAHRS);

    // Ok old bf method
    /*
    quaternion qBuff, qGyro;
    qGyro.w = 0;
    qGyro.x = gx * 0.5f * dt;
    qGyro.y = gy * 0.5f * dt;
    qGyro.z = gz * 0.5f * dt;
    quaternionMultiply(&qMahonyAHRS, &qGyro, &qBuff);
    quaternionAdd(&qMahonyAHRS, &qBuff, &qMahonyAHRS);*/

    // test only GyroAHRS in acro mode
    if ((!FLIGHT_MODE(ANGLE_MODE) && (!FLIGHT_MODE(HORIZON_MODE)))) {
      quaternionNormalize(&qGyroAHRS);
      quaternionComputeProducts(&qGyroAHRS, &qP);
    } else {
      quaternionNormalize(&qMahonyAHRS);
      quaternionComputeProducts(&qMahonyAHRS, &qP);
    }



}

STATIC_UNIT_TESTED void imuUpdateEulerAngles(void){
    quaternionProducts buffer;

    if (FLIGHT_MODE(HEADFREE_MODE)) {
      quaternionComputeProducts(&qHeadfree, &buffer);
    } else {
      if ((!FLIGHT_MODE(ANGLE_MODE) && (!FLIGHT_MODE(HORIZON_MODE)))) {
        quaternionComputeProducts(&qGyroAHRS, &buffer);
      } else {
        quaternionComputeProducts(&qMahonyAHRS, &buffer);
      }
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

static bool isMagnetometerHealthy(void)
{
    return (mag.magADC[X] != 0) && (mag.magADC[Y] != 0) && (mag.magADC[Z] != 0);
}

static void imuCalculateEstimatedAttitude(timeUs_t currentTimeUs)
{
    static uint32_t previousIMUUpdateTime;
    float rawYawError = 0;
    bool useAcc = false;
    bool useMag = false;
    bool useYaw = false;

    uint32_t deltaT = currentTimeUs - previousIMUUpdateTime;
    previousIMUUpdateTime = currentTimeUs;

    if (imuIsAccelerometerHealthy()) {
        useAcc = true;
    }

    if (sensors(SENSOR_MAG) && isMagnetometerHealthy()) {
        useMag = true;
    }
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
    return (1.0f - 2.0f * (qP.xx + qP.yy));
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

    qMahonyAHRS.w = w;
    qMahonyAHRS.x = x;
    qMahonyAHRS.y = y;
    qMahonyAHRS.z = z;

    imuComputeRotationMatrix();
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

void quaternionComputeProducts(quaternion *quat, quaternionProducts *quatProd) {
    quatProd->ww = quat->w * quat->w;
    quatProd->wx = quat->w * quat->x;
    quatProd->wy = quat->w * quat->y;
    quatProd->wz = quat->w * quat->z;
    quatProd->xx = quat->x * quat->x;
    quatProd->xy = quat->x * quat->y;
    quatProd->xz = quat->x * quat->z;
    quatProd->yy = quat->y * quat->y;
    quatProd->yz = quat->y * quat->z;
    quatProd->zz = quat->z * quat->z;
}

bool imuQuaternionHeadfreeOffsetSet(void) {

  if ((!FLIGHT_MODE(ANGLE_MODE) && (!FLIGHT_MODE(HORIZON_MODE)))) {
     quaternionCopy(&qMahonyAHRS, &qOffset);
     quaternionInverse(&qOffset, &qOffset);
     return(true);
  } else {

    if ((ABS(attitude.values.roll) < 450)  && (ABS(attitude.values.pitch) < 450)) {
        //const float yaw = -atan2_approx((+2.0f * (qP.wz + qP.xy)), (+1.0f - 2.0f * (qP.yy + qP.zz)));
        const float yaw = atan2_approx((+2.0f * (qP.wz + qP.xy)), (+1.0f - 2.0f * (qP.yy + qP.zz)));

        qOffset.w = cos_approx(yaw/2);
        qOffset.x = 0;
        qOffset.y = 0;
        qOffset.z = sin_approx(yaw/2);

        quaternionInverse(&qOffset, &qOffset);

        return(true);
    } else {
        return(false);
    }
  }
}

void imuQuaternionHeadfreeTransformVectorEarthToBody(t_fp_vector_def *v) {
    quaternionProducts buffer;

    if ((!FLIGHT_MODE(ANGLE_MODE) && (!FLIGHT_MODE(HORIZON_MODE)))) {
      quaternionMultiply(&qOffset, &qGyroAHRS, &qHeadfree);
    } else {
      quaternionMultiply(&qOffset, &qMahonyAHRS, &qHeadfree);
    }

    quaternionComputeProducts(&qHeadfree, &buffer);

    const float x = (buffer.ww + buffer.xx - buffer.yy - buffer.zz) * v->X + 2.0f * (buffer.xy + buffer.wz) * v->Y + 2.0f * (buffer.xz - buffer.wy) * v->Z;
    const float y = 2.0f * (buffer.xy - buffer.wz) * v->X + (buffer.ww - buffer.xx + buffer.yy - buffer.zz) * v->Y + 2.0f * (buffer.yz + buffer.wx) * v->Z;
    const float z = 2.0f * (buffer.xz + buffer.wy) * v->X + 2.0f * (buffer.yz - buffer.wx) * v->Y + (buffer.ww - buffer.xx - buffer.yy + buffer.zz) * v->Z;

    v->X = x;
    v->Y = y;
    v->Z = z;
}

void quaternionMultiply(quaternion *l, quaternion *r, quaternion *o) {
    const float w = l->w * r->w - l->x * r->x - l->y * r->y - l->z * r->z;
    const float x = l->x * r->w + l->w * r->x + l->y * r->z - l->z * r->y;
    const float y = l->w * r->y - l->x * r->z + l->y * r->w + l->z * r->x;
    const float z = l->w * r->z + l->x * r->y - l->y * r->x + l->z * r->w;
    o->w = w;
    o->x = x;
    o->y = y;
    o->z = z;
}

void quaternionNormalize(quaternion *q) {
    float norm = sqrtf(q->w * q->w + q->x * q->x + q->y * q->y + q->z * q->z);
    if (norm == 0) {
      norm = 0.0000001;
    }
    q->w /= norm;
    q->x /= norm;
    q->y /= norm;
    q->z /= norm;
}

void quaternionAdd(quaternion *l, quaternion *r, quaternion *o) {
    o->w = l->w + r->w;
    o->x = l->x + r->x;
    o->y = l->y + r->y;
    o->z = l->z + r->z;
}

void quaternionCopy(quaternion *s, quaternion *d) {
    d->w = s->w;
    d->x = s->x;
    d->y = s->y;
    d->z = s->z;
}

void quaternionInverse(quaternion *i, quaternion *o) {
    float norm_squared = i->w * i->w + i->x * i->x + i->y * i->y + i->z * i->z;
    if (norm_squared == 0) {
      norm_squared = 0.0000001;
    }
    o->w = i->w / norm_squared;
    o->x = i->x * -1 / norm_squared;
    o->y = i->y * -1 / norm_squared;
    o->z = i->z * -1 / norm_squared;
}
