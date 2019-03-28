#pragma once

// FIXME some of these are flight modes, some of these are general status indicators
typedef enum {
    ARMED                       = (1 << 0),
    WAS_EVER_ARMED              = (1 << 1),
    WAS_ARMED_WITH_PREARM       = (1 << 2)
} armingFlag_e;

extern uint8_t armingFlags;

#define DISABLE_ARMING_FLAG(mask) (armingFlags &= ~(mask))
#define ENABLE_ARMING_FLAG(mask) (armingFlags |= (mask))
#define ARMING_FLAG(mask) (armingFlags & (mask))

/*
 * Arming disable flags are listed in the order of criticalness.
 * (Beeper code can notify the most critical reason.)
 */
typedef enum {
    ARMING_DISABLED_NO_GYRO         = (1 << 0),
    ARMING_DISABLED_FAILSAFE        = (1 << 1),
    ARMING_DISABLED_RX_FAILSAFE     = (1 << 2),
    ARMING_DISABLED_BAD_RX_RECOVERY = (1 << 3),
    ARMING_DISABLED_BOXFAILSAFE     = (1 << 4),
    ARMING_DISABLED_RUNAWAY_TAKEOFF = (1 << 5),
    ARMING_DISABLED_THROTTLE        = (1 << 6),
    ARMING_DISABLED_ANGLE           = (1 << 7),
    ARMING_DISABLED_BOOT_GRACE_TIME = (1 << 8),
    ARMING_DISABLED_NOPREARM        = (1 << 9),
    ARMING_DISABLED_LOAD            = (1 << 10),
    ARMING_DISABLED_CALIBRATING     = (1 << 11),
    ARMING_DISABLED_CLI             = (1 << 12),
    ARMING_DISABLED_CMS_MENU        = (1 << 13),
    ARMING_DISABLED_OSD_MENU        = (1 << 14),
    ARMING_DISABLED_BST             = (1 << 15),
    ARMING_DISABLED_MSP             = (1 << 16),
    ARMING_DISABLED_PARALYZE        = (1 << 17),
    ARMING_DISABLED_GPS             = (1 << 18),
    ARMING_DISABLED_ARM_SWITCH      = (1 << 19), // Needs to be the last element, since it's always activated if one of the others is active when arming
} armingDisableFlags_e;

#define ARMING_DISABLE_FLAGS_COUNT 20

extern const char *armingDisableFlagNames[ARMING_DISABLE_FLAGS_COUNT];

void setArmingDisabled(armingDisableFlags_e flag);
void unsetArmingDisabled(armingDisableFlags_e flag);
bool getArmingDisabled(armingDisableFlags_e flag);
bool isArmingDisabled(void);
armingDisableFlags_e getArmingDisableFlags(void);

typedef enum {
    ANGLE_MODE      = (1 << 0),
    HORIZON_MODE    = (1 << 1),
    MAG_MODE        = (1 << 2),
    BARO_MODE       = (1 << 3),
    GPS_HOME_MODE   = (1 << 4),
    GPS_HOLD_MODE   = (1 << 5),
    HEADFREE_MODE   = (1 << 6),
//    UNUSED_MODE     = (1 << 7), // old autotune
    PASSTHRU_MODE   = (1 << 8),
//    RANGEFINDER_MODE= (1 << 9),
    FAILSAFE_MODE   = (1 << 10),
    GPS_RESCUE_MODE = (1 << 11)
} flightModeFlags_e;

extern uint16_t flightModeFlags;

#define DISABLE_FLIGHT_MODE(mask) disableFlightMode(mask)
#define ENABLE_FLIGHT_MODE(mask) enableFlightMode(mask)
#define FLIGHT_MODE(mask) (flightModeFlags & (mask))

// macro to initialize map from boxId_e to log2(flightModeFlags). Keep it in sync with flightModeFlags_e enum.
// [BOXARM] is left unpopulated
#define BOXID_TO_FLIGHT_MODE_MAP_INITIALIZER {           \
   [BOXANGLE]       = LOG2(ANGLE_MODE),                  \
   [BOXHORIZON]     = LOG2(HORIZON_MODE),                \
   [BOXMAG]         = LOG2(MAG_MODE),                    \
   [BOXBARO]        = LOG2(BARO_MODE),                   \
   [BOXGPSHOME]     = LOG2(GPS_HOME_MODE),               \
   [BOXGPSHOLD]     = LOG2( GPS_HOLD_MODE),              \
   [BOXHEADFREE]    = LOG2(HEADFREE_MODE),               \
   [BOXPASSTHRU]    = LOG2(PASSTHRU_MODE),               \
   [BOXFAILSAFE]    = LOG2(FAILSAFE_MODE),               \
   [BOXGPSRESCUE]   = LOG2(GPS_RESCUE_MODE),             \
}                                                        \
/**/

typedef enum {
    GPS_FIX_HOME   = (1 << 0),
    GPS_FIX        = (1 << 1),
    CALIBRATE_MAG  = (1 << 2),
    SMALL_ANGLE    = (1 << 3),
    FIXED_WING     = (1 << 4)                    // set when in flying_wing or airplane mode. currently used by althold selection code
} stateFlags_t;

#define DISABLE_STATE(mask) (stateFlags &= ~(mask))
#define ENABLE_STATE(mask) (stateFlags |= (mask))
#define STATE(mask) (stateFlags & (mask))

extern uint8_t stateFlags;

uint16_t enableFlightMode(flightModeFlags_e mask);
uint16_t disableFlightMode(flightModeFlags_e mask);

bool sensors(uint32_t mask);
void sensorsSet(uint32_t mask);
void sensorsClear(uint32_t mask);
uint32_t sensorsMask(void);

void mwDisarm(void);
