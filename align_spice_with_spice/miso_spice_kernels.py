import spiceypy
import datetime
import math
import os

class settings:
    KERNELS_FOLDER = os.environ['SPICE_KERNELS_SOLO']

class SpiceManager:
    # SOLAR ORBITER naif identifier
    solar_orbiter_naif_id = -144 # SOLO identifier: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/naif_ids.html

    def __init__(self, tls_filename, sclk_filename):
        try:
            spiceypy.furnsh(tls_filename)
        except:
            return
        try:
            spiceypy.furnsh(sclk_filename)
        except:
            return
    def obt2utc(self, obt_string):
        # Obt to Ephemeris time (seconds past J2000)
        ephemeris_time = spiceypy.scs2e(self.solar_orbiter_naif_id, obt_string)
        # Ephemeris time to Utc
        # Format of output epoch: ISOC (ISO Calendar format, UTC)
        # Digits of precision in fractional seconds: 3
        return spiceypy.et2utc(ephemeris_time, "ISOC", 3)

    def utc2obt(self, utc_string):
        # Utc to Ephemeris time (seconds past J2000)
        ephemeris_time = spiceypy.utc2et(utc_string)
        # Ephemeris time to Obt
        return spiceypy.sce2s(self.solar_orbiter_naif_id, ephemeris_time)


class SpiceKernels():
    UA_in_km = 149597870
    max_pos_per_object = 1e6
    utc_to_tdb = -69.1853861
    mk_file = None
    spice = None
    spice_manager = None
    kernels_folder = 'SPICE_kernels'

    def __init__(self):
        self.spice = spiceypy
        tls_filename = settings.KERNELS_FOLDER + '/lsk/naif0012.tls'

        # The file.tm loads each kernel file we need
        if self.SPICE_kernels_exists():
            self.mk_file = os.popen('ls ' + settings.KERNELS_FOLDER + '/mk/solo_ANC_soc-flown-mk_*').read().split()[0]
            self.mk_file_bepi = settings.KERNELS_FOLDER + '/BEPICOLOMBO/mk/bc_plan_v260_20201214_001.tm'
            sclk_filename = os.popen('ls -1 ' + settings.KERNELS_FOLDER + '/sclk/solo_ANC_soc-sclk_* | sort -r | head -n 1').read().split()[0]
            self.spice_manager = SpiceManager(tls_filename, sclk_filename)
        else:
            self.mk_file = ''
            self.mk_file_bepi = ''
            self.spice_manager = None

    def SPICE_kernels_exists(self, test_ready=False):
        ready_str = ''
        if test_ready:
            ready_str = '_ready'
        return os.path.exists(settings.KERNELS_FOLDER + ready_str)

    def getMkFile(self, object_id):
        if object_id == 'MPO' or object_id == 'MMO':
            return self.mk_file_bepi

        if self.SPICE_kernels_exists():
            # The path needs to be recovered again, otherwise it can yield a location path error because of an obsolete one
            return os.popen('ls ' + settings.KERNELS_FOLDER + '/mk/solo_ANC_soc-flown-mk_*').read().split()[0]
        return ""

    def getVersion(self):
        if self.SPICE_kernels_exists():
            cmd_ls = 'cd ' + settings.KERNELS_FOLDER + '/mk/ && ls solo_ANC_soc-flown-mk_*'
            result = os.popen(cmd_ls).read().split()[0]
            result = result.split('.')[0]
            result = result.split('-mk_')[1]
            return '"flown" ' + result.replace('_', ' ')
        else:
            settings.LOGGER.error(' SPICE_kernels/ is missing. MISO will run without the sun viewer visualization.')
            return ""

    def thereIsNewVersion(self):
        # check if in SPICE_kernels_ready/ new delta files were loaded (i.e. the mk/flown... file changed)
        if self.SPICE_kernels_exists() and self.SPICE_kernels_exists(True):
            cmd_ls = 'cd ' + settings.KERNELS_FOLDER + '_ready/mk/ && ls solo_ANC_soc-flown-mk_*'  # FIXME
            result = os.popen(cmd_ls).read().split()[0]
            ready_version = result.split('.')[0]
            cmd_ls = 'cd ' + settings.KERNELS_FOLDER + '/mk/ && ls solo_ANC_soc-flown-mk_*'
            result = os.popen(cmd_ls).read().split()[0]
            current_version = result.split('.')[0]
            return ready_version != current_version
        else:
            return False

    def thereIsBackup(self):
        # check if in SPICE_kernels_old/ exists
        return os.path.exists(settings.KERNELS_FOLDER + '_old/')  # FIXME

    def restoreUTC(self, et):
        # convert et in TDB J2000 (spice kernels) to string 1970 UTC
        et += self.utc_to_tdb
        j2000_to_1970 = datetime.datetime(2000, 1, 1, 12) - datetime.datetime(1970, 1, 1)
        dt = datetime.datetime.fromtimestamp(et) + j2000_to_1970
        return dt.strftime("%Y-%m-%d %H:%M:%S")

    def getOrbit(self, object_id, timestamps, obs_id='SUN'):
        ''' Get the position of an object in the HCI frame

        Parameters
        ==========
        object_id : str
            Target body name. Can be 'EARTH', 'SOLO', ...
        timestamps : list of str
            List of timestamps at which to get the orbit. These timestamps are
            strings describing an epoch, parsable by spiceypy.str2et (e.g.
            YYYY-MM-DDTHH:MM:SS.SSS)..
        obs_id : str (default: 'SUN')
            The reference observing body name.

        Returns
        =======
        positions : list of dict
            List of positions for the object in the HCI frame, in km.
            Each item of the list is a dict with keys 'x', 'y', and 'z'.
        '''
        if self.SPICE_kernels_exists():
            self.clear_session()
            self.spice.furnsh(self.getMkFile(object_id))

            points = []

            for et in self.spice.str2et(timestamps):
                points_km, lightTimes = self.spice.spkpos(object_id, et, 'HCI', 'NONE', obs_id)

                pos_x_y_z = {'x': points_km[0], 'y': points_km[2], 'z': -points_km[1]}
                points.append(pos_x_y_z)

            self.clear_session()

            return points
        else:
            return [{}]

    def getObjectPositions(self, object_id, timestamps, obs_id='SUN'):
        """ Get the positon and rotation of an object in the HCI frame.

        Parameters
        ==========
        object_id : str
            Target body name. Can be 'EARTH', 'SOLO', ...
        timestamps : list of str
            List of timestamps at which to get the orbit. These timestamps are
            strings describing an epoch, parsable by spiceypy.str2et (e.g.
            YYYY-MM-DDTHH:MM:SS.SSS)..
        obs_id : str (default: 'SUN')
            The reference observing body name.

        Returns
        =======
        positions : list of dict
            List of positions for the object in the HCI frame, in km.
            Each item of the list is a dict with keys 'x', 'y', and 'z'.
        rotations : list of dict
            List of rotation vectors for the object.
            Each item of the list is a dict with keys 'x', 'y', and 'z'.
        ecliptic :
            If the object is 'EARTH', return a list of rotations of the
            ecliptic.
        """
        positions = []
        rotations = []
        ecliptic = []
        if not self.SPICE_kernels_exists():
            return positions, rotations, ecliptic

        apply_rotations_on = ['SUN', 'EARTH']
        self.clear_session()
        self.spice.furnsh(self.getMkFile(object_id))

        for et in self.spice.str2et(timestamps):
            position_km, lightTimes = self.spice.spkpos(object_id, et, 'HCI', 'NONE', obs_id)

            pos_x_y_z = {'x': position_km[0], 'y': position_km[2], 'z': -position_km[1]}
            positions.append(pos_x_y_z)

            if object_id == 'SOLO':
                rotations_vector = self.getSoloRotation(et)
            elif object_id in apply_rotations_on:
                rotations_vector = self.getObjectRotation('IAU_' + object_id, et)
            else:
                rotations_vector = [0., 0., 0.]

            rotations.append(rotations_vector)

            if object_id == 'EARTH': # we don't want other planets
                ecliptic.append(self.getEarthEclipticPlane(et))

        self.clear_session()

        return positions, rotations, ecliptic

    def clear_session(self):
        self.spice.kclear()

    def getObjectRotation(self, frame_name, et):
        axis = [1., 0., 0.]
        try:
            pform = self.spice.pxform(frame_name, "HCI", et)
            axis_result = self.spice.mxv(pform, axis)
        except:
            axis_result = [0, 0, 0]

        return {'x': axis_result[0], 'y': axis_result[2], 'z': -axis_result[1]}

    def getSoloRotation(self, et = None):
        if et is None: # set manually the date, used only from debug.py
            self.spice.furnsh(self.mk_file)
            utc = "January 20, 2020"
            et = self.spice.str2et(utc)

        x = 1.
        z = 0.
        # make an unitary vector:
        # 1 = sqrt(x² + y² + z²)
        # 1 = x² + y² + z²
        # y = sqrt(1 - x² - z²)
        y = math.sqrt(1. - x*x - z*z)

        solo_axis = [x, y, z]

        try:
            pform = self.spice.pxform("SOLO_SRF", "HCI", et)
            new_solo_axis = self.spice.mxv(pform, solo_axis)
            solo_final_axis = {'x': new_solo_axis[0], 'y': new_solo_axis[2], 'z': -new_solo_axis[1]}
        except:
            try:
                pform = self.spice.pxform("SOLO_SUN_RTN", "HCI", et)
                # example of working frames: HCI, IAU_SUN, IAU_EARTH, SOLO_SRF, SOLO_HCI

                new_solo_axis = self.spice.mxv(pform, solo_axis)

                solo_final_axis = {'x': new_solo_axis[0], 'y': new_solo_axis[2], 'z': -new_solo_axis[1]}
            except:
                solo_final_axis = {'x': 0., 'y': 0., 'z': 0.}

        return solo_final_axis

    def getEarthEclipticPlane(self, et):
        try:
            pform = self.spice.pxform("HCI", "HEE", et)
            # HEE: Heliocentric Earth Ecliptic.
            # X is the Sun-Earth line
            # Z is the north pole for the ecliptic of date.
            ecliptic_axis_x = self.spice.mxv(pform, [1., 0., 0.])
            ecliptic_axis_y = self.spice.mxv(pform, [0., 1., 0.])
            ecliptic_axis_x = {'x': -ecliptic_axis_x[0], 'y': ecliptic_axis_x[2], 'z': -ecliptic_axis_x[1]}
            ecliptic_axis_y = {'x': -ecliptic_axis_y[0], 'y': ecliptic_axis_y[2], 'z': -ecliptic_axis_y[1]}

            return [ecliptic_axis_x, ecliptic_axis_y] # y will be -z in ThreeJS
        except:
            return [{'x': 0., 'y': 0., 'z': 0.}, {'x': 0., 'y': 0., 'z': 0.}]

    # Examples / Debugs ------------------------------------------------------------------------------------------------
    def ckgp_test(self, et):
        solo_id = self.spice_manager.solar_orbiter_naif_id
        self.spice.furnsh(self.mk_file)
        # need to convert "et" (J2000) to SCLK (encoded spacecraft clock) with sce2c()
        clock_sclk = self.spice.sce2c(solo_id, et)
        time_tolerance = self.spice.sctiks(solo_id, "10:10")
        self.spice.ckgp(solo_id, clock_sclk, time_tolerance, 'HCI')

    def showCassiniPosition(self):
        print('Example with Cassini position:')
        self.spice.furnsh(settings.KERNELS_FOLDER + "/cassini_example/cassMetaK.txt")
        step = 4000
        utc = ['Jun 20, 2004', 'Dec 1, 2005']
        etOne = self.spice.str2et(utc[0])
        etTwo = self.spice.str2et(utc[1])
        times = [x * (etTwo - etOne) / step + etOne for x in range(step)]
        positions, lightTimes = self.spice.spkpos('Cassini', times, 'J2000', 'NONE', 'SATURN BARYCENTER')
        print(positions[0])
        print(positions[1])
        print("xform:")
        xform = self.spice.pxform('CASSINI_HGA', 'J2000', self.spice.str2et("2004-05-01 00:00:00"))
        print(xform)
        print('----- End of Example with Cassini position\n')

    def showConversions(self):
        # Execute conversions
        obt = "0"
        utc = "2000-01-01T00:00:00.000"
        print("OBT {0} -> UTC {1}".format(obt, self.spice_manager.obt2utc(obt)))
        # Returns OBT 0 -> UTC 2000-01-01T00:00:00.000
        print("UTC {0} -> OBT {1}".format(utc, self.spice_manager.utc2obt(utc)))
        # Returns UTC 2000-01-01T00:00:00.000 -> OBT 1/0000000000:00000
        obt = "521651623:37539"
        utc = "2016-194T15:13:46.381"
        print("OBT {0} -> UTC {1}".format(obt, self.spice_manager.obt2utc(obt)))
        # Returns: OBT 521651623:37539 -> UTC 2016-07-12T15:13:46.381
        print("UTC {0} -> OBT {1}".format(utc, self.spice_manager.utc2obt(utc)))
        # Returns: UTC 2016-194T15:13:46.381 -> OBT 1/0521651623:3753
