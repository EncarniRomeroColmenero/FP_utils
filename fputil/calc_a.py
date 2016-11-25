import collections
import csv
import math
import os
import shutil
import sys

from enum import Enum

from astropy.io import fits
from PyQt4.QtGui import *
from PyQt4.QtCore import *


class ObservationType(Enum):
    arc = 'arc'
    bias = 'bias'
    flat = 'flat'
    science = 'science'
    other = 'other'


class ObservationDetails:
    def __init__(self, observation_type, lamp, filter, binning, gain, readout_speed):
        self.observation_type = observation_type.name
        self.lamp = lamp
        self.filter = filter
        self.binning = binning
        self.gain = gain
        self.readout_speed = readout_speed

    def hashcode(self, observation_type=None):
        """Return a hashcode for these observation details.

        The hashcode is unique in the sense that the hashcode of two ObservationDetails instances will be the same if
        and only if the details are the same for both instances.

        If you pass a non-None observation type, that observation type will be used instead of this instance's one.
        Otherwise the instance's observation type is used.

        Params:
        -------
        observation_type: ObservationType
            The observation type to use. None means that this instance's observation type is used.

        Returns:
        --------
        str:
            The hashcode.
        """

        return '{observation_type}::{lamp}::{filter}::{binning}::{gain}::{readout_speed}' \
            .format(observation_type=self.observation_type,
                    lamp=self.lamp,
                    filter=self.filter,
                    binning=self.binning,
                    gain=self.gain,
                    readout_speed=self.readout_speed)


class FlatSelectionDialog(QDialog):
    """Dialog for selecting a flat file.

    The user may select a file from a combo box, but they may also opt to use no flat. After the dialog has been closed
    the selected file may be retrieved by the get_selected_flat method.

    Params:
    -------
    flats: list of str
        The file paths for the flats to choose from. These can (and should) be absolute paths. The combo box will only
        show the filenames.
    parent: QWidget
        The parent component for this dialog. May be null.
    """

    PLEASE_SELECT = '--- Please select ---'

    NO_FLAT = 'Don\'t choose any flat'

    def __init__(self, flats, parent = None):
        super(FlatSelectionDialog, self).__init__(parent)

        self.flats = [None]
        self.flats.extend(flats)
        self.flats.append(None)
        self.items = [FlatSelectionDialog.PLEASE_SELECT]
        self.items.extend([os.path.basename(flat) for flat in self.flats[1:len(self.flats) - 1]])
        self.items.append(FlatSelectionDialog.NO_FLAT)
        flat_label = QLabel('Flat')
        self.flat_combo_box = QComboBox()
        self.flat_combo_box.addItems(self.items)
        self.error_label = QLabel()
        ok_cancel_button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)

        self.connect(ok_cancel_button_box, SIGNAL('accepted()'), self, SLOT('accept()'))
        self.connect(ok_cancel_button_box, SIGNAL('rejected()'), self, SLOT('reject()'))

        layout = QGridLayout()
        layout.addWidget(flat_label, 0, 0)
        layout.addWidget(self.flat_combo_box, 0, 1)
        layout.addWidget(self.error_label, 1, 0, 1, 2)
        layout.addWidget(ok_cancel_button_box, 2, 0, 1, 2)

        self.setLayout(layout)

    def selected_flat(self):
        """Return the file path of the flat selected by the user.

        This method returns None until the dialog is closed. If no flat has been selected, None is returned. None is
        also returned if the dialog has been cancelled.

        Returns:
        --------
        str:
            The file path of the selected flat.
        """

        return self.flat

    def accept(self):
        current_index = self.flat_combo_box.currentIndex()
        current_item = self.items[current_index]
        if current_item == FlatSelectionDialog.PLEASE_SELECT:
            self.error_label.setText('<font color="red"><b>You must select a flat.</b></font>')
            return
        if current_item != FlatSelectionDialog.NO_FLAT:
            self.flat = self.flats[current_index]
        else:
            self.flat = None
        QDialog.accept(self)


class ACalculator:
    def __init__(self, ring_file=None, output_dir=None, template_dir=None):
        if not ring_file or not os.path.exists(ring_file):
            raise ValueError('The ring file does not exist: {ring_file}'
                             .format(ring_file=ring_file))
        input_dir = os.path.abspath(os.path.join(ring_file, os.pardir))
        if not output_dir or not os.path.exists(output_dir):
            raise ValueError('The output directory does not exist: {output_dir}'.format(output_dir=output_dir))
        if input_dir == output_dir:
            raise ValueError('The output directory must not be the same as the input directory: {output_dir}'
                             .format(output_dir=output_dir))
        self.output_dir = os.path.abspath(output_dir)
        self.ring_file = os.path.abspath(os.path.join(self.output_dir, os.path.basename(ring_file)))
        if not template_dir or not os.path.exists(template_dir):
            raise ValueError('The template directory does not exist: {template_dir}'
                             .format(template_dir=template_dir))
        if input_dir == template_dir:
            raise ValueError('The template directory must not be the same as the input directory: {template_dir}'
                             .format(template_dir=template_dir))
        if output_dir == template_dir:
            raise ValueError('The template directory must not be the same as the output directory: {template_dir}'
                             .format(template_dir=template_dir))
        self.template_dir = os.path.abspath(template_dir)

        self.copy_input_files(input_dir)
        raw_hashcodes_cache = os.path.join(self.output_dir, 'raw_hashcodes.csv')
        self.input_file_hashcodes = self.find_observation_details_hashcodes(self.output_dir, raw_hashcodes_cache)
        template_hashcodes_cache = os.path.join(self.output_dir, 'template_hashcodes.csv')
        self.template_hashcodes = self.find_observation_details_hashcodes(self.template_dir, template_hashcodes_cache)
        self.ignored_files = self.find_ignored_files()
        self.ignored_warnings = self.find_ignored_warnings()

    def copy_input_files(self, input_dir):
        """Copy all the input files to the output directory.

        The input files are all the files located in the given input directory. If a file exists in the output
        directory already it isn't copied again.
        """

        for f in os.listdir(input_dir):
            input_file = os.path.abspath(os.path.join(input_dir, f))
            output_file = os.path.abspath(os.path.join(self.output_dir, f))
            if not os.path.exists(output_file):
                shutil.copy(input_file, output_file)

    def preprocess(self):
        # flat fielding
        flat = self.choose_flat()
        if flat is not None:
            pass  # perform flat fielding
        else:
            no_flat_warning = 'no flat for {hashcode}'.format(hashcode=self.input_file_hashcodes[self.ring_file])
            if no_flat_warning not in self.ignored_warnings:
                self.ignore_warning(no_flat_warning)

    def find_observation_details_hashcodes(self, fits_dir, cache_file):
        """Return the lookup dictionary for the observation detail hashcodes of all the FITS files.

        It is first attempted to read the frame types from the cache file. Any missing values are then obtained from
        the FITS files.

        Params:
        -------
        fits_dir: str
            The directory with the FITS files to consider.
        cache_file: str
            The path of the cache file to use.

        Return:
        -------
        dict:
            The dictionary of file paths and corresponding hashcodes.
        """

        # get cached observation details hashcodes
        observation_details_hashcodes = dict()
        try:
            with open(cache_file, 'r') as cf:
                reader = csv.DictReader(cf)
                for row in reader:
                    observation_details_hashcodes[row['File']] = row['Hashcode']
        except IOError:
            pass

        # get (missing) observation details from the FITS files
        fits_files = [os.path.join(fits_dir, f) for f in os.listdir(fits_dir)
                      if f.lower().endswith('.fits')]
        for f in fits_files:
            if f not in observation_details_hashcodes:
                observation_details_hashcodes[f] = self.observation_details_from_fits(os.path.abspath(f)).hashcode()

        # cache observation details
        with open(cache_file, 'w') as cf:
            writer = csv.DictWriter(cf, ['File', 'Hashcode'])
            writer.writeheader()
            for f in sorted(observation_details_hashcodes.keys()):
                writer.writerow({'File': f,
                                 'Hashcode': observation_details_hashcodes[f]})

        return observation_details_hashcodes

    def find_ignored_files(self):
        """Return the list of ignored files.

        Return:
        -------
        list:
            The ignored files.
        """

        try:
            with open(self.ignored_files_file(), 'r') as f:
                return [l.strip() for l in f.readlines()]
        except IOError:
            return []

    def ignore_files(self, files):
        """Mark a file to be ignored.

        Params:
        -------
        files: list of str
            The path of the file which shall be ignored.
        """

        try:
            with open(self.ignored_files_file(), 'r') as f:
                ignored_files = [l.strip() for l in f.readline()]
        except IOError:
            ignored_files = []

        for f in files:
            if f not in ignored_files:
                ignored_files.append(os.path.abspath(f))

        ignored_files.sort()
        with open(self.ignored_files_file(), 'w') as f:
            for ignored_file in ignored_files:
                f.write(ignored_file + '\n')

    def ignored_files_file(self):
        """Return the path of the file storing the files to ignore.

        Returns:
        --------
        str:
            The file path.
        """

        return os.path.join(self.output_dir, 'ignored_files.txt')

    def ignore_warning(self, warning):
        """Mark a warning as to be ignored.


        Params:
        -------
        warning: str
            A string uniquely identifying the warning to ignore.
        """

        try:
            with open(self.ignored_warnings_file(), 'r') as f:
                ignored_warnings = [l.strip() for l in f.readlines()]
        except IOError:
            ignored_warnings = []

        if warning not in ignored_warnings:
            ignored_warnings.append(warning)

        ignored_warnings.sort()
        with open(self.ignored_warnings_file(), 'w') as f:
            for ignored_warning in ignored_warnings:
                f.write(ignored_warning)
                QMessageBox.warning(None,
                                    "No flat",
                                    "No flat could be found, or you have decided to use no flat. Hence no flat"
                                    "fielding is done.")

    def find_ignored_warnings(self):
        """Return the list of ignored files.

        Return:
        -------
        list:
            The ignored files.
        """

        try:
            with open(self.ignored_warnings_file(), 'r') as f:
                return [l.strip() for l in f.readlines()]
        except IOError:
            return []

    def ignored_warnings_file(self):
        """Return the path of the file storing the warnings to ignore.

        Returns:
        --------
        str:
            The file path.
        """

        return os.path.join(self.output_dir, 'ignored_warnings.txt')

    def observation_details_from_fits(self, f):
        """Obtain the observation details.

        The observation type is found using the following rules.

        1. The observation is an arc if the OBSTYPE header is ARC or FLAT and the ET-STATE header does not start with
           S1.

        2. The observation is a flat if the OBSTYPE header is FLAT and the ET-STATE header starts with S1.

        3. The observation is a bias if the OBSTYPE header is BIAS.

        4. The observation is science if the OBSTYPE header is SCIENCE.

        5. The observation type is Other in all other cases.

        Params:
        -------
        f: str
            Name of the FITS file

        Return:
        -------
        ObservationType:
           File type of the FITS file.
        """

        with fits.open(f) as hdulist:
            # read headers from file
            prihdr = hdulist[0].header

            # observation  type
            et_state_header = prihdr.get('ET-STATE')
            obs_type_header = prihdr.get('OBSTYPE')
            if obs_type_header == 'ARC' or\
                    (obs_type_header == 'FLAT' and et_state_header and not et_state_header.startswith('S1')):
                obs_type = ObservationType.arc
            elif obs_type_header == 'FLAT' and et_state_header and et_state_header.startswith('S1'):
                obs_type = ObservationType.flat
            elif obs_type_header == 'BIAS':
                obs_type = ObservationType.bias
            elif obs_type_header == 'OBJECT':
                obs_type = ObservationType.science
            else:
                obs_type = ObservationType.other

            # lamp
            lamp = prihdr.get('LAMPID')

            # filter
            filter = prihdr.get('FILTER')

            # binning
            binning = prihdr.get('CCDSUM')

            # gain
            gain = prihdr.get('GAINSET')

            # readout speed
            readout_speed = prihdr.get('ROSPEED')

            return ObservationDetails(observation_type=obs_type,
                                      lamp=lamp,
                                      filter=filter,
                                      binning=binning,
                                      gain=gain,
                                      readout_speed=readout_speed)

    def choose_flat(self):
        # search for potential flats
        hashcode = self.observation_details_from_fits(self.ring_file).hashcode(observation_type=ObservationType.flat)
        flats = [f for f, h in self.input_file_hashcodes.items() if h == hashcode and f not in self.ignored_files]
        flats.sort()
        flat = None

        # if there is no flat whatsoever we try the templates
        if len(flats) == 0:
            flats = [f for f, h in self.template_hashcodes.items() if h == hashcode and f not in self.ignored_files]
            flats.sort()

        # all is well if there is only one file
        if len(flats) == 1:
            flat = flats[0]

        # the user has to choose if there is more than one file
        if len(flats) > 1:
            flat_dialog = FlatSelectionDialog(flats)
            if flat_dialog.exec_():
                flat = flat_dialog.selected_flat()
                ignored_flats = [f for f in flats if f != flat]
                self.ignore_files(ignored_flats)
            else:
                sys.exit(1)

        # return the choice
        return flat

aa = QApplication(sys.argv)
a = ACalculator(ring_file='/Users/christian/IdeaProjects/FP_Utils/TestData/raw/P201605090023.fits',
                output_dir='/Users/christian/IdeaProjects/FP_Utils/TestDataOutput',
                template_dir='/Users/christian/IdeaProjects/FP_Utils/TestDataTemplates')
a.preprocess()


def A(ring_file, R, W_line):
    """Calculate the A value for ring.

    A is calculated for the ring parameters in the given FITS file and the given ring radius R as well as emission line
    wavelength W_line.

    Params:
    -------
    ring_file: file path or file object or file-like object
        A file-like object or a file path for the FITS file from which to get ring parameters.
    R: float
        Ring radius in unbinned pixels.
    W_line: float
        Lab wavelength of the the emission line used (in Angstroms).

    Returns:
    --------
    float:
        The value of A.
    """

    # read headers from file
    with fits.open(ring_file) as hdulist:
        prihdr = hdulist[0].header

    missing_card = lambda key: 'Card missing in FITS header: {key}'.format(key=key)

    # make sure an etalon state is defined
    ETALON_STATE = 'ET-STATE'
    if ETALON_STATE not in prihdr:
        raise ValueError(missing_card(ETALON_STATE))

    # parse the etalon state
    etalon_state = prihdr[ETALON_STATE]
    if etalon_state == 'S2 - Etalon 1':
        etalon = 1
    elif etalon_state == 'S3 - Etalon 2':
        etalon = 2
    elif etalon_state == 'S4 - Etalon 1 & 2':
        etalon = 1
    else:
        raise ValueError('Unsupported value for {key} card: {value}'.format(key=ETALON_STATE, value=etalon_state))

    # make sure all ring parameters are defined
    ETALON_B = 'ET{etalon}B'.format(etalon=etalon)
    ETALON_F = 'ET{etalon}F'.format(etalon=etalon)
    ETALON_Z = 'ET{etalon}Z'.format(etalon=etalon)
    for key in [ETALON_B, ETALON_F, ETALON_Z]:
        if key not in prihdr:
            raise ValueError(missing_card(key))

    # collect the ring parameters
    B = prihdr[ETALON_B]
    F = prihdr[ETALON_F]
    Z = prihdr[ETALON_Z]

    # calculate A
    W_central = W_line * math.sqrt(1 + (float(R) / F) ** 2)
    return W_central - B * Z