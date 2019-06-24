# This script is for importing a set of pyControl behavioural data files and 
# pyPhotometry data files contained in seperate data folders.  The script automatically 
# finds correponding data files in each folder based on the subject ID and date.  
# If your subject IDs have a different format from those in the example data you may need
# to modify the function _get_file_info() that extracts the name and date.

# The import code generates a list called 'sessions' which contains behavioural session objects.
# Photometry data attached to each session as the session.photo_data attribute.
# The Rsync_aligner object which is used to convert behavioural event times into 
# photometry sample numbers is stored on each session as session.photo_data['aligner']

# The function average_response() shows a simple analysis using the imported data.

# The only preprocessing performed on the photometry data is bandpass filtering
# which occurs in the import_ppd function.  To add additional pre-processing steps it 
# is recomended to modify the import_ppd function.

# (c) Thomas Akam 2019. Released under the GPL-3 open source licence.

import os
import numpy as np
import pylab as plt

from pyControl_import import Session
from pyPhotometry_import import import_ppd
from rsync import Rsync_aligner, RsyncError

# Specify data folders.

photometry_data_folder = os.path.join('..','data','photometry')
behaviour_data_folder  = os.path.join('..','data','behaviour')

# Get all file names in each folder.

photometry_filenames = os.listdir(photometry_data_folder)
behaviour_filenames  = os.listdir(behaviour_data_folder) 

# Find correponding filenames.

def _get_file_info(filename):
    '''Given a filename return the subject ID and date.'''
    subject_ID = filename[:4]
    date = filename.split('-',1)[1][:10]
    return subject_ID, date

coresponding_files = [] # [(behaviour_filename, photometry_filename)]

for b_file_name in behaviour_filenames:
    b_file_info = _get_file_info(b_file_name)
    p_file_name = next((f_name for f_name in photometry_filenames 
                        if b_file_info == _get_file_info(f_name)), None)
    coresponding_files.append((b_file_name, p_file_name))

# Load data and setup synchronisation.

sessions = []     # List of sessions with attached photometry data.
sync_errors = []  # List of filenames where the sync pulses do not align.

for b_file_name, p_file_name in coresponding_files:
    # Load behavioural data.
    session = Session(os.path.join(behaviour_data_folder, b_file_name))
    if p_file_name: # If a corresponding photometry data file was found.
        # Load photometry data.
        photo_data = import_ppd(os.path.join(photometry_data_folder, p_file_name))
        # Setup Rsync aligner.
        sync_pulse_samples = photo_data['pulse_inds_2'] # Photometry samples when sync pulses occured.
        sync_pulse_times   = session.times['rsync']     # Sync pulse times recorded by behavioural system.
        try: 
            # Instantiate aligner.
            aligner = Rsync_aligner(pulse_times_A=sync_pulse_times , pulse_times_B=sync_pulse_samples, 
                                    units_B=1000/photo_data['sampling_rate'])
        except RsyncError:
            sync_errors.append((b_file_name, p_file_name))
            continue # Go to next file pair.
        # Attach photometry data to session and append to list of good sessions.
        session.photo_data = photo_data
        session.photo_data['aligner'] = aligner
        sessions.append(session)

print('\nMatching behaviour and photometry data found for {} sessions.'.format(len(sessions)))

missing_sessions = [b_file_name for (b_file_name, p_file_name) in coresponding_files
                             if p_file_name is None]

if missing_sessions:
    print('\n No matching photometry data file found for behaviour sessions:')
    for b_file_name in missing_sessions:
        print('    ' + b_file_name)

if sync_errors:
    print('\nSync pulses could not be aligned for sessions:')
    for b_file_name, p_file_name in sync_errors:
        print('    ' + b_file_name + ' : ' + p_file_name)

# ----------------------------------------------------------------------------------------------

def average_response(sessions, event='reward_right_available', window_ms=[-1000, 4000]):
    '''Calculate the average photometry response around the the specified event.
    Arguments:
        sessions  : the list of sessions to analyse.
        event     : the name of the event to extract traces around.
        window_ms : the time window to plot around the event in ms. 
    '''
    traces = []

    for session in sessions:
        window_size = (np.array(window_ms)*session.photo_data['sampling_rate']/1000).astype(int) # Time window around events to plot in samples.
        event_times = session.times[event] # Times when events occured.
        event_samples = session.photo_data['aligner'].A_to_B(event_times).astype(int) # Convert event times to photometry sample numbers.
        for es in event_samples:
            traces.append(session.photo_data['analog_1_filt'][es+window_size[0]:es+window_size[1]])

    traces = np.vstack(traces) # Convert traces to an array.
    mean_trace = np.mean(traces,0)
    rel_time = np.arange(*window_size) * 1000/session.photo_data['sampling_rate'] # Time relative to event

    plt.figure(1).clf()
    plt.plot(rel_time, mean_trace, 'g', linewidth=1.2)
    plt.xlim(*window_ms)
    plt.xlabel('Time relative to {} (ms)'.format(event))    
    plt.ylabel('GCaMP signal (V)')
