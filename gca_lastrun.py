#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This experiment was created using PsychoPy3 Experiment Builder (v2022.2.4),
    on Oktober 31, 2023, at 14:19
If you publish work using this script the most relevant publication is:

    Peirce J, Gray JR, Simpson S, MacAskill M, Höchenberger R, Sogo H, Kastman E, Lindeløv JK. (2019) 
        PsychoPy2: Experiments in behavior made easy Behav Res 51: 195. 
        https://doi.org/10.3758/s13428-018-01193-y

"""

# --- Import packages ---
from psychopy import locale_setup
from psychopy import prefs
prefs.hardware['audioLib'] = 'ptb'
from psychopy import sound, gui, visual, core, data, event, logging, clock, colors, layout, iohub, hardware, parallel
from psychopy.constants import (NOT_STARTED, STARTED, PLAYING, PAUSED,
                                STOPPED, FINISHED, PRESSED, RELEASED, FOREVER)

import numpy as np  # whole numpy lib is available, prepend 'np.'
from numpy import (sin, cos, tan, log, log10, pi, average,
                   sqrt, std, deg2rad, rad2deg, linspace, asarray)
from numpy.random import random, randint, normal, shuffle, choice as randchoice
import os  # handy system and path functions
import sys  # to get file system encoding

import psychopy.iohub as io
from psychopy.hardware import keyboard

# Run 'Before Experiment' code from codeTrial
feedback_color = (0,0,0)
feedback_opacity = 0
feedback_audio = "audio/silence.wav"
cursorcolor="white"
port_msg = 0
log_msg = ""

rectsize = (0.2, 0.2)
imagesize_test = [1, 1]

score = 0
points = 0

looked_at = False
paused = False
# Run 'Before Experiment' code from codeFeedback
soundFeedback = sound.Sound('A', secs=1, stereo=True, hamming=True, name='soundFeedback')


# Ensure that relative paths start from the same directory as this script
_thisDir = os.path.dirname(os.path.abspath(__file__))
os.chdir(_thisDir)
# Store info about the experiment session
psychopyVersion = '2022.2.4'
expName = 'gca'  # from the Builder filename that created this script
expInfo = {
    'participant': '',
    'group': ['A', 'B', 'C', 'D'],
}
# --- Show participant info dialog --
dlg = gui.DlgFromDict(dictionary=expInfo, sortKeys=False, title=expName)
if dlg.OK == False:
    core.quit()  # user pressed cancel
expInfo['date'] = data.getDateStr()  # add a simple timestamp
expInfo['expName'] = expName
expInfo['psychopyVersion'] = psychopyVersion

# Data file name stem = absolute path + name; later add .psyexp, .csv, .log, etc
filename = _thisDir + os.sep + u'data/%s_%s_%s' % (expName, expInfo['participant'], expInfo['date'])

# An ExperimentHandler isn't essential but helps with data saving
thisExp = data.ExperimentHandler(name=expName, version='',
    extraInfo=expInfo, runtimeInfo=None,
    originPath='C:\\Users\\sag22id\\Documents\\Projects\\GCA\\gaze_avoidance\\gca_lastrun.py',
    savePickle=True, saveWideText=True,
    dataFileName=filename)
# save a log file for detail verbose info
logFile = logging.LogFile(filename+'.log', level=logging.INFO)
logging.console.setLevel(logging.WARNING)  # this outputs to the screen, not a file

endExpNow = False  # flag for 'escape' or other condition => quit the exp
frameTolerance = 0.001  # how close to onset before 'same' frame

# Start Code - component code to be run after the window creation

# --- Setup the Window ---
win = visual.Window(
    size=[2560, 1440], fullscr=True, screen=1, 
    winType='pyglet', allowStencil=False,
    monitor='labMonitor', color=[0,0,0], colorSpace='rgb',
    blendMode='avg', useFBO=True, 
    units='norm')
win.mouseVisible = False
# store frame rate of monitor if we can measure it
expInfo['frameRate'] = win.getActualFrameRate()
if expInfo['frameRate'] != None:
    frameDur = 1.0 / round(expInfo['frameRate'])
else:
    frameDur = 1.0 / 60.0  # could not measure, so guess
# --- Setup input devices ---
ioConfig = {}

# Setup eyetracking
ioConfig['eyetracker.hw.sr_research.eyelink.EyeTracker'] = {
    'name': 'tracker',
    'model_name': 'EYELINK 1000 TOWER',
    'simulation_mode': False,
    'network_settings': '100.1.1.1',
    'default_native_data_file_name': 'EXPFILE',
    'runtime_settings': {
        'sampling_rate': 1000.0,
        'track_eyes': 'RIGHT_EYE',
        'sample_filtering': {
            'sample_filtering': 'FILTER_LEVEL_2',
            'elLiveFiltering': 'FILTER_LEVEL_OFF',
        },
        'vog_settings': {
            'pupil_measure_types': 'PUPIL_AREA',
            'tracking_mode': 'PUPIL_CR_TRACKING',
            'pupil_center_algorithm': 'CENTROID_FIT',
        }
    }
}

# Setup iohub keyboard
ioConfig['Keyboard'] = dict(use_keymap='psychopy')

ioSession = '1'
if 'session' in expInfo:
    ioSession = str(expInfo['session'])
ioServer = io.launchHubServer(window=win, experiment_code='gca', session_code=ioSession, datastore_name=filename, **ioConfig)
eyetracker = ioServer.getDevice('tracker')

# create a default keyboard (e.g. to check for escape)
defaultKeyboard = keyboard.Keyboard(backend='iohub')

# --- Initialize components for Routine "welcome" ---
textStart = visual.TextStim(win=win, name='textStart',
    text='Sehr geehrter Teilnehmer, sehr geehrte Teilnehmerin,\nvielen Dank für Ihre Teilnahme an unserem Experiment.\n\nBitte drücken Sie die Leertaste, um zu starten.',
    font='Open Sans',
    pos=(0, 0), height=0.06, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=0.0);
spaceStart = keyboard.Keyboard()

# --- Initialize components for Routine "startPainThreshold" ---
textPainThreshold = visual.TextStim(win=win, name='textPainThreshold',
    text='Wir möchten zuerst den elektrischen Reiz kalibrieren.\n\nBitte drücken Sie die Leertaste, um zu starten.',
    font='Open Sans',
    pos=(0, 0), height=0.06, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=0.0);
spacePainThreshold = keyboard.Keyboard()

# --- Initialize components for Routine "painThreshold" ---
textRateT = visual.TextStim(win=win, name='textRateT',
    text='Wie schmerzhaft fanden Sie den elektrischen Reiz?',
    font='Open Sans',
    pos=(0, 0.6), height=0.06, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=0.0);
text_noPainT = visual.TextStim(win=win, name='text_noPainT',
    text='nichts\ngespürt',
    font='Open Sans',
    pos=(-0.5, 0.15), height=0.05, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-1.0);
text_slightPainT = visual.TextStim(win=win, name='text_slightPainT',
    text='eben \nwahrnehmbarer\nSchmerz',
    font='Open Sans',
    pos=(-0.1, 0.15), height=0.05, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-2.0);
text_highPainT = visual.TextStim(win=win, name='text_highPainT',
    text='unerträglicher\nSchmerz',
    font='Open Sans',
    pos=(0.5, 0.15), height=0.05, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-3.0);
sliderPainT = visual.Slider(win=win, name='sliderPainT',
    startValue=None, size=(1.0, 0.1), pos=(0, 0), units=None,
    labels=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10), ticks=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10), granularity=1.0,
    style='rating', styleTweaks=('labels45',), opacity=None,
    labelColor='White', markerColor='Red', lineColor='White', colorSpace='rgb',
    font='Open Sans', labelHeight=0.03,
    flip=False, ori=0.0, depth=-4, readOnly=False)
keys_PainT = keyboard.Keyboard()
# Run 'Begin Experiment' code from codePainRatingT
skipShock = False

# --- Initialize components for Routine "shock" ---
portShock1T = parallel.ParallelPort(address='0x0378')
portShock2T = parallel.ParallelPort(address='0x0378')
portShock3T = parallel.ParallelPort(address='0x0378')
fixcrossT = visual.TextStim(win=win, name='fixcrossT',
    text='+',
    font='Open Sans',
    pos=(0, 0), height=0.1, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-4.0);

# --- Initialize components for Routine "startETCalibration" ---
text_ETCalibration = visual.TextStim(win=win, name='text_ETCalibration',
    text='Wir machen nun mit der Kalibrierung des Eye-Trackers weiter.\n\nBitte drücken Sie wieder die Leertaste, um diese zu starten.',
    font='Open Sans',
    pos=(0, 0), height=0.06, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=0.0);
spaceETCalibration = keyboard.Keyboard()

# --- Initialize components for Routine "startExp" ---
textStartExp = visual.TextStim(win=win, name='textStartExp',
    text='Die Kalibrierung ist abgeschlossen.\n\nSie werden während des Experiments mehrere Bilder sehen.\nBitte geben Sie nun an, wie Sie sich bei der Betrachtung der Bilder fühlen.\n\nDrücken Sie die Leertaste, um zu starten.',
    font='Open Sans',
    pos=(0, 0), height=0.06, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=0.0);
spaceStartExp = keyboard.Keyboard()

# --- Initialize components for Routine "stimRating" ---
textRateStim = visual.TextStim(win=win, name='textRateStim',
    text='Wie fühlen Sie sich bei der Betrachtung des Bildes?',
    font='Open Sans',
    pos=(0, 0.6), height=0.06, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-1.0);
imageRating = visual.ImageStim(
    win=win,
    name='imageRating', 
    image='sin', mask=None, anchor='center',
    ori=0.0, pos=[0,0], size=None,
    color=[1,1,1], colorSpace='rgb', opacity=None,
    flipHoriz=False, flipVert=False,
    texRes=128.0, interpolate=True, depth=-2.0)
textUnpleasant = visual.TextStim(win=win, name='textUnpleasant',
    text='sehr \nunwohl',
    font='Open Sans',
    pos=(-0.5, -0.15), height=0.05, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-3.0);
textPleasant = visual.TextStim(win=win, name='textPleasant',
    text='sehr\nwohl',
    font='Open Sans',
    pos=(0.5, -0.15), height=0.05, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-4.0);
sliderStim = visual.Slider(win=win, name='sliderStim',
    startValue=None, size=(1.0, 0.1), pos=(0, -0.3), units=None,
    labels=(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), ticks=(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), granularity=1.0,
    style='rating', styleTweaks=('labels45',), opacity=None,
    labelColor='White', markerColor='Red', lineColor='White', colorSpace='rgb',
    font='Open Sans', labelHeight=0.03,
    flip=False, ori=0.0, depth=-5, readOnly=False)
textSpaceStim = visual.TextStim(win=win, name='textSpaceStim',
    text='Bitte drücken Sie die Leertaste, um die Antwort zu bestätigen.',
    font='Open Sans',
    pos=(0, -0.6), height=0.06, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-6.0);
spaceStim = keyboard.Keyboard()

# --- Initialize components for Routine "startTask" ---
textStartTask = visual.TextStim(win=win, name='textStartTask',
    text='Wir starten nun mit dem Experiment.\n\nSie werden in jedem Durchgang in einer der vier Ecken des Bildschirms ein Bild präsentiert bekommen. Über Ihr Blickverhalten können Sie eine Belohnung in Form von Punkten erhalten. Allerdings lauert auch die Gefahr eines elektrischen Reizes. Ihr Ziel ist es Ihren Belohnungsscore zu maximieren und die elektrischen Reize zu vermeiden.\n\nBitte fixieren Sie am Anfang jedes Durchgangs das Fixationskreuz.\n\nDrücken Sie die Leertaste, um zu starten.',
    font='Open Sans',
    pos=(0, 0), height=0.06, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=0.0);
spaceStartTask = keyboard.Keyboard()

# --- Initialize components for Routine "blank" ---
blankScreen = visual.TextStim(win=win, name='blankScreen',
    text=None,
    font='Open Sans',
    pos=(0, 0), height=0.05, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=0.0);

# --- Initialize components for Routine "cross" ---
fixcross = visual.TextStim(win=win, name='fixcross',
    text='+',
    font='Open Sans',
    pos=(0, 0), height=0.1, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-1.0);

# --- Initialize components for Routine "blank" ---
blankScreen = visual.TextStim(win=win, name='blankScreen',
    text=None,
    font='Open Sans',
    pos=(0, 0), height=0.05, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=0.0);

# --- Initialize components for Routine "trial" ---
image = visual.ImageStim(
    win=win,
    name='image', 
    image='sin', mask=None, anchor='center',
    ori=0.0, pos=[0,0], size=None,
    color=[1,1,1], colorSpace='rgb', opacity=None,
    flipHoriz=False, flipVert=False,
    texRes=128.0, interpolate=True, depth=-1.0)
roi = visual.ROI(win, name='roi', device=eyetracker,
    debug=False,
    shape='rectangle',
    pos=[0,0], size=1.0, anchor='center', ori=0.0)
gazeCursor = visual.ShapeStim(
    win=win, name='gazeCursor', vertices='star7',
    size=(0.02, 0.02),
    ori=0.0, pos=[0,0], anchor='center',
    lineWidth=1.0,     colorSpace='rgb',  lineColor=None, fillColor='white',
    opacity=1.0, depth=-3.0, interpolate=True)
portImage = parallel.ParallelPort(address='0x0378')

# --- Initialize components for Routine "feedback" ---
textPoints = visual.TextStim(win=win, name='textPoints',
    text='',
    font='Open Sans',
    pos=(0, 0.3), height=0.06, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-1.0);
textScore = visual.TextStim(win=win, name='textScore',
    text='',
    font='Open Sans',
    pos=(0, 0.24), height=0.06, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-2.0);
polygon = visual.Rect(
    win=win, name='polygon',
    width=[1.0, 1.0][0], height=[1.0, 1.0][1],
    ori=0.0, pos=[0,0], anchor='center',
    lineWidth=100.0,     colorSpace='rgb',  lineColor='white', fillColor=[0.0000, 0.0000, 0.0000],
    opacity=1.0, depth=-3.0, interpolate=True)
imageFeedback = visual.ImageStim(
    win=win,
    name='imageFeedback', 
    image='sin', mask=None, anchor='center',
    ori=0.0, pos=[0,0], size=None,
    color=[1,1,1], colorSpace='rgb', opacity=None,
    flipHoriz=False, flipVert=False,
    texRes=128.0, interpolate=True, depth=-4.0)
gazeCursor_Feedback = visual.ShapeStim(
    win=win, name='gazeCursor_Feedback', vertices='star7',
    size=(0.02, 0.02),
    ori=0.0, pos=[0,0], anchor='center',
    lineWidth=1.0,     colorSpace='rgb',  lineColor=None, fillColor='white',
    opacity=1.0, depth=-5.0, interpolate=True)
soundFeedback = sound.Sound('A', secs=1, stereo=True, hamming=True,
    name='soundFeedback')
soundFeedback.setVolume(1.0)
portShock1 = parallel.ParallelPort(address='0x0378')
portShock2 = parallel.ParallelPort(address='0x0378')
portShock3 = parallel.ParallelPort(address='0x0378')

# --- Initialize components for Routine "crossEnd" ---
fixcrossE = visual.TextStim(win=win, name='fixcrossE',
    text='+',
    font='Open Sans',
    pos=(0, 0), height=0.1, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-1.0);

# --- Initialize components for Routine "painRating" ---
textRate = visual.TextStim(win=win, name='textRate',
    text='Wie schmerzhaft fanden Sie den elektrischen Reiz?',
    font='Open Sans',
    pos=(0, 0.6), height=0.06, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=0.0);
text_noPain = visual.TextStim(win=win, name='text_noPain',
    text='nichts\ngespürt',
    font='Open Sans',
    pos=(-0.5, 0.15), height=0.05, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-1.0);
text_slightPain = visual.TextStim(win=win, name='text_slightPain',
    text='eben \nwahrnehmbarer\nSchmerz',
    font='Open Sans',
    pos=(-0.1, 0.15), height=0.05, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-2.0);
text_highPain = visual.TextStim(win=win, name='text_highPain',
    text='unerträglicher\nSchmerz',
    font='Open Sans',
    pos=(0.5, 0.15), height=0.05, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-3.0);
sliderPain = visual.Slider(win=win, name='sliderPain',
    startValue=None, size=(1.0, 0.1), pos=(0, 0), units=None,
    labels=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10), ticks=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10), granularity=1.0,
    style='rating', styleTweaks=('labels45',), opacity=None,
    labelColor='White', markerColor='Red', lineColor='White', colorSpace='rgb',
    font='Open Sans', labelHeight=0.03,
    flip=False, ori=0.0, depth=-4, readOnly=False)
textSpacePain = visual.TextStim(win=win, name='textSpacePain',
    text='Bitte drücken Sie die Leertaste, um die Antwort zu bestätigen.',
    font='Open Sans',
    pos=(0, -0.3), height=0.06, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-5.0);
spacePain = keyboard.Keyboard()

# --- Initialize components for Routine "stimRating" ---
textRateStim = visual.TextStim(win=win, name='textRateStim',
    text='Wie fühlen Sie sich bei der Betrachtung des Bildes?',
    font='Open Sans',
    pos=(0, 0.6), height=0.06, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-1.0);
imageRating = visual.ImageStim(
    win=win,
    name='imageRating', 
    image='sin', mask=None, anchor='center',
    ori=0.0, pos=[0,0], size=None,
    color=[1,1,1], colorSpace='rgb', opacity=None,
    flipHoriz=False, flipVert=False,
    texRes=128.0, interpolate=True, depth=-2.0)
textUnpleasant = visual.TextStim(win=win, name='textUnpleasant',
    text='sehr \nunwohl',
    font='Open Sans',
    pos=(-0.5, -0.15), height=0.05, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-3.0);
textPleasant = visual.TextStim(win=win, name='textPleasant',
    text='sehr\nwohl',
    font='Open Sans',
    pos=(0.5, -0.15), height=0.05, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-4.0);
sliderStim = visual.Slider(win=win, name='sliderStim',
    startValue=None, size=(1.0, 0.1), pos=(0, -0.3), units=None,
    labels=(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), ticks=(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), granularity=1.0,
    style='rating', styleTweaks=('labels45',), opacity=None,
    labelColor='White', markerColor='Red', lineColor='White', colorSpace='rgb',
    font='Open Sans', labelHeight=0.03,
    flip=False, ori=0.0, depth=-5, readOnly=False)
textSpaceStim = visual.TextStim(win=win, name='textSpaceStim',
    text='Bitte drücken Sie die Leertaste, um die Antwort zu bestätigen.',
    font='Open Sans',
    pos=(0, -0.6), height=0.06, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-6.0);
spaceStim = keyboard.Keyboard()

# --- Initialize components for Routine "testInstr" ---
textTestInstr = visual.TextStim(win=win, name='textTestInstr',
    text='Im nächsten Block sehen Sie noch einmal die unterschiedlichen Bilder.\n\nBitte drücken Sie die Leertaste, um zu starten.',
    font='Open Sans',
    pos=(0, 0), height=0.06, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=0.0);
spaceTestInstr = keyboard.Keyboard()

# --- Initialize components for Routine "cross" ---
fixcross = visual.TextStim(win=win, name='fixcross',
    text='+',
    font='Open Sans',
    pos=(0, 0), height=0.1, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-1.0);

# --- Initialize components for Routine "testtrial" ---
imageTest = visual.ImageStim(
    win=win,
    name='imageTest', 
    image='sin', mask=None, anchor='center',
    ori=0.0, pos=[0,0], size=1.0,
    color=[1,1,1], colorSpace='rgb', opacity=None,
    flipHoriz=False, flipVert=False,
    texRes=128.0, interpolate=True, depth=-1.0)
roiTest = visual.ROI(win, name='roiTest', device=eyetracker,
    debug=False,
    shape='rectangle',
    pos=(0, 0), size=1.0, anchor='center', ori=0.0)
portTestImage = parallel.ParallelPort(address='0x0378')

# --- Initialize components for Routine "crossEnd" ---
fixcrossE = visual.TextStim(win=win, name='fixcrossE',
    text='+',
    font='Open Sans',
    pos=(0, 0), height=0.1, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-1.0);

# --- Initialize components for Routine "stimRating" ---
textRateStim = visual.TextStim(win=win, name='textRateStim',
    text='Wie fühlen Sie sich bei der Betrachtung des Bildes?',
    font='Open Sans',
    pos=(0, 0.6), height=0.06, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-1.0);
imageRating = visual.ImageStim(
    win=win,
    name='imageRating', 
    image='sin', mask=None, anchor='center',
    ori=0.0, pos=[0,0], size=None,
    color=[1,1,1], colorSpace='rgb', opacity=None,
    flipHoriz=False, flipVert=False,
    texRes=128.0, interpolate=True, depth=-2.0)
textUnpleasant = visual.TextStim(win=win, name='textUnpleasant',
    text='sehr \nunwohl',
    font='Open Sans',
    pos=(-0.5, -0.15), height=0.05, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-3.0);
textPleasant = visual.TextStim(win=win, name='textPleasant',
    text='sehr\nwohl',
    font='Open Sans',
    pos=(0.5, -0.15), height=0.05, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-4.0);
sliderStim = visual.Slider(win=win, name='sliderStim',
    startValue=None, size=(1.0, 0.1), pos=(0, -0.3), units=None,
    labels=(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), ticks=(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), granularity=1.0,
    style='rating', styleTweaks=('labels45',), opacity=None,
    labelColor='White', markerColor='Red', lineColor='White', colorSpace='rgb',
    font='Open Sans', labelHeight=0.03,
    flip=False, ori=0.0, depth=-5, readOnly=False)
textSpaceStim = visual.TextStim(win=win, name='textSpaceStim',
    text='Bitte drücken Sie die Leertaste, um die Antwort zu bestätigen.',
    font='Open Sans',
    pos=(0, -0.6), height=0.06, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-6.0);
spaceStim = keyboard.Keyboard()

# --- Initialize components for Routine "end" ---
textEnd = visual.TextStim(win=win, name='textEnd',
    text='Sehr geehrter Teilnehmer, sehr geehrte Teilnehmerin,\nvielen Dank für Ihre Teilnahme an unserem Experiment.\n\nDas Experiment ist beendet. \nSie können der Versuchsleitung Bescheid sagen.',
    font='Open Sans',
    pos=(0, 0), height=0.06, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=0.0);

# Create some handy timers
globalClock = core.Clock()  # to track the time since experiment started
routineTimer = core.Clock()  # to track time remaining of each (possibly non-slip) routine 

# --- Prepare to start Routine "welcome" ---
continueRoutine = True
routineForceEnded = False
# update component parameters for each repeat
spaceStart.keys = []
spaceStart.rt = []
_spaceStart_allKeys = []
# Run 'Begin Routine' code from codeWelcome
print("Welcome Screen. Press Space to continue.")
# keep track of which components have finished
welcomeComponents = [textStart, spaceStart]
for thisComponent in welcomeComponents:
    thisComponent.tStart = None
    thisComponent.tStop = None
    thisComponent.tStartRefresh = None
    thisComponent.tStopRefresh = None
    if hasattr(thisComponent, 'status'):
        thisComponent.status = NOT_STARTED
# reset timers
t = 0
_timeToFirstFrame = win.getFutureFlipTime(clock="now")
frameN = -1

# --- Run Routine "welcome" ---
while continueRoutine:
    # get current time
    t = routineTimer.getTime()
    tThisFlip = win.getFutureFlipTime(clock=routineTimer)
    tThisFlipGlobal = win.getFutureFlipTime(clock=None)
    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
    # update/draw components on each frame
    
    # *textStart* updates
    if textStart.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        # keep track of start time/frame for later
        textStart.frameNStart = frameN  # exact frame index
        textStart.tStart = t  # local t and not account for scr refresh
        textStart.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(textStart, 'tStartRefresh')  # time at next scr refresh
        # add timestamp to datafile
        thisExp.timestampOnFlip(win, 'textStart.started')
        textStart.setAutoDraw(True)
    
    # *spaceStart* updates
    waitOnFlip = False
    if spaceStart.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        # keep track of start time/frame for later
        spaceStart.frameNStart = frameN  # exact frame index
        spaceStart.tStart = t  # local t and not account for scr refresh
        spaceStart.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(spaceStart, 'tStartRefresh')  # time at next scr refresh
        # add timestamp to datafile
        thisExp.timestampOnFlip(win, 'spaceStart.started')
        spaceStart.status = STARTED
        # keyboard checking is just starting
        waitOnFlip = True
        win.callOnFlip(spaceStart.clock.reset)  # t=0 on next screen flip
        win.callOnFlip(spaceStart.clearEvents, eventType='keyboard')  # clear events on next screen flip
    if spaceStart.status == STARTED and not waitOnFlip:
        theseKeys = spaceStart.getKeys(keyList=['space', 'enter'], waitRelease=False)
        _spaceStart_allKeys.extend(theseKeys)
        if len(_spaceStart_allKeys):
            spaceStart.keys = _spaceStart_allKeys[-1].name  # just the last key pressed
            spaceStart.rt = _spaceStart_allKeys[-1].rt
            # a response ends the routine
            continueRoutine = False
    
    # check for quit (typically the Esc key)
    if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
        core.quit()
    
    # check if all components have finished
    if not continueRoutine:  # a component has requested a forced-end of Routine
        routineForceEnded = True
        break
    continueRoutine = False  # will revert to True if at least one component still running
    for thisComponent in welcomeComponents:
        if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
            continueRoutine = True
            break  # at least one component has not yet finished
    
    # refresh the screen
    if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
        win.flip()

# --- Ending Routine "welcome" ---
for thisComponent in welcomeComponents:
    if hasattr(thisComponent, "setAutoDraw"):
        thisComponent.setAutoDraw(False)
# check responses
if spaceStart.keys in ['', [], None]:  # No response was made
    spaceStart.keys = None
thisExp.addData('spaceStart.keys',spaceStart.keys)
if spaceStart.keys != None:  # we had a response
    thisExp.addData('spaceStart.rt', spaceStart.rt)
thisExp.nextEntry()
# the Routine "welcome" was not non-slip safe, so reset the non-slip timer
routineTimer.reset()

# --- Prepare to start Routine "startPainThreshold" ---
continueRoutine = True
routineForceEnded = False
# update component parameters for each repeat
spacePainThreshold.keys = []
spacePainThreshold.rt = []
_spacePainThreshold_allKeys = []
# Run 'Begin Routine' code from codeStartPainThreshold
print("Pain Calibration Screen. Press Space to continue.")
# keep track of which components have finished
startPainThresholdComponents = [textPainThreshold, spacePainThreshold]
for thisComponent in startPainThresholdComponents:
    thisComponent.tStart = None
    thisComponent.tStop = None
    thisComponent.tStartRefresh = None
    thisComponent.tStopRefresh = None
    if hasattr(thisComponent, 'status'):
        thisComponent.status = NOT_STARTED
# reset timers
t = 0
_timeToFirstFrame = win.getFutureFlipTime(clock="now")
frameN = -1

# --- Run Routine "startPainThreshold" ---
while continueRoutine:
    # get current time
    t = routineTimer.getTime()
    tThisFlip = win.getFutureFlipTime(clock=routineTimer)
    tThisFlipGlobal = win.getFutureFlipTime(clock=None)
    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
    # update/draw components on each frame
    
    # *textPainThreshold* updates
    if textPainThreshold.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        # keep track of start time/frame for later
        textPainThreshold.frameNStart = frameN  # exact frame index
        textPainThreshold.tStart = t  # local t and not account for scr refresh
        textPainThreshold.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(textPainThreshold, 'tStartRefresh')  # time at next scr refresh
        # add timestamp to datafile
        thisExp.timestampOnFlip(win, 'textPainThreshold.started')
        textPainThreshold.setAutoDraw(True)
    
    # *spacePainThreshold* updates
    waitOnFlip = False
    if spacePainThreshold.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        # keep track of start time/frame for later
        spacePainThreshold.frameNStart = frameN  # exact frame index
        spacePainThreshold.tStart = t  # local t and not account for scr refresh
        spacePainThreshold.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(spacePainThreshold, 'tStartRefresh')  # time at next scr refresh
        # add timestamp to datafile
        thisExp.timestampOnFlip(win, 'spacePainThreshold.started')
        spacePainThreshold.status = STARTED
        # keyboard checking is just starting
        waitOnFlip = True
        win.callOnFlip(spacePainThreshold.clock.reset)  # t=0 on next screen flip
        win.callOnFlip(spacePainThreshold.clearEvents, eventType='keyboard')  # clear events on next screen flip
    if spacePainThreshold.status == STARTED and not waitOnFlip:
        theseKeys = spacePainThreshold.getKeys(keyList=['space', 'enter'], waitRelease=False)
        _spacePainThreshold_allKeys.extend(theseKeys)
        if len(_spacePainThreshold_allKeys):
            spacePainThreshold.keys = _spacePainThreshold_allKeys[-1].name  # just the last key pressed
            spacePainThreshold.rt = _spacePainThreshold_allKeys[-1].rt
            # a response ends the routine
            continueRoutine = False
    
    # check for quit (typically the Esc key)
    if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
        core.quit()
    
    # check if all components have finished
    if not continueRoutine:  # a component has requested a forced-end of Routine
        routineForceEnded = True
        break
    continueRoutine = False  # will revert to True if at least one component still running
    for thisComponent in startPainThresholdComponents:
        if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
            continueRoutine = True
            break  # at least one component has not yet finished
    
    # refresh the screen
    if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
        win.flip()

# --- Ending Routine "startPainThreshold" ---
for thisComponent in startPainThresholdComponents:
    if hasattr(thisComponent, "setAutoDraw"):
        thisComponent.setAutoDraw(False)
# check responses
if spacePainThreshold.keys in ['', [], None]:  # No response was made
    spacePainThreshold.keys = None
thisExp.addData('spacePainThreshold.keys',spacePainThreshold.keys)
if spacePainThreshold.keys != None:  # we had a response
    thisExp.addData('spacePainThreshold.rt', spacePainThreshold.rt)
thisExp.nextEntry()
# the Routine "startPainThreshold" was not non-slip safe, so reset the non-slip timer
routineTimer.reset()

# set up handler to look after randomisation of conditions etc
trials_threshold = data.TrialHandler(nReps=99999.0, method='random', 
    extraInfo=expInfo, originPath=-1,
    trialList=[None],
    seed=None, name='trials_threshold')
thisExp.addLoop(trials_threshold)  # add the loop to the experiment
thisTrials_threshold = trials_threshold.trialList[0]  # so we can initialise stimuli with some values
# abbreviate parameter names if possible (e.g. rgb = thisTrials_threshold.rgb)
if thisTrials_threshold != None:
    for paramName in thisTrials_threshold:
        exec('{} = thisTrials_threshold[paramName]'.format(paramName))

for thisTrials_threshold in trials_threshold:
    currentLoop = trials_threshold
    # abbreviate parameter names if possible (e.g. rgb = thisTrials_threshold.rgb)
    if thisTrials_threshold != None:
        for paramName in thisTrials_threshold:
            exec('{} = thisTrials_threshold[paramName]'.format(paramName))
    
    # --- Prepare to start Routine "painThreshold" ---
    continueRoutine = True
    routineForceEnded = False
    # update component parameters for each repeat
    sliderPainT.reset()
    keys_PainT.keys = []
    keys_PainT.rt = []
    _keys_PainT_allKeys = []
    # Run 'Begin Routine' code from codePainRatingT
    logging.log(level=logging.INFO, msg=f'PainThreshold')
    
    if trials_threshold.thisN == 0:
        print("Press 's' to administer shock.")
    else:
        print(f"Press 's' to administer shock. Press 'g' to continue when pain threshold was found.")
    
    values = [None]
    # keep track of which components have finished
    painThresholdComponents = [textRateT, text_noPainT, text_slightPainT, text_highPainT, sliderPainT, keys_PainT]
    for thisComponent in painThresholdComponents:
        thisComponent.tStart = None
        thisComponent.tStop = None
        thisComponent.tStartRefresh = None
        thisComponent.tStopRefresh = None
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    # reset timers
    t = 0
    _timeToFirstFrame = win.getFutureFlipTime(clock="now")
    frameN = -1
    
    # --- Run Routine "painThreshold" ---
    while continueRoutine:
        # get current time
        t = routineTimer.getTime()
        tThisFlip = win.getFutureFlipTime(clock=routineTimer)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        
        # *textRateT* updates
        if textRateT.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            textRateT.frameNStart = frameN  # exact frame index
            textRateT.tStart = t  # local t and not account for scr refresh
            textRateT.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(textRateT, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'textRateT.started')
            textRateT.setAutoDraw(True)
        
        # *text_noPainT* updates
        if text_noPainT.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            text_noPainT.frameNStart = frameN  # exact frame index
            text_noPainT.tStart = t  # local t and not account for scr refresh
            text_noPainT.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(text_noPainT, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'text_noPainT.started')
            text_noPainT.setAutoDraw(True)
        
        # *text_slightPainT* updates
        if text_slightPainT.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            text_slightPainT.frameNStart = frameN  # exact frame index
            text_slightPainT.tStart = t  # local t and not account for scr refresh
            text_slightPainT.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(text_slightPainT, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'text_slightPainT.started')
            text_slightPainT.setAutoDraw(True)
        
        # *text_highPainT* updates
        if text_highPainT.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            text_highPainT.frameNStart = frameN  # exact frame index
            text_highPainT.tStart = t  # local t and not account for scr refresh
            text_highPainT.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(text_highPainT, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'text_highPainT.started')
            text_highPainT.setAutoDraw(True)
        
        # *sliderPainT* updates
        if sliderPainT.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            sliderPainT.frameNStart = frameN  # exact frame index
            sliderPainT.tStart = t  # local t and not account for scr refresh
            sliderPainT.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(sliderPainT, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'sliderPainT.started')
            sliderPainT.setAutoDraw(True)
        
        # *keys_PainT* updates
        waitOnFlip = False
        if keys_PainT.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
            # keep track of start time/frame for later
            keys_PainT.frameNStart = frameN  # exact frame index
            keys_PainT.tStart = t  # local t and not account for scr refresh
            keys_PainT.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(keys_PainT, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'keys_PainT.started')
            keys_PainT.status = STARTED
            # keyboard checking is just starting
            waitOnFlip = True
            win.callOnFlip(keys_PainT.clock.reset)  # t=0 on next screen flip
            win.callOnFlip(keys_PainT.clearEvents, eventType='keyboard')  # clear events on next screen flip
        if keys_PainT.status == STARTED and not waitOnFlip:
            theseKeys = keys_PainT.getKeys(keyList=['s', 'g'], waitRelease=False)
            _keys_PainT_allKeys.extend(theseKeys)
            if len(_keys_PainT_allKeys):
                keys_PainT.keys = _keys_PainT_allKeys[-1].name  # just the last key pressed
                keys_PainT.rt = _keys_PainT_allKeys[-1].rt
                # a response ends the routine
                continueRoutine = False
        # Run 'Each Frame' code from codePainRatingT
        if 'g' in keys_PainT.keys:
            continueRoutine = False
            trials_threshold.finished = True
            skipShock = True
            
        if sliderPainT.getRating() != values[-1]:
            values.append(sliderPainT.getRating())
            print(f"Value: {values[-1]}")
        
        # check for quit (typically the Esc key)
        if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
            core.quit()
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            routineForceEnded = True
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in painThresholdComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # --- Ending Routine "painThreshold" ---
    for thisComponent in painThresholdComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    trials_threshold.addData('sliderPainT.response', sliderPainT.getRating())
    trials_threshold.addData('sliderPainT.rt', sliderPainT.getRT())
    # check responses
    if keys_PainT.keys in ['', [], None]:  # No response was made
        keys_PainT.keys = None
    trials_threshold.addData('keys_PainT.keys',keys_PainT.keys)
    if keys_PainT.keys != None:  # we had a response
        trials_threshold.addData('keys_PainT.rt', keys_PainT.rt)
    # the Routine "painThreshold" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset()
    
    # --- Prepare to start Routine "shock" ---
    continueRoutine = True
    routineForceEnded = False
    # update component parameters for each repeat
    # Run 'Begin Routine' code from code
        
    # keep track of which components have finished
    shockComponents = [portShock1T, portShock2T, portShock3T, fixcrossT]
    for thisComponent in shockComponents:
        thisComponent.tStart = None
        thisComponent.tStop = None
        thisComponent.tStartRefresh = None
        thisComponent.tStopRefresh = None
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    # reset timers
    t = 0
    _timeToFirstFrame = win.getFutureFlipTime(clock="now")
    frameN = -1
    
    # --- Run Routine "shock" ---
    while continueRoutine and routineTimer.getTime() < 0.5:
        # get current time
        t = routineTimer.getTime()
        tThisFlip = win.getFutureFlipTime(clock=routineTimer)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        # Run 'Each Frame' code from code
        if skipShock:
            continueRoutine = False
        # *portShock1T* updates
        if portShock1T.status == NOT_STARTED and t >= 0.1-frameTolerance:
            # keep track of start time/frame for later
            portShock1T.frameNStart = frameN  # exact frame index
            portShock1T.tStart = t  # local t and not account for scr refresh
            portShock1T.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(portShock1T, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.addData('portShock1T.started', t)
            portShock1T.status = STARTED
            win.callOnFlip(portShock1T.setData, int(132))
        if portShock1T.status == STARTED:
            # is it time to stop? (based on global clock, using actual start)
            if tThisFlipGlobal > portShock1T.tStartRefresh + 0.005-frameTolerance:
                # keep track of stop time/frame for later
                portShock1T.tStop = t  # not accounting for scr refresh
                portShock1T.frameNStop = frameN  # exact frame index
                # add timestamp to datafile
                thisExp.addData('portShock1T.stopped', t)
                portShock1T.status = FINISHED
                win.callOnFlip(portShock1T.setData, int(0))
        # *portShock2T* updates
        if portShock2T.status == NOT_STARTED and t >= 0.15-frameTolerance:
            # keep track of start time/frame for later
            portShock2T.frameNStart = frameN  # exact frame index
            portShock2T.tStart = t  # local t and not account for scr refresh
            portShock2T.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(portShock2T, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.addData('portShock2T.started', t)
            portShock2T.status = STARTED
            win.callOnFlip(portShock2T.setData, int(132))
        if portShock2T.status == STARTED:
            # is it time to stop? (based on global clock, using actual start)
            if tThisFlipGlobal > portShock2T.tStartRefresh + 0.005-frameTolerance:
                # keep track of stop time/frame for later
                portShock2T.tStop = t  # not accounting for scr refresh
                portShock2T.frameNStop = frameN  # exact frame index
                # add timestamp to datafile
                thisExp.addData('portShock2T.stopped', t)
                portShock2T.status = FINISHED
                win.callOnFlip(portShock2T.setData, int(0))
        # *portShock3T* updates
        if portShock3T.status == NOT_STARTED and t >= 0.2-frameTolerance:
            # keep track of start time/frame for later
            portShock3T.frameNStart = frameN  # exact frame index
            portShock3T.tStart = t  # local t and not account for scr refresh
            portShock3T.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(portShock3T, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.addData('portShock3T.started', t)
            portShock3T.status = STARTED
            win.callOnFlip(portShock3T.setData, int(132))
        if portShock3T.status == STARTED:
            # is it time to stop? (based on global clock, using actual start)
            if tThisFlipGlobal > portShock3T.tStartRefresh + 0.005-frameTolerance:
                # keep track of stop time/frame for later
                portShock3T.tStop = t  # not accounting for scr refresh
                portShock3T.frameNStop = frameN  # exact frame index
                # add timestamp to datafile
                thisExp.addData('portShock3T.stopped', t)
                portShock3T.status = FINISHED
                win.callOnFlip(portShock3T.setData, int(0))
        
        # *fixcrossT* updates
        if fixcrossT.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            fixcrossT.frameNStart = frameN  # exact frame index
            fixcrossT.tStart = t  # local t and not account for scr refresh
            fixcrossT.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(fixcrossT, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'fixcrossT.started')
            fixcrossT.setAutoDraw(True)
        if fixcrossT.status == STARTED:
            # is it time to stop? (based on global clock, using actual start)
            if tThisFlipGlobal > fixcrossT.tStartRefresh + 0.5-frameTolerance:
                # keep track of stop time/frame for later
                fixcrossT.tStop = t  # not accounting for scr refresh
                fixcrossT.frameNStop = frameN  # exact frame index
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'fixcrossT.stopped')
                fixcrossT.setAutoDraw(False)
        
        # check for quit (typically the Esc key)
        if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
            core.quit()
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            routineForceEnded = True
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in shockComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # --- Ending Routine "shock" ---
    for thisComponent in shockComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    if portShock1T.status == STARTED:
        win.callOnFlip(portShock1T.setData, int(0))
    if portShock2T.status == STARTED:
        win.callOnFlip(portShock2T.setData, int(0))
    if portShock3T.status == STARTED:
        win.callOnFlip(portShock3T.setData, int(0))
    # using non-slip timing so subtract the expected duration of this Routine (unless ended on request)
    if routineForceEnded:
        routineTimer.reset()
    else:
        routineTimer.addTime(-0.500000)
    thisExp.nextEntry()
    
# completed 99999.0 repeats of 'trials_threshold'


# --- Prepare to start Routine "startETCalibration" ---
continueRoutine = True
routineForceEnded = False
# update component parameters for each repeat
spaceETCalibration.keys = []
spaceETCalibration.rt = []
_spaceETCalibration_allKeys = []
# Run 'Begin Routine' code from codeStartETCalibration
print("ET Calibration Screen. Press Space to continue.")
# keep track of which components have finished
startETCalibrationComponents = [text_ETCalibration, spaceETCalibration]
for thisComponent in startETCalibrationComponents:
    thisComponent.tStart = None
    thisComponent.tStop = None
    thisComponent.tStartRefresh = None
    thisComponent.tStopRefresh = None
    if hasattr(thisComponent, 'status'):
        thisComponent.status = NOT_STARTED
# reset timers
t = 0
_timeToFirstFrame = win.getFutureFlipTime(clock="now")
frameN = -1

# --- Run Routine "startETCalibration" ---
while continueRoutine:
    # get current time
    t = routineTimer.getTime()
    tThisFlip = win.getFutureFlipTime(clock=routineTimer)
    tThisFlipGlobal = win.getFutureFlipTime(clock=None)
    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
    # update/draw components on each frame
    
    # *text_ETCalibration* updates
    if text_ETCalibration.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        # keep track of start time/frame for later
        text_ETCalibration.frameNStart = frameN  # exact frame index
        text_ETCalibration.tStart = t  # local t and not account for scr refresh
        text_ETCalibration.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(text_ETCalibration, 'tStartRefresh')  # time at next scr refresh
        # add timestamp to datafile
        thisExp.timestampOnFlip(win, 'text_ETCalibration.started')
        text_ETCalibration.setAutoDraw(True)
    
    # *spaceETCalibration* updates
    waitOnFlip = False
    if spaceETCalibration.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        # keep track of start time/frame for later
        spaceETCalibration.frameNStart = frameN  # exact frame index
        spaceETCalibration.tStart = t  # local t and not account for scr refresh
        spaceETCalibration.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(spaceETCalibration, 'tStartRefresh')  # time at next scr refresh
        # add timestamp to datafile
        thisExp.timestampOnFlip(win, 'spaceETCalibration.started')
        spaceETCalibration.status = STARTED
        # keyboard checking is just starting
        waitOnFlip = True
        win.callOnFlip(spaceETCalibration.clock.reset)  # t=0 on next screen flip
        win.callOnFlip(spaceETCalibration.clearEvents, eventType='keyboard')  # clear events on next screen flip
    if spaceETCalibration.status == STARTED and not waitOnFlip:
        theseKeys = spaceETCalibration.getKeys(keyList=['space', 'enter'], waitRelease=False)
        _spaceETCalibration_allKeys.extend(theseKeys)
        if len(_spaceETCalibration_allKeys):
            spaceETCalibration.keys = _spaceETCalibration_allKeys[-1].name  # just the last key pressed
            spaceETCalibration.rt = _spaceETCalibration_allKeys[-1].rt
            # a response ends the routine
            continueRoutine = False
    
    # check for quit (typically the Esc key)
    if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
        core.quit()
    
    # check if all components have finished
    if not continueRoutine:  # a component has requested a forced-end of Routine
        routineForceEnded = True
        break
    continueRoutine = False  # will revert to True if at least one component still running
    for thisComponent in startETCalibrationComponents:
        if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
            continueRoutine = True
            break  # at least one component has not yet finished
    
    # refresh the screen
    if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
        win.flip()

# --- Ending Routine "startETCalibration" ---
for thisComponent in startETCalibrationComponents:
    if hasattr(thisComponent, "setAutoDraw"):
        thisComponent.setAutoDraw(False)
# check responses
if spaceETCalibration.keys in ['', [], None]:  # No response was made
    spaceETCalibration.keys = None
thisExp.addData('spaceETCalibration.keys',spaceETCalibration.keys)
if spaceETCalibration.keys != None:  # we had a response
    thisExp.addData('spaceETCalibration.rt', spaceETCalibration.rt)
thisExp.nextEntry()
# the Routine "startETCalibration" was not non-slip safe, so reset the non-slip timer
routineTimer.reset()
# define target for ETCalibration
ETCalibrationTarget = visual.TargetStim(win, 
    name='ETCalibrationTarget',
    radius=0.01, fillColor='', borderColor='black', lineWidth=2.0,
    innerRadius=0.0035, innerFillColor='white', innerBorderColor='black', innerLineWidth=2.0,
    colorSpace='rgb', units=None
)
# define parameters for ETCalibration
ETCalibration = hardware.eyetracker.EyetrackerCalibration(win, 
    eyetracker, ETCalibrationTarget,
    units=None, colorSpace='rgb',
    progressMode='time', targetDur=1.0, expandScale=1.2,
    targetLayout='NINE_POINTS', randomisePos=True, textColor='white',
    movementAnimation=True, targetDelay=1.0
)
# run calibration
ETCalibration.run()
# clear any keypresses from during ETCalibration so they don't interfere with the experiment
defaultKeyboard.clearEvents()
# the Routine "ETCalibration" was not non-slip safe, so reset the non-slip timer
routineTimer.reset()

# set up handler to look after randomisation of conditions etc
blocks = data.TrialHandler(nReps=1.0, method='random', 
    extraInfo=expInfo, originPath=-1,
    trialList=data.importConditions("chooseBlock"+expInfo['group']+".xlsx"),
    seed=None, name='blocks')
thisExp.addLoop(blocks)  # add the loop to the experiment
thisBlock = blocks.trialList[0]  # so we can initialise stimuli with some values
# abbreviate parameter names if possible (e.g. rgb = thisBlock.rgb)
if thisBlock != None:
    for paramName in thisBlock:
        exec('{} = thisBlock[paramName]'.format(paramName))

for thisBlock in blocks:
    currentLoop = blocks
    # abbreviate parameter names if possible (e.g. rgb = thisBlock.rgb)
    if thisBlock != None:
        for paramName in thisBlock:
            exec('{} = thisBlock[paramName]'.format(paramName))
    
    # --- Prepare to start Routine "startExp" ---
    continueRoutine = True
    routineForceEnded = False
    # update component parameters for each repeat
    spaceStartExp.keys = []
    spaceStartExp.rt = []
    _spaceStartExp_allKeys = []
    # Run 'Begin Routine' code from codeStartExp
    print("Start Experiment Screen. Press Space to continue.")
    # keep track of which components have finished
    startExpComponents = [textStartExp, spaceStartExp]
    for thisComponent in startExpComponents:
        thisComponent.tStart = None
        thisComponent.tStop = None
        thisComponent.tStartRefresh = None
        thisComponent.tStopRefresh = None
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    # reset timers
    t = 0
    _timeToFirstFrame = win.getFutureFlipTime(clock="now")
    frameN = -1
    
    # --- Run Routine "startExp" ---
    while continueRoutine:
        # get current time
        t = routineTimer.getTime()
        tThisFlip = win.getFutureFlipTime(clock=routineTimer)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        
        # *textStartExp* updates
        if textStartExp.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            textStartExp.frameNStart = frameN  # exact frame index
            textStartExp.tStart = t  # local t and not account for scr refresh
            textStartExp.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(textStartExp, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'textStartExp.started')
            textStartExp.setAutoDraw(True)
        
        # *spaceStartExp* updates
        waitOnFlip = False
        if spaceStartExp.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            spaceStartExp.frameNStart = frameN  # exact frame index
            spaceStartExp.tStart = t  # local t and not account for scr refresh
            spaceStartExp.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(spaceStartExp, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'spaceStartExp.started')
            spaceStartExp.status = STARTED
            # keyboard checking is just starting
            waitOnFlip = True
            win.callOnFlip(spaceStartExp.clock.reset)  # t=0 on next screen flip
            win.callOnFlip(spaceStartExp.clearEvents, eventType='keyboard')  # clear events on next screen flip
        if spaceStartExp.status == STARTED and not waitOnFlip:
            theseKeys = spaceStartExp.getKeys(keyList=['space', 'enter'], waitRelease=False)
            _spaceStartExp_allKeys.extend(theseKeys)
            if len(_spaceStartExp_allKeys):
                spaceStartExp.keys = _spaceStartExp_allKeys[-1].name  # just the last key pressed
                spaceStartExp.rt = _spaceStartExp_allKeys[-1].rt
                # a response ends the routine
                continueRoutine = False
        
        # check for quit (typically the Esc key)
        if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
            core.quit()
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            routineForceEnded = True
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in startExpComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # --- Ending Routine "startExp" ---
    for thisComponent in startExpComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    # check responses
    if spaceStartExp.keys in ['', [], None]:  # No response was made
        spaceStartExp.keys = None
    blocks.addData('spaceStartExp.keys',spaceStartExp.keys)
    if spaceStartExp.keys != None:  # we had a response
        blocks.addData('spaceStartExp.rt', spaceStartExp.rt)
    # the Routine "startExp" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset()
    
    # set up handler to look after randomisation of conditions etc
    ratingtrials1 = data.TrialHandler(nReps=1.0, method='random', 
        extraInfo=expInfo, originPath=-1,
        trialList=data.importConditions(stimFile),
        seed=None, name='ratingtrials1')
    thisExp.addLoop(ratingtrials1)  # add the loop to the experiment
    thisRatingtrials1 = ratingtrials1.trialList[0]  # so we can initialise stimuli with some values
    # abbreviate parameter names if possible (e.g. rgb = thisRatingtrials1.rgb)
    if thisRatingtrials1 != None:
        for paramName in thisRatingtrials1:
            exec('{} = thisRatingtrials1[paramName]'.format(paramName))
    
    for thisRatingtrials1 in ratingtrials1:
        currentLoop = ratingtrials1
        # abbreviate parameter names if possible (e.g. rgb = thisRatingtrials1.rgb)
        if thisRatingtrials1 != None:
            for paramName in thisRatingtrials1:
                exec('{} = thisRatingtrials1[paramName]'.format(paramName))
        
        # --- Prepare to start Routine "stimRating" ---
        continueRoutine = True
        routineForceEnded = False
        # update component parameters for each repeat
        # Run 'Begin Routine' code from codeStimRating
        logging.log(level=logging.INFO, msg=f'StimRating_{stimtype}')
        print("Rating of Stimulus: %s"%(stimtype))
        win.mouseVisible = True
        imageRating.setPos((0, 0.2))
        imageRating.setImage(eval(stimtype))
        sliderStim.reset()
        spaceStim.keys = []
        spaceStim.rt = []
        _spaceStim_allKeys = []
        # keep track of which components have finished
        stimRatingComponents = [textRateStim, imageRating, textUnpleasant, textPleasant, sliderStim, textSpaceStim, spaceStim]
        for thisComponent in stimRatingComponents:
            thisComponent.tStart = None
            thisComponent.tStop = None
            thisComponent.tStartRefresh = None
            thisComponent.tStopRefresh = None
            if hasattr(thisComponent, 'status'):
                thisComponent.status = NOT_STARTED
        # reset timers
        t = 0
        _timeToFirstFrame = win.getFutureFlipTime(clock="now")
        frameN = -1
        
        # --- Run Routine "stimRating" ---
        while continueRoutine:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # *textRateStim* updates
            if textRateStim.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                textRateStim.frameNStart = frameN  # exact frame index
                textRateStim.tStart = t  # local t and not account for scr refresh
                textRateStim.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(textRateStim, 'tStartRefresh')  # time at next scr refresh
                textRateStim.setAutoDraw(True)
            
            # *imageRating* updates
            if imageRating.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                # keep track of start time/frame for later
                imageRating.frameNStart = frameN  # exact frame index
                imageRating.tStart = t  # local t and not account for scr refresh
                imageRating.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(imageRating, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'imageRating.started')
                imageRating.setAutoDraw(True)
            
            # *textUnpleasant* updates
            if textUnpleasant.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                textUnpleasant.frameNStart = frameN  # exact frame index
                textUnpleasant.tStart = t  # local t and not account for scr refresh
                textUnpleasant.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(textUnpleasant, 'tStartRefresh')  # time at next scr refresh
                textUnpleasant.setAutoDraw(True)
            
            # *textPleasant* updates
            if textPleasant.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                textPleasant.frameNStart = frameN  # exact frame index
                textPleasant.tStart = t  # local t and not account for scr refresh
                textPleasant.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(textPleasant, 'tStartRefresh')  # time at next scr refresh
                textPleasant.setAutoDraw(True)
            
            # *sliderStim* updates
            if sliderStim.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                sliderStim.frameNStart = frameN  # exact frame index
                sliderStim.tStart = t  # local t and not account for scr refresh
                sliderStim.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sliderStim, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'sliderStim.started')
                sliderStim.setAutoDraw(True)
            
            # *textSpaceStim* updates
            if textSpaceStim.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                textSpaceStim.frameNStart = frameN  # exact frame index
                textSpaceStim.tStart = t  # local t and not account for scr refresh
                textSpaceStim.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(textSpaceStim, 'tStartRefresh')  # time at next scr refresh
                textSpaceStim.setAutoDraw(True)
            
            # *spaceStim* updates
            waitOnFlip = False
            if spaceStim.status == NOT_STARTED and sliderStim.rating:
                # keep track of start time/frame for later
                spaceStim.frameNStart = frameN  # exact frame index
                spaceStim.tStart = t  # local t and not account for scr refresh
                spaceStim.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(spaceStim, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'spaceStim.started')
                spaceStim.status = STARTED
                # keyboard checking is just starting
                waitOnFlip = True
                win.callOnFlip(spaceStim.clock.reset)  # t=0 on next screen flip
                win.callOnFlip(spaceStim.clearEvents, eventType='keyboard')  # clear events on next screen flip
            if spaceStim.status == STARTED and not waitOnFlip:
                theseKeys = spaceStim.getKeys(keyList=['space', 'enter'], waitRelease=False)
                _spaceStim_allKeys.extend(theseKeys)
                if len(_spaceStim_allKeys):
                    spaceStim.keys = _spaceStim_allKeys[-1].name  # just the last key pressed
                    spaceStim.rt = _spaceStim_allKeys[-1].rt
                    # a response ends the routine
                    continueRoutine = False
            
            # check for quit (typically the Esc key)
            if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
                core.quit()
            
            # check if all components have finished
            if not continueRoutine:  # a component has requested a forced-end of Routine
                routineForceEnded = True
                break
            continueRoutine = False  # will revert to True if at least one component still running
            for thisComponent in stimRatingComponents:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # --- Ending Routine "stimRating" ---
        for thisComponent in stimRatingComponents:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        ratingtrials1.addData('sliderStim.response', sliderStim.getRating())
        ratingtrials1.addData('sliderStim.rt', sliderStim.getRT())
        # check responses
        if spaceStim.keys in ['', [], None]:  # No response was made
            spaceStim.keys = None
        ratingtrials1.addData('spaceStim.keys',spaceStim.keys)
        if spaceStim.keys != None:  # we had a response
            ratingtrials1.addData('spaceStim.rt', spaceStim.rt)
        # the Routine "stimRating" was not non-slip safe, so reset the non-slip timer
        routineTimer.reset()
        thisExp.nextEntry()
        
    # completed 1.0 repeats of 'ratingtrials1'
    
    
    # --- Prepare to start Routine "startTask" ---
    continueRoutine = True
    routineForceEnded = False
    # update component parameters for each repeat
    spaceStartTask.keys = []
    spaceStartTask.rt = []
    _spaceStartTask_allKeys = []
    # Run 'Begin Routine' code from codeStartTask
    print("Start Task Screen. Press Space to continue.")
    win.mouseVisible = False
    # keep track of which components have finished
    startTaskComponents = [textStartTask, spaceStartTask]
    for thisComponent in startTaskComponents:
        thisComponent.tStart = None
        thisComponent.tStop = None
        thisComponent.tStartRefresh = None
        thisComponent.tStopRefresh = None
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    # reset timers
    t = 0
    _timeToFirstFrame = win.getFutureFlipTime(clock="now")
    frameN = -1
    
    # --- Run Routine "startTask" ---
    while continueRoutine:
        # get current time
        t = routineTimer.getTime()
        tThisFlip = win.getFutureFlipTime(clock=routineTimer)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        
        # *textStartTask* updates
        if textStartTask.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            textStartTask.frameNStart = frameN  # exact frame index
            textStartTask.tStart = t  # local t and not account for scr refresh
            textStartTask.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(textStartTask, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'textStartTask.started')
            textStartTask.setAutoDraw(True)
        
        # *spaceStartTask* updates
        waitOnFlip = False
        if spaceStartTask.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            spaceStartTask.frameNStart = frameN  # exact frame index
            spaceStartTask.tStart = t  # local t and not account for scr refresh
            spaceStartTask.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(spaceStartTask, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'spaceStartTask.started')
            spaceStartTask.status = STARTED
            # keyboard checking is just starting
            waitOnFlip = True
            win.callOnFlip(spaceStartTask.clock.reset)  # t=0 on next screen flip
            win.callOnFlip(spaceStartTask.clearEvents, eventType='keyboard')  # clear events on next screen flip
        if spaceStartTask.status == STARTED and not waitOnFlip:
            theseKeys = spaceStartTask.getKeys(keyList=['space', 'enter'], waitRelease=False)
            _spaceStartTask_allKeys.extend(theseKeys)
            if len(_spaceStartTask_allKeys):
                spaceStartTask.keys = _spaceStartTask_allKeys[-1].name  # just the last key pressed
                spaceStartTask.rt = _spaceStartTask_allKeys[-1].rt
                # a response ends the routine
                continueRoutine = False
        
        # check for quit (typically the Esc key)
        if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
            core.quit()
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            routineForceEnded = True
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in startTaskComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # --- Ending Routine "startTask" ---
    for thisComponent in startTaskComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    # check responses
    if spaceStartTask.keys in ['', [], None]:  # No response was made
        spaceStartTask.keys = None
    blocks.addData('spaceStartTask.keys',spaceStartTask.keys)
    if spaceStartTask.keys != None:  # we had a response
        blocks.addData('spaceStartTask.rt', spaceStartTask.rt)
    # the Routine "startTask" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset()
    
    # --- Prepare to start Routine "blank" ---
    continueRoutine = True
    routineForceEnded = False
    # update component parameters for each repeat
    # keep track of which components have finished
    blankComponents = [blankScreen]
    for thisComponent in blankComponents:
        thisComponent.tStart = None
        thisComponent.tStop = None
        thisComponent.tStartRefresh = None
        thisComponent.tStopRefresh = None
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    # reset timers
    t = 0
    _timeToFirstFrame = win.getFutureFlipTime(clock="now")
    frameN = -1
    
    # --- Run Routine "blank" ---
    while continueRoutine:
        # get current time
        t = routineTimer.getTime()
        tThisFlip = win.getFutureFlipTime(clock=routineTimer)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        
        # *blankScreen* updates
        if blankScreen.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            blankScreen.frameNStart = frameN  # exact frame index
            blankScreen.tStart = t  # local t and not account for scr refresh
            blankScreen.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(blankScreen, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'blankScreen.started')
            blankScreen.setAutoDraw(True)
        if blankScreen.status == STARTED:
            # is it time to stop? (based on global clock, using actual start)
            if tThisFlipGlobal > blankScreen.tStartRefresh + 1.0 + random()-frameTolerance:
                # keep track of stop time/frame for later
                blankScreen.tStop = t  # not accounting for scr refresh
                blankScreen.frameNStop = frameN  # exact frame index
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'blankScreen.stopped')
                blankScreen.setAutoDraw(False)
        
        # check for quit (typically the Esc key)
        if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
            core.quit()
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            routineForceEnded = True
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in blankComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # --- Ending Routine "blank" ---
    for thisComponent in blankComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    # the Routine "blank" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset()
    
    # set up handler to look after randomisation of conditions etc
    trials = data.TrialHandler(nReps=5.0, method='random', 
        extraInfo=expInfo, originPath=-1,
        trialList=data.importConditions(posFile),
        seed=None, name='trials')
    thisExp.addLoop(trials)  # add the loop to the experiment
    thisTrial = trials.trialList[0]  # so we can initialise stimuli with some values
    # abbreviate parameter names if possible (e.g. rgb = thisTrial.rgb)
    if thisTrial != None:
        for paramName in thisTrial:
            exec('{} = thisTrial[paramName]'.format(paramName))
    
    for thisTrial in trials:
        currentLoop = trials
        # abbreviate parameter names if possible (e.g. rgb = thisTrial.rgb)
        if thisTrial != None:
            for paramName in thisTrial:
                exec('{} = thisTrial[paramName]'.format(paramName))
        
        # --- Prepare to start Routine "cross" ---
        continueRoutine = True
        routineForceEnded = False
        # update component parameters for each repeat
        # Run 'Begin Routine' code from codeFixate
        logging.log(level=logging.INFO, msg=f'FixationCross')
        eyetracker.setRecordingState(True)
        
        if currentLoop.name == "trials":
            eyetracker.sendMessage("Trial " + str(trials.thisN+1) + ", Condition " + trialtype)
            print("VP %s TRIALID %d CONDITION %s"%(expInfo['participant'], trials.thisN+1, trialtype))
            
        if currentLoop.name == "testtrials":
            eyetracker.sendMessage("Test-Trial " + str(testtrials.thisN+1) + ", Condition " + stimtype)
            print("VP %s TESTTRIALID %d CONDITION %s"%(expInfo['participant'], testtrials.thisN+1, stimtype))
        # keep track of which components have finished
        crossComponents = [fixcross]
        for thisComponent in crossComponents:
            thisComponent.tStart = None
            thisComponent.tStop = None
            thisComponent.tStartRefresh = None
            thisComponent.tStopRefresh = None
            if hasattr(thisComponent, 'status'):
                thisComponent.status = NOT_STARTED
        # reset timers
        t = 0
        _timeToFirstFrame = win.getFutureFlipTime(clock="now")
        frameN = -1
        
        # --- Run Routine "cross" ---
        while continueRoutine and routineTimer.getTime() < 1.0:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # *fixcross* updates
            if fixcross.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                fixcross.frameNStart = frameN  # exact frame index
                fixcross.tStart = t  # local t and not account for scr refresh
                fixcross.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(fixcross, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'fixcross.started')
                fixcross.setAutoDraw(True)
            if fixcross.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > fixcross.tStartRefresh + 1-frameTolerance:
                    # keep track of stop time/frame for later
                    fixcross.tStop = t  # not accounting for scr refresh
                    fixcross.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'fixcross.stopped')
                    fixcross.setAutoDraw(False)
            
            # check for quit (typically the Esc key)
            if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
                core.quit()
            
            # check if all components have finished
            if not continueRoutine:  # a component has requested a forced-end of Routine
                routineForceEnded = True
                break
            continueRoutine = False  # will revert to True if at least one component still running
            for thisComponent in crossComponents:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # --- Ending Routine "cross" ---
        for thisComponent in crossComponents:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        # using non-slip timing so subtract the expected duration of this Routine (unless ended on request)
        if routineForceEnded:
            routineTimer.reset()
        else:
            routineTimer.addTime(-1.000000)
        
        # --- Prepare to start Routine "blank" ---
        continueRoutine = True
        routineForceEnded = False
        # update component parameters for each repeat
        # keep track of which components have finished
        blankComponents = [blankScreen]
        for thisComponent in blankComponents:
            thisComponent.tStart = None
            thisComponent.tStop = None
            thisComponent.tStartRefresh = None
            thisComponent.tStopRefresh = None
            if hasattr(thisComponent, 'status'):
                thisComponent.status = NOT_STARTED
        # reset timers
        t = 0
        _timeToFirstFrame = win.getFutureFlipTime(clock="now")
        frameN = -1
        
        # --- Run Routine "blank" ---
        while continueRoutine:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # *blankScreen* updates
            if blankScreen.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                blankScreen.frameNStart = frameN  # exact frame index
                blankScreen.tStart = t  # local t and not account for scr refresh
                blankScreen.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(blankScreen, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'blankScreen.started')
                blankScreen.setAutoDraw(True)
            if blankScreen.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > blankScreen.tStartRefresh + 1.0 + random()-frameTolerance:
                    # keep track of stop time/frame for later
                    blankScreen.tStop = t  # not accounting for scr refresh
                    blankScreen.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'blankScreen.stopped')
                    blankScreen.setAutoDraw(False)
            
            # check for quit (typically the Esc key)
            if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
                core.quit()
            
            # check if all components have finished
            if not continueRoutine:  # a component has requested a forced-end of Routine
                routineForceEnded = True
                break
            continueRoutine = False  # will revert to True if at least one component still running
            for thisComponent in blankComponents:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # --- Ending Routine "blank" ---
        for thisComponent in blankComponents:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        # the Routine "blank" was not non-slip safe, so reset the non-slip timer
        routineTimer.reset()
        
        # --- Prepare to start Routine "trial" ---
        continueRoutine = True
        routineForceEnded = False
        # update component parameters for each repeat
        # Run 'Begin Routine' code from codeTrial
        looked_at = False
        cursorcolor="white"
        
        logging.log(level=logging.INFO, msg=f'ImageOnset_{trialtype}')
        ioServer.sendMessageEvent(text='ImageOnset')
        eyetracker.sendMessage('ImageOnset')
        image.setPos(position)
        image.setImage(eval(trialtype))
        roi.setPos(position)
        roi.setSize([image.size])
        # clear any previous roi data
        roi.reset()
        # keep track of which components have finished
        trialComponents = [image, roi, gazeCursor, portImage]
        for thisComponent in trialComponents:
            thisComponent.tStart = None
            thisComponent.tStop = None
            thisComponent.tStartRefresh = None
            thisComponent.tStopRefresh = None
            if hasattr(thisComponent, 'status'):
                thisComponent.status = NOT_STARTED
        # reset timers
        t = 0
        _timeToFirstFrame = win.getFutureFlipTime(clock="now")
        frameN = -1
        
        # --- Run Routine "trial" ---
        while continueRoutine and routineTimer.getTime() < 5.0:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            # Run 'Each Frame' code from codeTrial
            roi.size = image.size
            roi.pos = position
            
            if roi.isLookedIn:
                looked_at = True
                continueRoutine = False
                
            if event.getKeys(keyList=["p"]) and not paused:
                print("Experiment Paused - Press 'p' to continue.")
                paused = True
                event.waitKeys(keyList=["p"], clearEvents=True)
                
                print("Experiment Continued")
                paused = False
                event.clearEvents(eventType="keyboard")
                
            
            # *image* updates
            if image.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                # keep track of start time/frame for later
                image.frameNStart = frameN  # exact frame index
                image.tStart = t  # local t and not account for scr refresh
                image.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(image, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'image.started')
                image.setAutoDraw(True)
            if image.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > image.tStartRefresh + 5-frameTolerance:
                    # keep track of stop time/frame for later
                    image.tStop = t  # not accounting for scr refresh
                    image.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'image.stopped')
                    image.setAutoDraw(False)
            if roi.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                roi.frameNStart = frameN  # exact frame index
                roi.tStart = t  # local t and not account for scr refresh
                roi.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(roi, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'roi.started')
                roi.status = STARTED
            if roi.status == STARTED:
                # check whether roi has been looked in
                if roi.isLookedIn:
                    if not roi.wasLookedIn:
                        roi.timesOn.append(routineTimer.getTime()) # store time of first look
                        roi.timesOff.append(routineTimer.getTime()) # store time looked until
                    else:
                        roi.timesOff[-1] = routineTimer.getTime() # update time looked until
                    roi.wasLookedIn = True  # if roi is still looked at next frame, it is not a new look
                else:
                    if roi.wasLookedIn:
                        roi.timesOff[-1] = routineTimer.getTime() # update time looked until
                    roi.wasLookedIn = False  # if roi is looked at next frame, it is a new look
            else:
                roi.clock.reset() # keep clock at 0 if roi hasn't started / has finished
                roi.wasLookedIn = False  # if roi is looked at next frame, it is a new look
            if roi.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > roi.tStartRefresh + 5-frameTolerance:
                    # keep track of stop time/frame for later
                    roi.tStop = t  # not accounting for scr refresh
                    roi.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'roi.stopped')
                    roi.status = FINISHED
            
            # *gazeCursor* updates
            if gazeCursor.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                gazeCursor.frameNStart = frameN  # exact frame index
                gazeCursor.tStart = t  # local t and not account for scr refresh
                gazeCursor.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(gazeCursor, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'gazeCursor.started')
                gazeCursor.setAutoDraw(True)
            if gazeCursor.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > gazeCursor.tStartRefresh + 5-frameTolerance:
                    # keep track of stop time/frame for later
                    gazeCursor.tStop = t  # not accounting for scr refresh
                    gazeCursor.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'gazeCursor.stopped')
                    gazeCursor.setAutoDraw(False)
            if gazeCursor.status == STARTED:  # only update if drawing
                gazeCursor.setFillColor(cursorcolor, log=False)
                gazeCursor.setOpacity(0.0, log=False)
                gazeCursor.setPos([eyetracker.getPos()], log=False)
            # *portImage* updates
            if portImage.status == NOT_STARTED and t >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                portImage.frameNStart = frameN  # exact frame index
                portImage.tStart = t  # local t and not account for scr refresh
                portImage.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(portImage, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.addData('portImage.started', t)
                portImage.status = STARTED
                win.callOnFlip(portImage.setData, int(1))
            if portImage.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > portImage.tStartRefresh + 0.005-frameTolerance:
                    # keep track of stop time/frame for later
                    portImage.tStop = t  # not accounting for scr refresh
                    portImage.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.addData('portImage.stopped', t)
                    portImage.status = FINISHED
                    win.callOnFlip(portImage.setData, int(0))
            
            # check for quit (typically the Esc key)
            if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
                core.quit()
            
            # check if all components have finished
            if not continueRoutine:  # a component has requested a forced-end of Routine
                routineForceEnded = True
                break
            continueRoutine = False  # will revert to True if at least one component still running
            for thisComponent in trialComponents:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # --- Ending Routine "trial" ---
        for thisComponent in trialComponents:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        # Run 'End Routine' code from codeTrial
        if looked_at & ("plus" in trialtype):
            feedback_color = "red"
            feedback_opacity = 1
            cursorcolor="red"
            feedback_audio = "audio/error.wav"
            port_msg = 128 + 4  # Shock1
            log_msg = "lookedat_cs_plus"
            points = 0  # 5
            # score -= points
            feedback_points = ""  # f"- {points}"
            feedback_score = ""  # f"Score: {score}"
        elif looked_at & ("minus" in trialtype):
            feedback_color = "green"
            feedback_opacity = 1
            cursorcolor="green"
            feedback_audio = "audio/win.wav"
            port_msg = 8  # Reward
            log_msg = "lookedat_cs_minus"
            points = 5
            score += points
            feedback_points = f"+ {points}"
            feedback_score = f"Score: {score}"
        else:
            feedback_color = (0,0,0)
            feedback_opacity = 0
            feedback_audio = "audio/silence.wav"
            port_msg = 16  # No Feedback
            log_msg = "no_look"
            points = 0
            feedback_points = ""
            feedback_score = ""
        
        rectsize = [item * 1.02 for item in image.size]
        
        image_w = image.size[0]
        image_h = image.size[1]
        imagesize_test = [image_w * 1.2, image_h * 1.2]
        trials.addData('roi.numLooks', roi.numLooks)
        if roi.numLooks:
           trials.addData('roi.timesOn', roi.timesOn)
           trials.addData('roi.timesOff', roi.timesOff)
        else:
           trials.addData('roi.timesOn', "")
           trials.addData('roi.timesOff', "")
        if portImage.status == STARTED:
            win.callOnFlip(portImage.setData, int(0))
        # using non-slip timing so subtract the expected duration of this Routine (unless ended on request)
        if routineForceEnded:
            routineTimer.reset()
        else:
            routineTimer.addTime(-5.000000)
        
        # --- Prepare to start Routine "feedback" ---
        continueRoutine = True
        routineForceEnded = False
        # update component parameters for each repeat
        # Run 'Begin Routine' code from codeFeedback
        logging.log(level=logging.INFO, msg=f'FeedbackOnset_{log_msg}')
        ioServer.sendMessageEvent(text='FeedbackOnset')
        eyetracker.sendMessage('FeedbackOnset')
        polygon.setOpacity(feedback_opacity)
        polygon.setPos(position)
        polygon.setLineColor(feedback_color)
        imageFeedback.setPos(position)
        imageFeedback.setImage(eval(trialtype))
        soundFeedback.setSound(feedback_audio, secs=1, hamming=True)
        soundFeedback.setVolume(1.0, log=False)
        # keep track of which components have finished
        feedbackComponents = [textPoints, textScore, polygon, imageFeedback, gazeCursor_Feedback, soundFeedback, portShock1, portShock2, portShock3]
        for thisComponent in feedbackComponents:
            thisComponent.tStart = None
            thisComponent.tStop = None
            thisComponent.tStartRefresh = None
            thisComponent.tStopRefresh = None
            if hasattr(thisComponent, 'status'):
                thisComponent.status = NOT_STARTED
        # reset timers
        t = 0
        _timeToFirstFrame = win.getFutureFlipTime(clock="now")
        frameN = -1
        
        # --- Run Routine "feedback" ---
        while continueRoutine and routineTimer.getTime() < 3.0:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # *textPoints* updates
            if textPoints.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                # keep track of start time/frame for later
                textPoints.frameNStart = frameN  # exact frame index
                textPoints.tStart = t  # local t and not account for scr refresh
                textPoints.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(textPoints, 'tStartRefresh')  # time at next scr refresh
                textPoints.setAutoDraw(True)
            if textPoints.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > textPoints.tStartRefresh + 3-frameTolerance:
                    # keep track of stop time/frame for later
                    textPoints.tStop = t  # not accounting for scr refresh
                    textPoints.frameNStop = frameN  # exact frame index
                    textPoints.setAutoDraw(False)
            if textPoints.status == STARTED:  # only update if drawing
                textPoints.setText(feedback_points, log=False)
            
            # *textScore* updates
            if textScore.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                # keep track of start time/frame for later
                textScore.frameNStart = frameN  # exact frame index
                textScore.tStart = t  # local t and not account for scr refresh
                textScore.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(textScore, 'tStartRefresh')  # time at next scr refresh
                textScore.setAutoDraw(True)
            if textScore.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > textScore.tStartRefresh + 3-frameTolerance:
                    # keep track of stop time/frame for later
                    textScore.tStop = t  # not accounting for scr refresh
                    textScore.frameNStop = frameN  # exact frame index
                    textScore.setAutoDraw(False)
            if textScore.status == STARTED:  # only update if drawing
                textScore.setText(feedback_score, log=False)
            
            # *polygon* updates
            if polygon.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                polygon.frameNStart = frameN  # exact frame index
                polygon.tStart = t  # local t and not account for scr refresh
                polygon.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(polygon, 'tStartRefresh')  # time at next scr refresh
                polygon.setAutoDraw(True)
            if polygon.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > polygon.tStartRefresh + 3-frameTolerance:
                    # keep track of stop time/frame for later
                    polygon.tStop = t  # not accounting for scr refresh
                    polygon.frameNStop = frameN  # exact frame index
                    polygon.setAutoDraw(False)
            if polygon.status == STARTED:  # only update if drawing
                polygon.setSize(rectsize, log=False)
            
            # *imageFeedback* updates
            if imageFeedback.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                # keep track of start time/frame for later
                imageFeedback.frameNStart = frameN  # exact frame index
                imageFeedback.tStart = t  # local t and not account for scr refresh
                imageFeedback.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(imageFeedback, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'imageFeedback.started')
                imageFeedback.setAutoDraw(True)
            if imageFeedback.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > imageFeedback.tStartRefresh + 3-frameTolerance:
                    # keep track of stop time/frame for later
                    imageFeedback.tStop = t  # not accounting for scr refresh
                    imageFeedback.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'imageFeedback.stopped')
                    imageFeedback.setAutoDraw(False)
            
            # *gazeCursor_Feedback* updates
            if gazeCursor_Feedback.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                gazeCursor_Feedback.frameNStart = frameN  # exact frame index
                gazeCursor_Feedback.tStart = t  # local t and not account for scr refresh
                gazeCursor_Feedback.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(gazeCursor_Feedback, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'gazeCursor_Feedback.started')
                gazeCursor_Feedback.setAutoDraw(True)
            if gazeCursor_Feedback.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > gazeCursor_Feedback.tStartRefresh + 3-frameTolerance:
                    # keep track of stop time/frame for later
                    gazeCursor_Feedback.tStop = t  # not accounting for scr refresh
                    gazeCursor_Feedback.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'gazeCursor_Feedback.stopped')
                    gazeCursor_Feedback.setAutoDraw(False)
            if gazeCursor_Feedback.status == STARTED:  # only update if drawing
                gazeCursor_Feedback.setFillColor(cursorcolor, log=False)
                gazeCursor_Feedback.setOpacity(0.0, log=False)
                gazeCursor_Feedback.setPos([eyetracker.getPos()], log=False)
            # start/stop soundFeedback
            if soundFeedback.status == NOT_STARTED and t >= 0-frameTolerance:
                # keep track of start time/frame for later
                soundFeedback.frameNStart = frameN  # exact frame index
                soundFeedback.tStart = t  # local t and not account for scr refresh
                soundFeedback.tStartRefresh = tThisFlipGlobal  # on global time
                # add timestamp to datafile
                thisExp.addData('soundFeedback.started', t)
                soundFeedback.play()  # start the sound (it finishes automatically)
            if soundFeedback.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > soundFeedback.tStartRefresh + 1-frameTolerance:
                    # keep track of stop time/frame for later
                    soundFeedback.tStop = t  # not accounting for scr refresh
                    soundFeedback.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.addData('soundFeedback.stopped', t)
                    soundFeedback.stop()
            # *portShock1* updates
            if portShock1.status == NOT_STARTED and t >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                portShock1.frameNStart = frameN  # exact frame index
                portShock1.tStart = t  # local t and not account for scr refresh
                portShock1.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(portShock1, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.addData('portShock1.started', t)
                portShock1.status = STARTED
                win.callOnFlip(portShock1.setData, int(port_msg))
            if portShock1.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > portShock1.tStartRefresh + 0.005-frameTolerance:
                    # keep track of stop time/frame for later
                    portShock1.tStop = t  # not accounting for scr refresh
                    portShock1.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.addData('portShock1.stopped', t)
                    portShock1.status = FINISHED
                    win.callOnFlip(portShock1.setData, int(0))
            # *portShock2* updates
            if portShock2.status == NOT_STARTED and t >= 0.05-frameTolerance:
                # keep track of start time/frame for later
                portShock2.frameNStart = frameN  # exact frame index
                portShock2.tStart = t  # local t and not account for scr refresh
                portShock2.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(portShock2, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.addData('portShock2.started', t)
                portShock2.status = STARTED
                win.callOnFlip(portShock2.setData, int(128 if port_msg == 132 else 0))
            if portShock2.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > portShock2.tStartRefresh + 0.005-frameTolerance:
                    # keep track of stop time/frame for later
                    portShock2.tStop = t  # not accounting for scr refresh
                    portShock2.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.addData('portShock2.stopped', t)
                    portShock2.status = FINISHED
                    win.callOnFlip(portShock2.setData, int(0))
            # *portShock3* updates
            if portShock3.status == NOT_STARTED and t >= 0.1-frameTolerance:
                # keep track of start time/frame for later
                portShock3.frameNStart = frameN  # exact frame index
                portShock3.tStart = t  # local t and not account for scr refresh
                portShock3.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(portShock3, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.addData('portShock3.started', t)
                portShock3.status = STARTED
                win.callOnFlip(portShock3.setData, int(128 if port_msg == 132 else 0))
            if portShock3.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > portShock3.tStartRefresh + 0.005-frameTolerance:
                    # keep track of stop time/frame for later
                    portShock3.tStop = t  # not accounting for scr refresh
                    portShock3.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.addData('portShock3.stopped', t)
                    portShock3.status = FINISHED
                    win.callOnFlip(portShock3.setData, int(0))
            
            # check for quit (typically the Esc key)
            if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
                core.quit()
            
            # check if all components have finished
            if not continueRoutine:  # a component has requested a forced-end of Routine
                routineForceEnded = True
                break
            continueRoutine = False  # will revert to True if at least one component still running
            for thisComponent in feedbackComponents:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # --- Ending Routine "feedback" ---
        for thisComponent in feedbackComponents:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        soundFeedback.stop()  # ensure sound has stopped at end of routine
        if portShock1.status == STARTED:
            win.callOnFlip(portShock1.setData, int(0))
        if portShock2.status == STARTED:
            win.callOnFlip(portShock2.setData, int(0))
        if portShock3.status == STARTED:
            win.callOnFlip(portShock3.setData, int(0))
        # using non-slip timing so subtract the expected duration of this Routine (unless ended on request)
        if routineForceEnded:
            routineTimer.reset()
        else:
            routineTimer.addTime(-3.000000)
        
        # --- Prepare to start Routine "crossEnd" ---
        continueRoutine = True
        routineForceEnded = False
        # update component parameters for each repeat
        # Run 'Begin Routine' code from code_end
        eyetracker.setRecordingState(False)
        # keep track of which components have finished
        crossEndComponents = [fixcrossE]
        for thisComponent in crossEndComponents:
            thisComponent.tStart = None
            thisComponent.tStop = None
            thisComponent.tStartRefresh = None
            thisComponent.tStopRefresh = None
            if hasattr(thisComponent, 'status'):
                thisComponent.status = NOT_STARTED
        # reset timers
        t = 0
        _timeToFirstFrame = win.getFutureFlipTime(clock="now")
        frameN = -1
        
        # --- Run Routine "crossEnd" ---
        while continueRoutine:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # *fixcrossE* updates
            if fixcrossE.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                fixcrossE.frameNStart = frameN  # exact frame index
                fixcrossE.tStart = t  # local t and not account for scr refresh
                fixcrossE.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(fixcrossE, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'fixcrossE.started')
                fixcrossE.setAutoDraw(True)
            if fixcrossE.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > fixcrossE.tStartRefresh +  3 + random()-frameTolerance:
                    # keep track of stop time/frame for later
                    fixcrossE.tStop = t  # not accounting for scr refresh
                    fixcrossE.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'fixcrossE.stopped')
                    fixcrossE.setAutoDraw(False)
            
            # check for quit (typically the Esc key)
            if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
                core.quit()
            
            # check if all components have finished
            if not continueRoutine:  # a component has requested a forced-end of Routine
                routineForceEnded = True
                break
            continueRoutine = False  # will revert to True if at least one component still running
            for thisComponent in crossEndComponents:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # --- Ending Routine "crossEnd" ---
        for thisComponent in crossEndComponents:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        # the Routine "crossEnd" was not non-slip safe, so reset the non-slip timer
        routineTimer.reset()
        thisExp.nextEntry()
        
    # completed 5.0 repeats of 'trials'
    
    
    # --- Prepare to start Routine "painRating" ---
    continueRoutine = True
    routineForceEnded = False
    # update component parameters for each repeat
    sliderPain.reset()
    spacePain.keys = []
    spacePain.rt = []
    _spacePain_allKeys = []
    # Run 'Begin Routine' code from codePainRating
    logging.log(level=logging.INFO, msg=f'PainRating')
    print("Pain Rating. Press Space to continue.")
    
    values = [None]
    
    win.mouseVisible = True
    # keep track of which components have finished
    painRatingComponents = [textRate, text_noPain, text_slightPain, text_highPain, sliderPain, textSpacePain, spacePain]
    for thisComponent in painRatingComponents:
        thisComponent.tStart = None
        thisComponent.tStop = None
        thisComponent.tStartRefresh = None
        thisComponent.tStopRefresh = None
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    # reset timers
    t = 0
    _timeToFirstFrame = win.getFutureFlipTime(clock="now")
    frameN = -1
    
    # --- Run Routine "painRating" ---
    while continueRoutine:
        # get current time
        t = routineTimer.getTime()
        tThisFlip = win.getFutureFlipTime(clock=routineTimer)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        
        # *textRate* updates
        if textRate.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            textRate.frameNStart = frameN  # exact frame index
            textRate.tStart = t  # local t and not account for scr refresh
            textRate.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(textRate, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'textRate.started')
            textRate.setAutoDraw(True)
        
        # *text_noPain* updates
        if text_noPain.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            text_noPain.frameNStart = frameN  # exact frame index
            text_noPain.tStart = t  # local t and not account for scr refresh
            text_noPain.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(text_noPain, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'text_noPain.started')
            text_noPain.setAutoDraw(True)
        
        # *text_slightPain* updates
        if text_slightPain.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            text_slightPain.frameNStart = frameN  # exact frame index
            text_slightPain.tStart = t  # local t and not account for scr refresh
            text_slightPain.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(text_slightPain, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'text_slightPain.started')
            text_slightPain.setAutoDraw(True)
        
        # *text_highPain* updates
        if text_highPain.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            text_highPain.frameNStart = frameN  # exact frame index
            text_highPain.tStart = t  # local t and not account for scr refresh
            text_highPain.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(text_highPain, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'text_highPain.started')
            text_highPain.setAutoDraw(True)
        
        # *sliderPain* updates
        if sliderPain.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            sliderPain.frameNStart = frameN  # exact frame index
            sliderPain.tStart = t  # local t and not account for scr refresh
            sliderPain.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(sliderPain, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'sliderPain.started')
            sliderPain.setAutoDraw(True)
        
        # *textSpacePain* updates
        if textSpacePain.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            textSpacePain.frameNStart = frameN  # exact frame index
            textSpacePain.tStart = t  # local t and not account for scr refresh
            textSpacePain.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(textSpacePain, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'textSpacePain.started')
            textSpacePain.setAutoDraw(True)
        
        # *spacePain* updates
        waitOnFlip = False
        if spacePain.status == NOT_STARTED and sliderPain.rating:
            # keep track of start time/frame for later
            spacePain.frameNStart = frameN  # exact frame index
            spacePain.tStart = t  # local t and not account for scr refresh
            spacePain.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(spacePain, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'spacePain.started')
            spacePain.status = STARTED
            # keyboard checking is just starting
            waitOnFlip = True
            win.callOnFlip(spacePain.clock.reset)  # t=0 on next screen flip
            win.callOnFlip(spacePain.clearEvents, eventType='keyboard')  # clear events on next screen flip
        if spacePain.status == STARTED and not waitOnFlip:
            theseKeys = spacePain.getKeys(keyList=['space', 'enter'], waitRelease=False)
            _spacePain_allKeys.extend(theseKeys)
            if len(_spacePain_allKeys):
                spacePain.keys = _spacePain_allKeys[-1].name  # just the last key pressed
                spacePain.rt = _spacePain_allKeys[-1].rt
                # a response ends the routine
                continueRoutine = False
        # Run 'Each Frame' code from codePainRating
        if sliderPainT.getRating() != values[-1]:
            values.append(sliderPainT.getRating())
            print(f"Value: {values[-1]}")
        
        # check for quit (typically the Esc key)
        if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
            core.quit()
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            routineForceEnded = True
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in painRatingComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # --- Ending Routine "painRating" ---
    for thisComponent in painRatingComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    blocks.addData('sliderPain.response', sliderPain.getRating())
    blocks.addData('sliderPain.rt', sliderPain.getRT())
    # check responses
    if spacePain.keys in ['', [], None]:  # No response was made
        spacePain.keys = None
    blocks.addData('spacePain.keys',spacePain.keys)
    if spacePain.keys != None:  # we had a response
        blocks.addData('spacePain.rt', spacePain.rt)
    # the Routine "painRating" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset()
    
    # set up handler to look after randomisation of conditions etc
    ratingtrials2 = data.TrialHandler(nReps=1.0, method='random', 
        extraInfo=expInfo, originPath=-1,
        trialList=data.importConditions(stimFile),
        seed=None, name='ratingtrials2')
    thisExp.addLoop(ratingtrials2)  # add the loop to the experiment
    thisRatingtrials2 = ratingtrials2.trialList[0]  # so we can initialise stimuli with some values
    # abbreviate parameter names if possible (e.g. rgb = thisRatingtrials2.rgb)
    if thisRatingtrials2 != None:
        for paramName in thisRatingtrials2:
            exec('{} = thisRatingtrials2[paramName]'.format(paramName))
    
    for thisRatingtrials2 in ratingtrials2:
        currentLoop = ratingtrials2
        # abbreviate parameter names if possible (e.g. rgb = thisRatingtrials2.rgb)
        if thisRatingtrials2 != None:
            for paramName in thisRatingtrials2:
                exec('{} = thisRatingtrials2[paramName]'.format(paramName))
        
        # --- Prepare to start Routine "stimRating" ---
        continueRoutine = True
        routineForceEnded = False
        # update component parameters for each repeat
        # Run 'Begin Routine' code from codeStimRating
        logging.log(level=logging.INFO, msg=f'StimRating_{stimtype}')
        print("Rating of Stimulus: %s"%(stimtype))
        win.mouseVisible = True
        imageRating.setPos((0, 0.2))
        imageRating.setImage(eval(stimtype))
        sliderStim.reset()
        spaceStim.keys = []
        spaceStim.rt = []
        _spaceStim_allKeys = []
        # keep track of which components have finished
        stimRatingComponents = [textRateStim, imageRating, textUnpleasant, textPleasant, sliderStim, textSpaceStim, spaceStim]
        for thisComponent in stimRatingComponents:
            thisComponent.tStart = None
            thisComponent.tStop = None
            thisComponent.tStartRefresh = None
            thisComponent.tStopRefresh = None
            if hasattr(thisComponent, 'status'):
                thisComponent.status = NOT_STARTED
        # reset timers
        t = 0
        _timeToFirstFrame = win.getFutureFlipTime(clock="now")
        frameN = -1
        
        # --- Run Routine "stimRating" ---
        while continueRoutine:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # *textRateStim* updates
            if textRateStim.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                textRateStim.frameNStart = frameN  # exact frame index
                textRateStim.tStart = t  # local t and not account for scr refresh
                textRateStim.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(textRateStim, 'tStartRefresh')  # time at next scr refresh
                textRateStim.setAutoDraw(True)
            
            # *imageRating* updates
            if imageRating.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                # keep track of start time/frame for later
                imageRating.frameNStart = frameN  # exact frame index
                imageRating.tStart = t  # local t and not account for scr refresh
                imageRating.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(imageRating, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'imageRating.started')
                imageRating.setAutoDraw(True)
            
            # *textUnpleasant* updates
            if textUnpleasant.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                textUnpleasant.frameNStart = frameN  # exact frame index
                textUnpleasant.tStart = t  # local t and not account for scr refresh
                textUnpleasant.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(textUnpleasant, 'tStartRefresh')  # time at next scr refresh
                textUnpleasant.setAutoDraw(True)
            
            # *textPleasant* updates
            if textPleasant.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                textPleasant.frameNStart = frameN  # exact frame index
                textPleasant.tStart = t  # local t and not account for scr refresh
                textPleasant.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(textPleasant, 'tStartRefresh')  # time at next scr refresh
                textPleasant.setAutoDraw(True)
            
            # *sliderStim* updates
            if sliderStim.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                sliderStim.frameNStart = frameN  # exact frame index
                sliderStim.tStart = t  # local t and not account for scr refresh
                sliderStim.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sliderStim, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'sliderStim.started')
                sliderStim.setAutoDraw(True)
            
            # *textSpaceStim* updates
            if textSpaceStim.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                textSpaceStim.frameNStart = frameN  # exact frame index
                textSpaceStim.tStart = t  # local t and not account for scr refresh
                textSpaceStim.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(textSpaceStim, 'tStartRefresh')  # time at next scr refresh
                textSpaceStim.setAutoDraw(True)
            
            # *spaceStim* updates
            waitOnFlip = False
            if spaceStim.status == NOT_STARTED and sliderStim.rating:
                # keep track of start time/frame for later
                spaceStim.frameNStart = frameN  # exact frame index
                spaceStim.tStart = t  # local t and not account for scr refresh
                spaceStim.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(spaceStim, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'spaceStim.started')
                spaceStim.status = STARTED
                # keyboard checking is just starting
                waitOnFlip = True
                win.callOnFlip(spaceStim.clock.reset)  # t=0 on next screen flip
                win.callOnFlip(spaceStim.clearEvents, eventType='keyboard')  # clear events on next screen flip
            if spaceStim.status == STARTED and not waitOnFlip:
                theseKeys = spaceStim.getKeys(keyList=['space', 'enter'], waitRelease=False)
                _spaceStim_allKeys.extend(theseKeys)
                if len(_spaceStim_allKeys):
                    spaceStim.keys = _spaceStim_allKeys[-1].name  # just the last key pressed
                    spaceStim.rt = _spaceStim_allKeys[-1].rt
                    # a response ends the routine
                    continueRoutine = False
            
            # check for quit (typically the Esc key)
            if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
                core.quit()
            
            # check if all components have finished
            if not continueRoutine:  # a component has requested a forced-end of Routine
                routineForceEnded = True
                break
            continueRoutine = False  # will revert to True if at least one component still running
            for thisComponent in stimRatingComponents:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # --- Ending Routine "stimRating" ---
        for thisComponent in stimRatingComponents:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        ratingtrials2.addData('sliderStim.response', sliderStim.getRating())
        ratingtrials2.addData('sliderStim.rt', sliderStim.getRT())
        # check responses
        if spaceStim.keys in ['', [], None]:  # No response was made
            spaceStim.keys = None
        ratingtrials2.addData('spaceStim.keys',spaceStim.keys)
        if spaceStim.keys != None:  # we had a response
            ratingtrials2.addData('spaceStim.rt', spaceStim.rt)
        # the Routine "stimRating" was not non-slip safe, so reset the non-slip timer
        routineTimer.reset()
        thisExp.nextEntry()
        
    # completed 1.0 repeats of 'ratingtrials2'
    
    
    # --- Prepare to start Routine "testInstr" ---
    continueRoutine = True
    routineForceEnded = False
    # update component parameters for each repeat
    spaceTestInstr.keys = []
    spaceTestInstr.rt = []
    _spaceTestInstr_allKeys = []
    # Run 'Begin Routine' code from codeTestInstruction
    print("Start Test Phase. Press Space to continue.")
    win.mouseVisible = False
    # keep track of which components have finished
    testInstrComponents = [textTestInstr, spaceTestInstr]
    for thisComponent in testInstrComponents:
        thisComponent.tStart = None
        thisComponent.tStop = None
        thisComponent.tStartRefresh = None
        thisComponent.tStopRefresh = None
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    # reset timers
    t = 0
    _timeToFirstFrame = win.getFutureFlipTime(clock="now")
    frameN = -1
    
    # --- Run Routine "testInstr" ---
    while continueRoutine:
        # get current time
        t = routineTimer.getTime()
        tThisFlip = win.getFutureFlipTime(clock=routineTimer)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        
        # *textTestInstr* updates
        if textTestInstr.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            textTestInstr.frameNStart = frameN  # exact frame index
            textTestInstr.tStart = t  # local t and not account for scr refresh
            textTestInstr.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(textTestInstr, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'textTestInstr.started')
            textTestInstr.setAutoDraw(True)
        
        # *spaceTestInstr* updates
        waitOnFlip = False
        if spaceTestInstr.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            spaceTestInstr.frameNStart = frameN  # exact frame index
            spaceTestInstr.tStart = t  # local t and not account for scr refresh
            spaceTestInstr.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(spaceTestInstr, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'spaceTestInstr.started')
            spaceTestInstr.status = STARTED
            # keyboard checking is just starting
            waitOnFlip = True
            win.callOnFlip(spaceTestInstr.clock.reset)  # t=0 on next screen flip
            win.callOnFlip(spaceTestInstr.clearEvents, eventType='keyboard')  # clear events on next screen flip
        if spaceTestInstr.status == STARTED and not waitOnFlip:
            theseKeys = spaceTestInstr.getKeys(keyList=['space', 'enter'], waitRelease=False)
            _spaceTestInstr_allKeys.extend(theseKeys)
            if len(_spaceTestInstr_allKeys):
                spaceTestInstr.keys = _spaceTestInstr_allKeys[-1].name  # just the last key pressed
                spaceTestInstr.rt = _spaceTestInstr_allKeys[-1].rt
                # a response ends the routine
                continueRoutine = False
        
        # check for quit (typically the Esc key)
        if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
            core.quit()
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            routineForceEnded = True
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in testInstrComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # --- Ending Routine "testInstr" ---
    for thisComponent in testInstrComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    # check responses
    if spaceTestInstr.keys in ['', [], None]:  # No response was made
        spaceTestInstr.keys = None
    blocks.addData('spaceTestInstr.keys',spaceTestInstr.keys)
    if spaceTestInstr.keys != None:  # we had a response
        blocks.addData('spaceTestInstr.rt', spaceTestInstr.rt)
    # the Routine "testInstr" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset()
    
    # set up handler to look after randomisation of conditions etc
    testtrials = data.TrialHandler(nReps=10.0, method='random', 
        extraInfo=expInfo, originPath=-1,
        trialList=data.importConditions(stimFile),
        seed=None, name='testtrials')
    thisExp.addLoop(testtrials)  # add the loop to the experiment
    thisTesttrial = testtrials.trialList[0]  # so we can initialise stimuli with some values
    # abbreviate parameter names if possible (e.g. rgb = thisTesttrial.rgb)
    if thisTesttrial != None:
        for paramName in thisTesttrial:
            exec('{} = thisTesttrial[paramName]'.format(paramName))
    
    for thisTesttrial in testtrials:
        currentLoop = testtrials
        # abbreviate parameter names if possible (e.g. rgb = thisTesttrial.rgb)
        if thisTesttrial != None:
            for paramName in thisTesttrial:
                exec('{} = thisTesttrial[paramName]'.format(paramName))
        
        # --- Prepare to start Routine "cross" ---
        continueRoutine = True
        routineForceEnded = False
        # update component parameters for each repeat
        # Run 'Begin Routine' code from codeFixate
        logging.log(level=logging.INFO, msg=f'FixationCross')
        eyetracker.setRecordingState(True)
        
        if currentLoop.name == "trials":
            eyetracker.sendMessage("Trial " + str(trials.thisN+1) + ", Condition " + trialtype)
            print("VP %s TRIALID %d CONDITION %s"%(expInfo['participant'], trials.thisN+1, trialtype))
            
        if currentLoop.name == "testtrials":
            eyetracker.sendMessage("Test-Trial " + str(testtrials.thisN+1) + ", Condition " + stimtype)
            print("VP %s TESTTRIALID %d CONDITION %s"%(expInfo['participant'], testtrials.thisN+1, stimtype))
        # keep track of which components have finished
        crossComponents = [fixcross]
        for thisComponent in crossComponents:
            thisComponent.tStart = None
            thisComponent.tStop = None
            thisComponent.tStartRefresh = None
            thisComponent.tStopRefresh = None
            if hasattr(thisComponent, 'status'):
                thisComponent.status = NOT_STARTED
        # reset timers
        t = 0
        _timeToFirstFrame = win.getFutureFlipTime(clock="now")
        frameN = -1
        
        # --- Run Routine "cross" ---
        while continueRoutine and routineTimer.getTime() < 1.0:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # *fixcross* updates
            if fixcross.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                fixcross.frameNStart = frameN  # exact frame index
                fixcross.tStart = t  # local t and not account for scr refresh
                fixcross.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(fixcross, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'fixcross.started')
                fixcross.setAutoDraw(True)
            if fixcross.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > fixcross.tStartRefresh + 1-frameTolerance:
                    # keep track of stop time/frame for later
                    fixcross.tStop = t  # not accounting for scr refresh
                    fixcross.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'fixcross.stopped')
                    fixcross.setAutoDraw(False)
            
            # check for quit (typically the Esc key)
            if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
                core.quit()
            
            # check if all components have finished
            if not continueRoutine:  # a component has requested a forced-end of Routine
                routineForceEnded = True
                break
            continueRoutine = False  # will revert to True if at least one component still running
            for thisComponent in crossComponents:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # --- Ending Routine "cross" ---
        for thisComponent in crossComponents:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        # using non-slip timing so subtract the expected duration of this Routine (unless ended on request)
        if routineForceEnded:
            routineTimer.reset()
        else:
            routineTimer.addTime(-1.000000)
        
        # --- Prepare to start Routine "testtrial" ---
        continueRoutine = True
        routineForceEnded = False
        # update component parameters for each repeat
        # Run 'Begin Routine' code from codeTesttrial
        from psychopy import logging
        logging.log(level=logging.INFO, msg=f'TestImageOnset_{stimtype}')
        ioServer.sendMessageEvent(text='TestImageOnset')
        eyetracker.sendMessage('TestImageOnset')
        imageTest.setPos((0, 0))
        imageTest.setImage(eval(stimtype))
        roiTest.setSize([imageTest.size])
        # clear any previous roi data
        roiTest.reset()
        # keep track of which components have finished
        testtrialComponents = [imageTest, roiTest, portTestImage]
        for thisComponent in testtrialComponents:
            thisComponent.tStart = None
            thisComponent.tStop = None
            thisComponent.tStartRefresh = None
            thisComponent.tStopRefresh = None
            if hasattr(thisComponent, 'status'):
                thisComponent.status = NOT_STARTED
        # reset timers
        t = 0
        _timeToFirstFrame = win.getFutureFlipTime(clock="now")
        frameN = -1
        
        # --- Run Routine "testtrial" ---
        while continueRoutine and routineTimer.getTime() < 10.0:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            # Run 'Each Frame' code from codeTesttrial
            roiTest.size = imageTest.size
            
            if defaultKeyboard.getKeys(keyList=["p"]) and not paused:
                print("Experiment Paused - Press 'p' to continue.")
                paused = True
                event.waitKeys()
            if defaultKeyboard.getKeys(keyList=["p"]) and paused:
                print("Experiment Continued")
                paused = False
            
            # *imageTest* updates
            if imageTest.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                # keep track of start time/frame for later
                imageTest.frameNStart = frameN  # exact frame index
                imageTest.tStart = t  # local t and not account for scr refresh
                imageTest.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(imageTest, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'imageTest.started')
                imageTest.setAutoDraw(True)
            if imageTest.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > imageTest.tStartRefresh + 10-frameTolerance:
                    # keep track of stop time/frame for later
                    imageTest.tStop = t  # not accounting for scr refresh
                    imageTest.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'imageTest.stopped')
                    imageTest.setAutoDraw(False)
            if imageTest.status == STARTED:  # only update if drawing
                imageTest.setSize(imagesize_test, log=False)
            if roiTest.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                roiTest.frameNStart = frameN  # exact frame index
                roiTest.tStart = t  # local t and not account for scr refresh
                roiTest.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(roiTest, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'roiTest.started')
                roiTest.status = STARTED
            if roiTest.status == STARTED:
                # check whether roiTest has been looked in
                if roiTest.isLookedIn:
                    if not roiTest.wasLookedIn:
                        roiTest.timesOn.append(routineTimer.getTime()) # store time of first look
                        roiTest.timesOff.append(routineTimer.getTime()) # store time looked until
                    else:
                        roiTest.timesOff[-1] = routineTimer.getTime() # update time looked until
                    roiTest.wasLookedIn = True  # if roiTest is still looked at next frame, it is not a new look
                else:
                    if roiTest.wasLookedIn:
                        roiTest.timesOff[-1] = routineTimer.getTime() # update time looked until
                    roiTest.wasLookedIn = False  # if roiTest is looked at next frame, it is a new look
            else:
                roiTest.clock.reset() # keep clock at 0 if roi hasn't started / has finished
                roiTest.wasLookedIn = False  # if roiTest is looked at next frame, it is a new look
            if roiTest.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > roiTest.tStartRefresh + 10-frameTolerance:
                    # keep track of stop time/frame for later
                    roiTest.tStop = t  # not accounting for scr refresh
                    roiTest.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'roiTest.stopped')
                    roiTest.status = FINISHED
            # *portTestImage* updates
            if portTestImage.status == NOT_STARTED and t >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                portTestImage.frameNStart = frameN  # exact frame index
                portTestImage.tStart = t  # local t and not account for scr refresh
                portTestImage.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(portTestImage, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.addData('portTestImage.started', t)
                portTestImage.status = STARTED
                win.callOnFlip(portTestImage.setData, int(2))
            if portTestImage.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > portTestImage.tStartRefresh + 0.005-frameTolerance:
                    # keep track of stop time/frame for later
                    portTestImage.tStop = t  # not accounting for scr refresh
                    portTestImage.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.addData('portTestImage.stopped', t)
                    portTestImage.status = FINISHED
                    win.callOnFlip(portTestImage.setData, int(0))
            
            # check for quit (typically the Esc key)
            if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
                core.quit()
            
            # check if all components have finished
            if not continueRoutine:  # a component has requested a forced-end of Routine
                routineForceEnded = True
                break
            continueRoutine = False  # will revert to True if at least one component still running
            for thisComponent in testtrialComponents:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # --- Ending Routine "testtrial" ---
        for thisComponent in testtrialComponents:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        testtrials.addData('roiTest.numLooks', roiTest.numLooks)
        if roiTest.numLooks:
           testtrials.addData('roiTest.timesOn', roiTest.timesOn)
           testtrials.addData('roiTest.timesOff', roiTest.timesOff)
        else:
           testtrials.addData('roiTest.timesOn', "")
           testtrials.addData('roiTest.timesOff', "")
        if portTestImage.status == STARTED:
            win.callOnFlip(portTestImage.setData, int(0))
        # using non-slip timing so subtract the expected duration of this Routine (unless ended on request)
        if routineForceEnded:
            routineTimer.reset()
        else:
            routineTimer.addTime(-10.000000)
        
        # --- Prepare to start Routine "crossEnd" ---
        continueRoutine = True
        routineForceEnded = False
        # update component parameters for each repeat
        # Run 'Begin Routine' code from code_end
        eyetracker.setRecordingState(False)
        # keep track of which components have finished
        crossEndComponents = [fixcrossE]
        for thisComponent in crossEndComponents:
            thisComponent.tStart = None
            thisComponent.tStop = None
            thisComponent.tStartRefresh = None
            thisComponent.tStopRefresh = None
            if hasattr(thisComponent, 'status'):
                thisComponent.status = NOT_STARTED
        # reset timers
        t = 0
        _timeToFirstFrame = win.getFutureFlipTime(clock="now")
        frameN = -1
        
        # --- Run Routine "crossEnd" ---
        while continueRoutine:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # *fixcrossE* updates
            if fixcrossE.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                fixcrossE.frameNStart = frameN  # exact frame index
                fixcrossE.tStart = t  # local t and not account for scr refresh
                fixcrossE.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(fixcrossE, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'fixcrossE.started')
                fixcrossE.setAutoDraw(True)
            if fixcrossE.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > fixcrossE.tStartRefresh +  3 + random()-frameTolerance:
                    # keep track of stop time/frame for later
                    fixcrossE.tStop = t  # not accounting for scr refresh
                    fixcrossE.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'fixcrossE.stopped')
                    fixcrossE.setAutoDraw(False)
            
            # check for quit (typically the Esc key)
            if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
                core.quit()
            
            # check if all components have finished
            if not continueRoutine:  # a component has requested a forced-end of Routine
                routineForceEnded = True
                break
            continueRoutine = False  # will revert to True if at least one component still running
            for thisComponent in crossEndComponents:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # --- Ending Routine "crossEnd" ---
        for thisComponent in crossEndComponents:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        # the Routine "crossEnd" was not non-slip safe, so reset the non-slip timer
        routineTimer.reset()
        thisExp.nextEntry()
        
    # completed 10.0 repeats of 'testtrials'
    
    
    # set up handler to look after randomisation of conditions etc
    ratingtrials3 = data.TrialHandler(nReps=1.0, method='random', 
        extraInfo=expInfo, originPath=-1,
        trialList=data.importConditions(stimFile),
        seed=None, name='ratingtrials3')
    thisExp.addLoop(ratingtrials3)  # add the loop to the experiment
    thisRatingtrials3 = ratingtrials3.trialList[0]  # so we can initialise stimuli with some values
    # abbreviate parameter names if possible (e.g. rgb = thisRatingtrials3.rgb)
    if thisRatingtrials3 != None:
        for paramName in thisRatingtrials3:
            exec('{} = thisRatingtrials3[paramName]'.format(paramName))
    
    for thisRatingtrials3 in ratingtrials3:
        currentLoop = ratingtrials3
        # abbreviate parameter names if possible (e.g. rgb = thisRatingtrials3.rgb)
        if thisRatingtrials3 != None:
            for paramName in thisRatingtrials3:
                exec('{} = thisRatingtrials3[paramName]'.format(paramName))
        
        # --- Prepare to start Routine "stimRating" ---
        continueRoutine = True
        routineForceEnded = False
        # update component parameters for each repeat
        # Run 'Begin Routine' code from codeStimRating
        logging.log(level=logging.INFO, msg=f'StimRating_{stimtype}')
        print("Rating of Stimulus: %s"%(stimtype))
        win.mouseVisible = True
        imageRating.setPos((0, 0.2))
        imageRating.setImage(eval(stimtype))
        sliderStim.reset()
        spaceStim.keys = []
        spaceStim.rt = []
        _spaceStim_allKeys = []
        # keep track of which components have finished
        stimRatingComponents = [textRateStim, imageRating, textUnpleasant, textPleasant, sliderStim, textSpaceStim, spaceStim]
        for thisComponent in stimRatingComponents:
            thisComponent.tStart = None
            thisComponent.tStop = None
            thisComponent.tStartRefresh = None
            thisComponent.tStopRefresh = None
            if hasattr(thisComponent, 'status'):
                thisComponent.status = NOT_STARTED
        # reset timers
        t = 0
        _timeToFirstFrame = win.getFutureFlipTime(clock="now")
        frameN = -1
        
        # --- Run Routine "stimRating" ---
        while continueRoutine:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # *textRateStim* updates
            if textRateStim.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                textRateStim.frameNStart = frameN  # exact frame index
                textRateStim.tStart = t  # local t and not account for scr refresh
                textRateStim.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(textRateStim, 'tStartRefresh')  # time at next scr refresh
                textRateStim.setAutoDraw(True)
            
            # *imageRating* updates
            if imageRating.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                # keep track of start time/frame for later
                imageRating.frameNStart = frameN  # exact frame index
                imageRating.tStart = t  # local t and not account for scr refresh
                imageRating.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(imageRating, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'imageRating.started')
                imageRating.setAutoDraw(True)
            
            # *textUnpleasant* updates
            if textUnpleasant.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                textUnpleasant.frameNStart = frameN  # exact frame index
                textUnpleasant.tStart = t  # local t and not account for scr refresh
                textUnpleasant.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(textUnpleasant, 'tStartRefresh')  # time at next scr refresh
                textUnpleasant.setAutoDraw(True)
            
            # *textPleasant* updates
            if textPleasant.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                textPleasant.frameNStart = frameN  # exact frame index
                textPleasant.tStart = t  # local t and not account for scr refresh
                textPleasant.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(textPleasant, 'tStartRefresh')  # time at next scr refresh
                textPleasant.setAutoDraw(True)
            
            # *sliderStim* updates
            if sliderStim.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                sliderStim.frameNStart = frameN  # exact frame index
                sliderStim.tStart = t  # local t and not account for scr refresh
                sliderStim.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sliderStim, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'sliderStim.started')
                sliderStim.setAutoDraw(True)
            
            # *textSpaceStim* updates
            if textSpaceStim.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                textSpaceStim.frameNStart = frameN  # exact frame index
                textSpaceStim.tStart = t  # local t and not account for scr refresh
                textSpaceStim.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(textSpaceStim, 'tStartRefresh')  # time at next scr refresh
                textSpaceStim.setAutoDraw(True)
            
            # *spaceStim* updates
            waitOnFlip = False
            if spaceStim.status == NOT_STARTED and sliderStim.rating:
                # keep track of start time/frame for later
                spaceStim.frameNStart = frameN  # exact frame index
                spaceStim.tStart = t  # local t and not account for scr refresh
                spaceStim.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(spaceStim, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'spaceStim.started')
                spaceStim.status = STARTED
                # keyboard checking is just starting
                waitOnFlip = True
                win.callOnFlip(spaceStim.clock.reset)  # t=0 on next screen flip
                win.callOnFlip(spaceStim.clearEvents, eventType='keyboard')  # clear events on next screen flip
            if spaceStim.status == STARTED and not waitOnFlip:
                theseKeys = spaceStim.getKeys(keyList=['space', 'enter'], waitRelease=False)
                _spaceStim_allKeys.extend(theseKeys)
                if len(_spaceStim_allKeys):
                    spaceStim.keys = _spaceStim_allKeys[-1].name  # just the last key pressed
                    spaceStim.rt = _spaceStim_allKeys[-1].rt
                    # a response ends the routine
                    continueRoutine = False
            
            # check for quit (typically the Esc key)
            if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
                core.quit()
            
            # check if all components have finished
            if not continueRoutine:  # a component has requested a forced-end of Routine
                routineForceEnded = True
                break
            continueRoutine = False  # will revert to True if at least one component still running
            for thisComponent in stimRatingComponents:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # --- Ending Routine "stimRating" ---
        for thisComponent in stimRatingComponents:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        ratingtrials3.addData('sliderStim.response', sliderStim.getRating())
        ratingtrials3.addData('sliderStim.rt', sliderStim.getRT())
        # check responses
        if spaceStim.keys in ['', [], None]:  # No response was made
            spaceStim.keys = None
        ratingtrials3.addData('spaceStim.keys',spaceStim.keys)
        if spaceStim.keys != None:  # we had a response
            ratingtrials3.addData('spaceStim.rt', spaceStim.rt)
        # the Routine "stimRating" was not non-slip safe, so reset the non-slip timer
        routineTimer.reset()
        thisExp.nextEntry()
        
    # completed 1.0 repeats of 'ratingtrials3'
    
# completed 1.0 repeats of 'blocks'


# --- Prepare to start Routine "end" ---
continueRoutine = True
routineForceEnded = False
# update component parameters for each repeat
# keep track of which components have finished
endComponents = [textEnd]
for thisComponent in endComponents:
    thisComponent.tStart = None
    thisComponent.tStop = None
    thisComponent.tStartRefresh = None
    thisComponent.tStopRefresh = None
    if hasattr(thisComponent, 'status'):
        thisComponent.status = NOT_STARTED
# reset timers
t = 0
_timeToFirstFrame = win.getFutureFlipTime(clock="now")
frameN = -1

# --- Run Routine "end" ---
while continueRoutine and routineTimer.getTime() < 5.0:
    # get current time
    t = routineTimer.getTime()
    tThisFlip = win.getFutureFlipTime(clock=routineTimer)
    tThisFlipGlobal = win.getFutureFlipTime(clock=None)
    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
    # update/draw components on each frame
    
    # *textEnd* updates
    if textEnd.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
        # keep track of start time/frame for later
        textEnd.frameNStart = frameN  # exact frame index
        textEnd.tStart = t  # local t and not account for scr refresh
        textEnd.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(textEnd, 'tStartRefresh')  # time at next scr refresh
        textEnd.setAutoDraw(True)
    if textEnd.status == STARTED:
        # is it time to stop? (based on global clock, using actual start)
        if tThisFlipGlobal > textEnd.tStartRefresh + 5-frameTolerance:
            # keep track of stop time/frame for later
            textEnd.tStop = t  # not accounting for scr refresh
            textEnd.frameNStop = frameN  # exact frame index
            textEnd.setAutoDraw(False)
    
    # check for quit (typically the Esc key)
    if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
        core.quit()
    
    # check if all components have finished
    if not continueRoutine:  # a component has requested a forced-end of Routine
        routineForceEnded = True
        break
    continueRoutine = False  # will revert to True if at least one component still running
    for thisComponent in endComponents:
        if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
            continueRoutine = True
            break  # at least one component has not yet finished
    
    # refresh the screen
    if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
        win.flip()

# --- Ending Routine "end" ---
for thisComponent in endComponents:
    if hasattr(thisComponent, "setAutoDraw"):
        thisComponent.setAutoDraw(False)
# Run 'End Routine' code from codeEnd
print("----------- END -----------")
# using non-slip timing so subtract the expected duration of this Routine (unless ended on request)
if routineForceEnded:
    routineTimer.reset()
else:
    routineTimer.addTime(-5.000000)

# --- End experiment ---
# Flip one final time so any remaining win.callOnFlip() 
# and win.timeOnFlip() tasks get executed before quitting
win.flip()

# these shouldn't be strictly necessary (should auto-save)
thisExp.saveAsWideText(filename+'.csv', delim='auto')
thisExp.saveAsPickle(filename)
logging.flush()
# make sure everything is closed down
if eyetracker:
    eyetracker.setConnectionState(False)
thisExp.abort()  # or data files will save again on exit
win.close()
core.quit()
