#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This experiment was created using PsychoPy3 Experiment Builder (v2023.1.1),
    on Juli 12, 2023, at 08:14
If you publish work using this script the most relevant publication is:

    Peirce J, Gray JR, Simpson S, MacAskill M, Höchenberger R, Sogo H, Kastman E, Lindeløv JK. (2019) 
        PsychoPy2: Experiments in behavior made easy Behav Res 51: 195. 
        https://doi.org/10.3758/s13428-018-01193-y

"""

# --- Import packages ---
from psychopy import locale_setup
from psychopy import prefs
from psychopy import plugins
plugins.activatePlugins()
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
imagesize_test = [0.2, 0.2]
imagesize_rating = [0.2, 0.2]

score = 0
points = 0

looked_at = False
# Run 'Before Experiment' code from codeFeedback
soundFeedback = sound.Sound('A', secs=1, stereo=True, hamming=True, name='soundFeedback')


# Ensure that relative paths start from the same directory as this script
_thisDir = os.path.dirname(os.path.abspath(__file__))
os.chdir(_thisDir)
# Store info about the experiment session
psychopyVersion = '2023.1.1'
expName = 'gca'  # from the Builder filename that created this script
expInfo = {
    'participant': '',
    'session': '001',
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
filename = _thisDir + os.sep + u'data/%s_%s_%s' % (expInfo['participant'], expName, expInfo['date'])

# An ExperimentHandler isn't essential but helps with data saving
thisExp = data.ExperimentHandler(name=expName, version='',
    extraInfo=expInfo, runtimeInfo=None,
    originPath='C:\\Users\\Public\\Documents\\Projects\\GCA\\gaze_avoidance\\gca_lastrun.py',
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
    size=[2194, 1234], fullscr=True, screen=1, 
    winType='pyglet', allowStencil=False,
    monitor='labMonitor', color=[0,0,0], colorSpace='rgb',
    backgroundImage='', backgroundFit='none',
    blendMode='avg', useFBO=True, 
    units='height')
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
ioConfig['eyetracker.hw.mouse.EyeTracker'] = {
    'name': 'tracker',
    'controls': {
        'move': [],
        'blink':('MIDDLE_BUTTON',),
        'saccade_threshold': 0.5,
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
    text='Sehr geehrter Teilnehmer, sehr geehrte Teilnehmerin,\nvielen Dank für Ihre Teilnahme an unserem Experiment.\n\nBitte drücken Sie die Leertaste um zu starten.',
    font='Open Sans',
    pos=(0, 0), height=0.03, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=0.0);
spaceStart = keyboard.Keyboard()

# --- Initialize components for Routine "startPainThreshold" ---
textPainThreshold = visual.TextStim(win=win, name='textPainThreshold',
    text='Wir möchten zuerst den elektrischen Reiz richtig einstellen.\n\nBitte drücken Sie die Leertaste um zu starten.',
    font='Open Sans',
    pos=(0, 0), height=0.03, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=0.0);
spacePainThreshold = keyboard.Keyboard()

# --- Initialize components for Routine "painThreshold" ---
textRateT = visual.TextStim(win=win, name='textRateT',
    text='Wie unangenehm fanden Sie den elektrischen Reiz?',
    font='Open Sans',
    pos=(0, 0.3), height=0.05, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=0.0);
text_noPainT = visual.TextStim(win=win, name='text_noPainT',
    text='nichts\ngespürt',
    font='Open Sans',
    pos=(-0.5, 0.1), height=0.03, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-1.0);
text_slightPainT = visual.TextStim(win=win, name='text_slightPainT',
    text='eben \nwahrnehmbarer\nSchmerz',
    font='Open Sans',
    pos=(-0.1, 0.1), height=0.03, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-2.0);
text_highPainT = visual.TextStim(win=win, name='text_highPainT',
    text='unerträglicher\nSchmerz',
    font='Open Sans',
    pos=(0.5, 0.1), height=0.03, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-3.0);
sliderPainT = visual.Slider(win=win, name='sliderPainT',
    startValue=None, size=(1.0, 0.05), pos=(0, 0), units=win.units,
    labels=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10), ticks=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10), granularity=1.0,
    style='rating', styleTweaks=('labels45',), opacity=None,
    labelColor='White', markerColor='Red', lineColor='White', colorSpace='rgb',
    font='Open Sans', labelHeight=0.02,
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
    pos=(0, 0), height=0.05, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-4.0);

# --- Initialize components for Routine "startETCalibration" ---
text_ETCalibration = visual.TextStim(win=win, name='text_ETCalibration',
    text='Wir machen nun mit der Kalibrierung des Eye-Trackers weiter.\n\nBitte drücken Sie wieder die Leertaste um zu starten.',
    font='Open Sans',
    pos=(0, 0), height=0.03, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=0.0);
spaceETCalibration = keyboard.Keyboard()

# --- Initialize components for Routine "startExp" ---
textStartExp = visual.TextStim(win=win, name='textStartExp',
    text='Die Kalibrierung ist nun abgeschlossen.\n\nSie werden in unserem Experiment unterschiedliche Bilder in den vier Ecken des Bildschirms präsentiert bekommen.\nVersuchen Sie Ihren Score zu maximieren.\n\nBitte drücken Sie die Leertaste um zu starten.',
    font='Open Sans',
    pos=(0, 0), height=0.03, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=0.0);
spaceStartExp = keyboard.Keyboard()

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
    pos=(0, 0), height=0.05, wrapWidth=None, ori=0.0, 
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
    image='default.png', mask=None, anchor='center',
    ori=0.0, pos=[0,0], size=None,
    color=[1,1,1], colorSpace='rgb', opacity=None,
    flipHoriz=False, flipVert=False,
    texRes=128.0, interpolate=True, depth=-1.0)
roi = visual.ROI(win, name='roi', device=eyetracker,
    debug=False,
    shape='rectangle',
    pos=[0,0], size=[image.size], 
    anchor='center', ori=0.0, depth=-2
    )
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
    image='default.png', mask=None, anchor='center',
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
    pos=(0, 0), height=0.05, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-1.0);

# --- Initialize components for Routine "painRating" ---
textRate = visual.TextStim(win=win, name='textRate',
    text='Wie unangenehm fanden Sie den letzten elektrischen Reiz?',
    font='Open Sans',
    pos=(0, 0.3), height=0.05, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=0.0);
text_noPain = visual.TextStim(win=win, name='text_noPain',
    text='nichts\ngespürt',
    font='Open Sans',
    pos=(-0.5, 0.1), height=0.03, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-1.0);
text_slightPain = visual.TextStim(win=win, name='text_slightPain',
    text='eben \nwahrnehmbarer\nSchmerz',
    font='Open Sans',
    pos=(-0.1, 0.1), height=0.03, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-2.0);
text_highPain = visual.TextStim(win=win, name='text_highPain',
    text='unerträglicher\nSchmerz',
    font='Open Sans',
    pos=(0.5, 0.1), height=0.03, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-3.0);
sliderPain = visual.Slider(win=win, name='sliderPain',
    startValue=None, size=(1.0, 0.05), pos=(0, 0), units=win.units,
    labels=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10), ticks=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10), granularity=1.0,
    style='rating', styleTweaks=('labels45',), opacity=None,
    labelColor='White', markerColor='Red', lineColor='White', colorSpace='rgb',
    font='Open Sans', labelHeight=0.02,
    flip=False, ori=0.0, depth=-4, readOnly=False)
textSpacePain = visual.TextStim(win=win, name='textSpacePain',
    text='Bitte drücken Sie die Leertaste, um die Antwort zu bestätigen.',
    font='Open Sans',
    pos=(0, -0.3), height=0.05, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-5.0);
spacePain = keyboard.Keyboard()

# --- Initialize components for Routine "stimRating" ---
textRateStim = visual.TextStim(win=win, name='textRateStim',
    text='Wie wohl fühlen Sie sich bei der Betrachtung des Bildes?',
    font='Open Sans',
    pos=(0, 0.4), height=0.05, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=0.0);
imageRating = visual.ImageStim(
    win=win,
    name='imageRating', 
    image='default.png', mask=None, anchor='center',
    ori=0.0, pos=[0,0], size=1.0,
    color=[1,1,1], colorSpace='rgb', opacity=None,
    flipHoriz=False, flipVert=False,
    texRes=128.0, interpolate=True, depth=-1.0)
textUnpleasant = visual.TextStim(win=win, name='textUnpleasant',
    text='sehr \nunwohl',
    font='Open Sans',
    pos=(-0.5, -0.1), height=0.03, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-2.0);
textPleasant = visual.TextStim(win=win, name='textPleasant',
    text='sehr\nwohl',
    font='Open Sans',
    pos=(0.5, -0.1), height=0.03, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-3.0);
sliderStim = visual.Slider(win=win, name='sliderStim',
    startValue=None, size=(1.0, 0.05), pos=(0, -0.2), units=win.units,
    labels=(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), ticks=(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), granularity=1.0,
    style='rating', styleTweaks=(), opacity=None,
    labelColor='White', markerColor='Red', lineColor='White', colorSpace='rgb',
    font='Open Sans', labelHeight=0.02,
    flip=False, ori=0.0, depth=-4, readOnly=False)
textSpaceStim = visual.TextStim(win=win, name='textSpaceStim',
    text='Bitte drücken Sie die Leertaste, um die Antwort zu bestätigen.',
    font='Open Sans',
    pos=(0, -0.4), height=0.05, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-5.0);
spaceStim = keyboard.Keyboard()

# --- Initialize components for Routine "testInstr" ---
textTestInstr = visual.TextStim(win=win, name='textTestInstr',
    text='Im nächsten Block sehen Sie noch einmal die unterschiedlichen Bilder.\n\nBitte drücken Sie die Leertaste um zu starten.',
    font='Open Sans',
    pos=(0, 0), height=0.03, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=0.0);
spaceTestInstr = keyboard.Keyboard()

# --- Initialize components for Routine "cross" ---
fixcross = visual.TextStim(win=win, name='fixcross',
    text='+',
    font='Open Sans',
    pos=(0, 0), height=0.05, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-1.0);

# --- Initialize components for Routine "testtrial" ---
imageTest = visual.ImageStim(
    win=win,
    name='imageTest', 
    image='default.png', mask=None, anchor='center',
    ori=0.0, pos=[0,0], size=1.0,
    color=[1,1,1], colorSpace='rgb', opacity=None,
    flipHoriz=False, flipVert=False,
    texRes=128.0, interpolate=True, depth=0.0)
portTestImage = parallel.ParallelPort(address='0x0378')

# --- Initialize components for Routine "crossEnd" ---
fixcrossE = visual.TextStim(win=win, name='fixcrossE',
    text='+',
    font='Open Sans',
    pos=(0, 0), height=0.05, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-1.0);

# --- Initialize components for Routine "stimRating" ---
textRateStim = visual.TextStim(win=win, name='textRateStim',
    text='Wie wohl fühlen Sie sich bei der Betrachtung des Bildes?',
    font='Open Sans',
    pos=(0, 0.4), height=0.05, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=0.0);
imageRating = visual.ImageStim(
    win=win,
    name='imageRating', 
    image='default.png', mask=None, anchor='center',
    ori=0.0, pos=[0,0], size=1.0,
    color=[1,1,1], colorSpace='rgb', opacity=None,
    flipHoriz=False, flipVert=False,
    texRes=128.0, interpolate=True, depth=-1.0)
textUnpleasant = visual.TextStim(win=win, name='textUnpleasant',
    text='sehr \nunwohl',
    font='Open Sans',
    pos=(-0.5, -0.1), height=0.03, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-2.0);
textPleasant = visual.TextStim(win=win, name='textPleasant',
    text='sehr\nwohl',
    font='Open Sans',
    pos=(0.5, -0.1), height=0.03, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-3.0);
sliderStim = visual.Slider(win=win, name='sliderStim',
    startValue=None, size=(1.0, 0.05), pos=(0, -0.2), units=win.units,
    labels=(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), ticks=(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), granularity=1.0,
    style='rating', styleTweaks=(), opacity=None,
    labelColor='White', markerColor='Red', lineColor='White', colorSpace='rgb',
    font='Open Sans', labelHeight=0.02,
    flip=False, ori=0.0, depth=-4, readOnly=False)
textSpaceStim = visual.TextStim(win=win, name='textSpaceStim',
    text='Bitte drücken Sie die Leertaste, um die Antwort zu bestätigen.',
    font='Open Sans',
    pos=(0, -0.4), height=0.05, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-5.0);
spaceStim = keyboard.Keyboard()

# --- Initialize components for Routine "end" ---
textEnd = visual.TextStim(win=win, name='textEnd',
    text='Sehr geehrter Teilnehmer, sehr geehrte Teilnehmerin,\nvielen Dank für Ihre Teilnahme an unserem Experiment.\n\nDas Experiment ist beendet. \nSie können der Versuchsleitung Bescheid sagen.',
    font='Open Sans',
    pos=(0, 0), height=0.03, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=0.0);

# Create some handy timers
globalClock = core.Clock()  # to track the time since experiment started
routineTimer = core.Clock()  # to track time remaining of each (possibly non-slip) routine 

# --- Prepare to start Routine "welcome" ---
continueRoutine = True
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
routineForceEnded = not continueRoutine
while continueRoutine:
    # get current time
    t = routineTimer.getTime()
    tThisFlip = win.getFutureFlipTime(clock=routineTimer)
    tThisFlipGlobal = win.getFutureFlipTime(clock=None)
    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
    # update/draw components on each frame
    
    # *textStart* updates
    
    # if textStart is starting this frame...
    if textStart.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        # keep track of start time/frame for later
        textStart.frameNStart = frameN  # exact frame index
        textStart.tStart = t  # local t and not account for scr refresh
        textStart.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(textStart, 'tStartRefresh')  # time at next scr refresh
        # add timestamp to datafile
        thisExp.timestampOnFlip(win, 'textStart.started')
        # update status
        textStart.status = STARTED
        textStart.setAutoDraw(True)
    
    # if textStart is active this frame...
    if textStart.status == STARTED:
        # update params
        pass
    
    # *spaceStart* updates
    waitOnFlip = False
    
    # if spaceStart is starting this frame...
    if spaceStart.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        # keep track of start time/frame for later
        spaceStart.frameNStart = frameN  # exact frame index
        spaceStart.tStart = t  # local t and not account for scr refresh
        spaceStart.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(spaceStart, 'tStartRefresh')  # time at next scr refresh
        # add timestamp to datafile
        thisExp.timestampOnFlip(win, 'spaceStart.started')
        # update status
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
        if eyetracker:
            eyetracker.setConnectionState(False)
    
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
routineForceEnded = not continueRoutine
while continueRoutine:
    # get current time
    t = routineTimer.getTime()
    tThisFlip = win.getFutureFlipTime(clock=routineTimer)
    tThisFlipGlobal = win.getFutureFlipTime(clock=None)
    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
    # update/draw components on each frame
    
    # *textPainThreshold* updates
    
    # if textPainThreshold is starting this frame...
    if textPainThreshold.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        # keep track of start time/frame for later
        textPainThreshold.frameNStart = frameN  # exact frame index
        textPainThreshold.tStart = t  # local t and not account for scr refresh
        textPainThreshold.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(textPainThreshold, 'tStartRefresh')  # time at next scr refresh
        # add timestamp to datafile
        thisExp.timestampOnFlip(win, 'textPainThreshold.started')
        # update status
        textPainThreshold.status = STARTED
        textPainThreshold.setAutoDraw(True)
    
    # if textPainThreshold is active this frame...
    if textPainThreshold.status == STARTED:
        # update params
        pass
    
    # *spacePainThreshold* updates
    waitOnFlip = False
    
    # if spacePainThreshold is starting this frame...
    if spacePainThreshold.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        # keep track of start time/frame for later
        spacePainThreshold.frameNStart = frameN  # exact frame index
        spacePainThreshold.tStart = t  # local t and not account for scr refresh
        spacePainThreshold.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(spacePainThreshold, 'tStartRefresh')  # time at next scr refresh
        # add timestamp to datafile
        thisExp.timestampOnFlip(win, 'spacePainThreshold.started')
        # update status
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
        if eyetracker:
            eyetracker.setConnectionState(False)
    
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
    # update component parameters for each repeat
    sliderPainT.reset()
    keys_PainT.keys = []
    keys_PainT.rt = []
    _keys_PainT_allKeys = []
    # Run 'Begin Routine' code from codePainRatingT
    logging.log(level=logging.INFO, msg=f'PainThreshold')
    print("Press 's' to administer shock. Press 'g' to continue when pain threshold was found.")
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
    routineForceEnded = not continueRoutine
    while continueRoutine:
        # get current time
        t = routineTimer.getTime()
        tThisFlip = win.getFutureFlipTime(clock=routineTimer)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        
        # *textRateT* updates
        
        # if textRateT is starting this frame...
        if textRateT.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            textRateT.frameNStart = frameN  # exact frame index
            textRateT.tStart = t  # local t and not account for scr refresh
            textRateT.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(textRateT, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'textRateT.started')
            # update status
            textRateT.status = STARTED
            textRateT.setAutoDraw(True)
        
        # if textRateT is active this frame...
        if textRateT.status == STARTED:
            # update params
            pass
        
        # *text_noPainT* updates
        
        # if text_noPainT is starting this frame...
        if text_noPainT.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            text_noPainT.frameNStart = frameN  # exact frame index
            text_noPainT.tStart = t  # local t and not account for scr refresh
            text_noPainT.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(text_noPainT, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'text_noPainT.started')
            # update status
            text_noPainT.status = STARTED
            text_noPainT.setAutoDraw(True)
        
        # if text_noPainT is active this frame...
        if text_noPainT.status == STARTED:
            # update params
            pass
        
        # *text_slightPainT* updates
        
        # if text_slightPainT is starting this frame...
        if text_slightPainT.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            text_slightPainT.frameNStart = frameN  # exact frame index
            text_slightPainT.tStart = t  # local t and not account for scr refresh
            text_slightPainT.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(text_slightPainT, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'text_slightPainT.started')
            # update status
            text_slightPainT.status = STARTED
            text_slightPainT.setAutoDraw(True)
        
        # if text_slightPainT is active this frame...
        if text_slightPainT.status == STARTED:
            # update params
            pass
        
        # *text_highPainT* updates
        
        # if text_highPainT is starting this frame...
        if text_highPainT.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            text_highPainT.frameNStart = frameN  # exact frame index
            text_highPainT.tStart = t  # local t and not account for scr refresh
            text_highPainT.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(text_highPainT, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'text_highPainT.started')
            # update status
            text_highPainT.status = STARTED
            text_highPainT.setAutoDraw(True)
        
        # if text_highPainT is active this frame...
        if text_highPainT.status == STARTED:
            # update params
            pass
        
        # *sliderPainT* updates
        
        # if sliderPainT is starting this frame...
        if sliderPainT.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            sliderPainT.frameNStart = frameN  # exact frame index
            sliderPainT.tStart = t  # local t and not account for scr refresh
            sliderPainT.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(sliderPainT, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'sliderPainT.started')
            # update status
            sliderPainT.status = STARTED
            sliderPainT.setAutoDraw(True)
        
        # if sliderPainT is active this frame...
        if sliderPainT.status == STARTED:
            # update params
            pass
        
        # *keys_PainT* updates
        waitOnFlip = False
        
        # if keys_PainT is starting this frame...
        if keys_PainT.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
            # keep track of start time/frame for later
            keys_PainT.frameNStart = frameN  # exact frame index
            keys_PainT.tStart = t  # local t and not account for scr refresh
            keys_PainT.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(keys_PainT, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'keys_PainT.started')
            # update status
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
        
        # check for quit (typically the Esc key)
        if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
            core.quit()
            if eyetracker:
                eyetracker.setConnectionState(False)
        
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
    routineForceEnded = not continueRoutine
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
        
        # if portShock1T is starting this frame...
        if portShock1T.status == NOT_STARTED and t >= 0.1-frameTolerance:
            # keep track of start time/frame for later
            portShock1T.frameNStart = frameN  # exact frame index
            portShock1T.tStart = t  # local t and not account for scr refresh
            portShock1T.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(portShock1T, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.addData('portShock1T.started', t)
            # update status
            portShock1T.status = STARTED
            portShock1T.status = STARTED
            win.callOnFlip(portShock1T.setData, int(132))
        
        # if portShock1T is stopping this frame...
        if portShock1T.status == STARTED:
            # is it time to stop? (based on global clock, using actual start)
            if tThisFlipGlobal > portShock1T.tStartRefresh + 0.005-frameTolerance:
                # keep track of stop time/frame for later
                portShock1T.tStop = t  # not accounting for scr refresh
                portShock1T.frameNStop = frameN  # exact frame index
                # add timestamp to datafile
                thisExp.addData('portShock1T.stopped', t)
                # update status
                portShock1T.status = FINISHED
                win.callOnFlip(portShock1T.setData, int(0))
        # *portShock2T* updates
        
        # if portShock2T is starting this frame...
        if portShock2T.status == NOT_STARTED and t >= 0.15-frameTolerance:
            # keep track of start time/frame for later
            portShock2T.frameNStart = frameN  # exact frame index
            portShock2T.tStart = t  # local t and not account for scr refresh
            portShock2T.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(portShock2T, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.addData('portShock2T.started', t)
            # update status
            portShock2T.status = STARTED
            portShock2T.status = STARTED
            win.callOnFlip(portShock2T.setData, int(132))
        
        # if portShock2T is stopping this frame...
        if portShock2T.status == STARTED:
            # is it time to stop? (based on global clock, using actual start)
            if tThisFlipGlobal > portShock2T.tStartRefresh + 0.005-frameTolerance:
                # keep track of stop time/frame for later
                portShock2T.tStop = t  # not accounting for scr refresh
                portShock2T.frameNStop = frameN  # exact frame index
                # add timestamp to datafile
                thisExp.addData('portShock2T.stopped', t)
                # update status
                portShock2T.status = FINISHED
                win.callOnFlip(portShock2T.setData, int(0))
        # *portShock3T* updates
        
        # if portShock3T is starting this frame...
        if portShock3T.status == NOT_STARTED and t >= 0.2-frameTolerance:
            # keep track of start time/frame for later
            portShock3T.frameNStart = frameN  # exact frame index
            portShock3T.tStart = t  # local t and not account for scr refresh
            portShock3T.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(portShock3T, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.addData('portShock3T.started', t)
            # update status
            portShock3T.status = STARTED
            portShock3T.status = STARTED
            win.callOnFlip(portShock3T.setData, int(132))
        
        # if portShock3T is stopping this frame...
        if portShock3T.status == STARTED:
            # is it time to stop? (based on global clock, using actual start)
            if tThisFlipGlobal > portShock3T.tStartRefresh + 0.005-frameTolerance:
                # keep track of stop time/frame for later
                portShock3T.tStop = t  # not accounting for scr refresh
                portShock3T.frameNStop = frameN  # exact frame index
                # add timestamp to datafile
                thisExp.addData('portShock3T.stopped', t)
                # update status
                portShock3T.status = FINISHED
                win.callOnFlip(portShock3T.setData, int(0))
        
        # *fixcrossT* updates
        
        # if fixcrossT is starting this frame...
        if fixcrossT.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            fixcrossT.frameNStart = frameN  # exact frame index
            fixcrossT.tStart = t  # local t and not account for scr refresh
            fixcrossT.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(fixcrossT, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'fixcrossT.started')
            # update status
            fixcrossT.status = STARTED
            fixcrossT.setAutoDraw(True)
        
        # if fixcrossT is active this frame...
        if fixcrossT.status == STARTED:
            # update params
            pass
        
        # if fixcrossT is stopping this frame...
        if fixcrossT.status == STARTED:
            # is it time to stop? (based on global clock, using actual start)
            if tThisFlipGlobal > fixcrossT.tStartRefresh + 0.5-frameTolerance:
                # keep track of stop time/frame for later
                fixcrossT.tStop = t  # not accounting for scr refresh
                fixcrossT.frameNStop = frameN  # exact frame index
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'fixcrossT.stopped')
                # update status
                fixcrossT.status = FINISHED
                fixcrossT.setAutoDraw(False)
        
        # check for quit (typically the Esc key)
        if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
            core.quit()
            if eyetracker:
                eyetracker.setConnectionState(False)
        
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
routineForceEnded = not continueRoutine
while continueRoutine:
    # get current time
    t = routineTimer.getTime()
    tThisFlip = win.getFutureFlipTime(clock=routineTimer)
    tThisFlipGlobal = win.getFutureFlipTime(clock=None)
    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
    # update/draw components on each frame
    
    # *text_ETCalibration* updates
    
    # if text_ETCalibration is starting this frame...
    if text_ETCalibration.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        # keep track of start time/frame for later
        text_ETCalibration.frameNStart = frameN  # exact frame index
        text_ETCalibration.tStart = t  # local t and not account for scr refresh
        text_ETCalibration.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(text_ETCalibration, 'tStartRefresh')  # time at next scr refresh
        # add timestamp to datafile
        thisExp.timestampOnFlip(win, 'text_ETCalibration.started')
        # update status
        text_ETCalibration.status = STARTED
        text_ETCalibration.setAutoDraw(True)
    
    # if text_ETCalibration is active this frame...
    if text_ETCalibration.status == STARTED:
        # update params
        pass
    
    # *spaceETCalibration* updates
    waitOnFlip = False
    
    # if spaceETCalibration is starting this frame...
    if spaceETCalibration.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        # keep track of start time/frame for later
        spaceETCalibration.frameNStart = frameN  # exact frame index
        spaceETCalibration.tStart = t  # local t and not account for scr refresh
        spaceETCalibration.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(spaceETCalibration, 'tStartRefresh')  # time at next scr refresh
        # add timestamp to datafile
        thisExp.timestampOnFlip(win, 'spaceETCalibration.started')
        # update status
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
        if eyetracker:
            eyetracker.setConnectionState(False)
    
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

# --- Prepare to start Routine "startExp" ---
continueRoutine = True
# update component parameters for each repeat
spaceStartExp.keys = []
spaceStartExp.rt = []
_spaceStartExp_allKeys = []
# Run 'Begin Routine' code from codeStartExperiment
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
routineForceEnded = not continueRoutine
while continueRoutine:
    # get current time
    t = routineTimer.getTime()
    tThisFlip = win.getFutureFlipTime(clock=routineTimer)
    tThisFlipGlobal = win.getFutureFlipTime(clock=None)
    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
    # update/draw components on each frame
    
    # *textStartExp* updates
    
    # if textStartExp is starting this frame...
    if textStartExp.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        # keep track of start time/frame for later
        textStartExp.frameNStart = frameN  # exact frame index
        textStartExp.tStart = t  # local t and not account for scr refresh
        textStartExp.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(textStartExp, 'tStartRefresh')  # time at next scr refresh
        # add timestamp to datafile
        thisExp.timestampOnFlip(win, 'textStartExp.started')
        # update status
        textStartExp.status = STARTED
        textStartExp.setAutoDraw(True)
    
    # if textStartExp is active this frame...
    if textStartExp.status == STARTED:
        # update params
        pass
    
    # *spaceStartExp* updates
    waitOnFlip = False
    
    # if spaceStartExp is starting this frame...
    if spaceStartExp.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        # keep track of start time/frame for later
        spaceStartExp.frameNStart = frameN  # exact frame index
        spaceStartExp.tStart = t  # local t and not account for scr refresh
        spaceStartExp.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(spaceStartExp, 'tStartRefresh')  # time at next scr refresh
        # add timestamp to datafile
        thisExp.timestampOnFlip(win, 'spaceStartExp.started')
        # update status
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
        if eyetracker:
            eyetracker.setConnectionState(False)
    
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
thisExp.addData('spaceStartExp.keys',spaceStartExp.keys)
if spaceStartExp.keys != None:  # we had a response
    thisExp.addData('spaceStartExp.rt', spaceStartExp.rt)
thisExp.nextEntry()
# the Routine "startExp" was not non-slip safe, so reset the non-slip timer
routineTimer.reset()

# --- Prepare to start Routine "blank" ---
continueRoutine = True
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
routineForceEnded = not continueRoutine
while continueRoutine:
    # get current time
    t = routineTimer.getTime()
    tThisFlip = win.getFutureFlipTime(clock=routineTimer)
    tThisFlipGlobal = win.getFutureFlipTime(clock=None)
    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
    # update/draw components on each frame
    
    # *blankScreen* updates
    
    # if blankScreen is starting this frame...
    if blankScreen.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        # keep track of start time/frame for later
        blankScreen.frameNStart = frameN  # exact frame index
        blankScreen.tStart = t  # local t and not account for scr refresh
        blankScreen.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(blankScreen, 'tStartRefresh')  # time at next scr refresh
        # add timestamp to datafile
        thisExp.timestampOnFlip(win, 'blankScreen.started')
        # update status
        blankScreen.status = STARTED
        blankScreen.setAutoDraw(True)
    
    # if blankScreen is active this frame...
    if blankScreen.status == STARTED:
        # update params
        pass
    
    # if blankScreen is stopping this frame...
    if blankScreen.status == STARTED:
        # is it time to stop? (based on global clock, using actual start)
        if tThisFlipGlobal > blankScreen.tStartRefresh + 1.0 + random()-frameTolerance:
            # keep track of stop time/frame for later
            blankScreen.tStop = t  # not accounting for scr refresh
            blankScreen.frameNStop = frameN  # exact frame index
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'blankScreen.stopped')
            # update status
            blankScreen.status = FINISHED
            blankScreen.setAutoDraw(False)
    
    # check for quit (typically the Esc key)
    if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
        core.quit()
        if eyetracker:
            eyetracker.setConnectionState(False)
    
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
        # update component parameters for each repeat
        # Run 'Begin Routine' code from codeFixate
        logging.log(level=logging.INFO, msg=f'FixationCross')
        eyetracker.setRecordingState(True)
        
        if currentLoop.name == "trials":
            eyetracker.sendMessage("VP %s TRIALID %d CONDITION %s"%(expInfo['participant'], trials.thisN+1, trialtype))
            print("VP %s TRIALID %d CONDITION %s"%(expInfo['participant'], trials.thisN+1, trialtype))
            
        if currentLoop.name == "testtrials":
            eyetracker.sendMessage("VP %s TESTTRIALID %d CONDITION %s"%(expInfo['participant'], testtrials.thisN+1, stimtype))
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
        routineForceEnded = not continueRoutine
        while continueRoutine and routineTimer.getTime() < 1.0:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # *fixcross* updates
            
            # if fixcross is starting this frame...
            if fixcross.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                fixcross.frameNStart = frameN  # exact frame index
                fixcross.tStart = t  # local t and not account for scr refresh
                fixcross.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(fixcross, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'fixcross.started')
                # update status
                fixcross.status = STARTED
                fixcross.setAutoDraw(True)
            
            # if fixcross is active this frame...
            if fixcross.status == STARTED:
                # update params
                pass
            
            # if fixcross is stopping this frame...
            if fixcross.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > fixcross.tStartRefresh + 1-frameTolerance:
                    # keep track of stop time/frame for later
                    fixcross.tStop = t  # not accounting for scr refresh
                    fixcross.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'fixcross.stopped')
                    # update status
                    fixcross.status = FINISHED
                    fixcross.setAutoDraw(False)
            
            # check for quit (typically the Esc key)
            if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
                core.quit()
                if eyetracker:
                    eyetracker.setConnectionState(False)
            
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
        routineForceEnded = not continueRoutine
        while continueRoutine:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # *blankScreen* updates
            
            # if blankScreen is starting this frame...
            if blankScreen.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                blankScreen.frameNStart = frameN  # exact frame index
                blankScreen.tStart = t  # local t and not account for scr refresh
                blankScreen.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(blankScreen, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'blankScreen.started')
                # update status
                blankScreen.status = STARTED
                blankScreen.setAutoDraw(True)
            
            # if blankScreen is active this frame...
            if blankScreen.status == STARTED:
                # update params
                pass
            
            # if blankScreen is stopping this frame...
            if blankScreen.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > blankScreen.tStartRefresh + 1.0 + random()-frameTolerance:
                    # keep track of stop time/frame for later
                    blankScreen.tStop = t  # not accounting for scr refresh
                    blankScreen.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'blankScreen.stopped')
                    # update status
                    blankScreen.status = FINISHED
                    blankScreen.setAutoDraw(False)
            
            # check for quit (typically the Esc key)
            if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
                core.quit()
                if eyetracker:
                    eyetracker.setConnectionState(False)
            
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
        # update component parameters for each repeat
        # Run 'Begin Routine' code from codeTrial
        looked_at = False
        cursorcolor="white"
        
        logging.log(level=logging.INFO, msg=f'ImageOnset_{trialtype}')
        image.setPos(position)
        image.setImage(eval(trialtype))
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
        routineForceEnded = not continueRoutine
        while continueRoutine and routineTimer.getTime() < 5.0:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            # Run 'Each Frame' code from codeTrial
            if roi.isLookedIn:
                looked_at = True
                continueRoutine = False
            
            
            # *image* updates
            
            # if image is starting this frame...
            if image.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                # keep track of start time/frame for later
                image.frameNStart = frameN  # exact frame index
                image.tStart = t  # local t and not account for scr refresh
                image.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(image, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'image.started')
                # update status
                image.status = STARTED
                image.setAutoDraw(True)
            
            # if image is active this frame...
            if image.status == STARTED:
                # update params
                pass
            
            # if image is stopping this frame...
            if image.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > image.tStartRefresh + 5-frameTolerance:
                    # keep track of stop time/frame for later
                    image.tStop = t  # not accounting for scr refresh
                    image.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'image.stopped')
                    # update status
                    image.status = FINISHED
                    image.setAutoDraw(False)
            
            # if roi is starting this frame...
            if roi.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                roi.frameNStart = frameN  # exact frame index
                roi.tStart = t  # local t and not account for scr refresh
                roi.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(roi, 'tStartRefresh')  # time at next scr refresh
                # update status
                roi.status = STARTED
                roi.setAutoDraw(True)
            
            # if roi is active this frame...
            if roi.status == STARTED:
                # update params
                roi.setPos(position, log=False)
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
            
            # if roi is stopping this frame...
            if roi.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > roi.tStartRefresh + 5-frameTolerance:
                    # keep track of stop time/frame for later
                    roi.tStop = t  # not accounting for scr refresh
                    roi.frameNStop = frameN  # exact frame index
                    # update status
                    roi.status = FINISHED
                    roi.setAutoDraw(False)
            
            # *gazeCursor* updates
            
            # if gazeCursor is starting this frame...
            if gazeCursor.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                gazeCursor.frameNStart = frameN  # exact frame index
                gazeCursor.tStart = t  # local t and not account for scr refresh
                gazeCursor.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(gazeCursor, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'gazeCursor.started')
                # update status
                gazeCursor.status = STARTED
                gazeCursor.setAutoDraw(True)
            
            # if gazeCursor is active this frame...
            if gazeCursor.status == STARTED:
                # update params
                gazeCursor.setFillColor(cursorcolor, log=False)
                gazeCursor.setOpacity(1.0, log=False)
                gazeCursor.setPos([eyetracker.getPos()], log=False)
            
            # if gazeCursor is stopping this frame...
            if gazeCursor.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > gazeCursor.tStartRefresh + 5-frameTolerance:
                    # keep track of stop time/frame for later
                    gazeCursor.tStop = t  # not accounting for scr refresh
                    gazeCursor.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'gazeCursor.stopped')
                    # update status
                    gazeCursor.status = FINISHED
                    gazeCursor.setAutoDraw(False)
            # *portImage* updates
            
            # if portImage is starting this frame...
            if portImage.status == NOT_STARTED and t >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                portImage.frameNStart = frameN  # exact frame index
                portImage.tStart = t  # local t and not account for scr refresh
                portImage.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(portImage, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.addData('portImage.started', t)
                # update status
                portImage.status = STARTED
                portImage.status = STARTED
                win.callOnFlip(portImage.setData, int(1))
            
            # if portImage is stopping this frame...
            if portImage.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > portImage.tStartRefresh + 0.005-frameTolerance:
                    # keep track of stop time/frame for later
                    portImage.tStop = t  # not accounting for scr refresh
                    portImage.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.addData('portImage.stopped', t)
                    # update status
                    portImage.status = FINISHED
                    win.callOnFlip(portImage.setData, int(0))
            
            # check for quit (typically the Esc key)
            if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
                core.quit()
                if eyetracker:
                    eyetracker.setConnectionState(False)
            
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
            points = 5
            score -= points
            feedback_points = f"- {points}"
            feedback_score = f"Score: {score}"
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
        
        rectsize = [item * 1.05 for item in image.size]
        
        image_w = image.size[0]
        image_h = image.size[1]
        image_h_spec_rating = 0.3
        imagesize_rating = [image_w * image_h_spec_rating / image_h, image_h_spec_rating]
        
        image_h_spec_test = 0.5
        imagesize_test = [image_w * image_h_spec_test / image_h, image_h_spec_test]
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
        # update component parameters for each repeat
        # Run 'Begin Routine' code from codeFeedback
        logging.log(level=logging.INFO, msg=f'FeedbackOnset_{log_msg}')
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
        routineForceEnded = not continueRoutine
        while continueRoutine and routineTimer.getTime() < 3.0:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # *textPoints* updates
            
            # if textPoints is starting this frame...
            if textPoints.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                # keep track of start time/frame for later
                textPoints.frameNStart = frameN  # exact frame index
                textPoints.tStart = t  # local t and not account for scr refresh
                textPoints.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(textPoints, 'tStartRefresh')  # time at next scr refresh
                # update status
                textPoints.status = STARTED
                textPoints.setAutoDraw(True)
            
            # if textPoints is active this frame...
            if textPoints.status == STARTED:
                # update params
                textPoints.setText(feedback_points, log=False)
            
            # if textPoints is stopping this frame...
            if textPoints.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > textPoints.tStartRefresh + 3-frameTolerance:
                    # keep track of stop time/frame for later
                    textPoints.tStop = t  # not accounting for scr refresh
                    textPoints.frameNStop = frameN  # exact frame index
                    # update status
                    textPoints.status = FINISHED
                    textPoints.setAutoDraw(False)
            
            # *textScore* updates
            
            # if textScore is starting this frame...
            if textScore.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                # keep track of start time/frame for later
                textScore.frameNStart = frameN  # exact frame index
                textScore.tStart = t  # local t and not account for scr refresh
                textScore.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(textScore, 'tStartRefresh')  # time at next scr refresh
                # update status
                textScore.status = STARTED
                textScore.setAutoDraw(True)
            
            # if textScore is active this frame...
            if textScore.status == STARTED:
                # update params
                textScore.setText(feedback_score, log=False)
            
            # if textScore is stopping this frame...
            if textScore.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > textScore.tStartRefresh + 3-frameTolerance:
                    # keep track of stop time/frame for later
                    textScore.tStop = t  # not accounting for scr refresh
                    textScore.frameNStop = frameN  # exact frame index
                    # update status
                    textScore.status = FINISHED
                    textScore.setAutoDraw(False)
            
            # *polygon* updates
            
            # if polygon is starting this frame...
            if polygon.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                polygon.frameNStart = frameN  # exact frame index
                polygon.tStart = t  # local t and not account for scr refresh
                polygon.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(polygon, 'tStartRefresh')  # time at next scr refresh
                # update status
                polygon.status = STARTED
                polygon.setAutoDraw(True)
            
            # if polygon is active this frame...
            if polygon.status == STARTED:
                # update params
                polygon.setSize(rectsize, log=False)
            
            # if polygon is stopping this frame...
            if polygon.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > polygon.tStartRefresh + 3-frameTolerance:
                    # keep track of stop time/frame for later
                    polygon.tStop = t  # not accounting for scr refresh
                    polygon.frameNStop = frameN  # exact frame index
                    # update status
                    polygon.status = FINISHED
                    polygon.setAutoDraw(False)
            
            # *imageFeedback* updates
            
            # if imageFeedback is starting this frame...
            if imageFeedback.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                # keep track of start time/frame for later
                imageFeedback.frameNStart = frameN  # exact frame index
                imageFeedback.tStart = t  # local t and not account for scr refresh
                imageFeedback.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(imageFeedback, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'imageFeedback.started')
                # update status
                imageFeedback.status = STARTED
                imageFeedback.setAutoDraw(True)
            
            # if imageFeedback is active this frame...
            if imageFeedback.status == STARTED:
                # update params
                pass
            
            # if imageFeedback is stopping this frame...
            if imageFeedback.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > imageFeedback.tStartRefresh + 3-frameTolerance:
                    # keep track of stop time/frame for later
                    imageFeedback.tStop = t  # not accounting for scr refresh
                    imageFeedback.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'imageFeedback.stopped')
                    # update status
                    imageFeedback.status = FINISHED
                    imageFeedback.setAutoDraw(False)
            
            # *gazeCursor_Feedback* updates
            
            # if gazeCursor_Feedback is starting this frame...
            if gazeCursor_Feedback.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                gazeCursor_Feedback.frameNStart = frameN  # exact frame index
                gazeCursor_Feedback.tStart = t  # local t and not account for scr refresh
                gazeCursor_Feedback.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(gazeCursor_Feedback, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'gazeCursor_Feedback.started')
                # update status
                gazeCursor_Feedback.status = STARTED
                gazeCursor_Feedback.setAutoDraw(True)
            
            # if gazeCursor_Feedback is active this frame...
            if gazeCursor_Feedback.status == STARTED:
                # update params
                gazeCursor_Feedback.setFillColor(cursorcolor, log=False)
                gazeCursor_Feedback.setOpacity(0.0, log=False)
                gazeCursor_Feedback.setPos([eyetracker.getPos()], log=False)
            
            # if gazeCursor_Feedback is stopping this frame...
            if gazeCursor_Feedback.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > gazeCursor_Feedback.tStartRefresh + 3-frameTolerance:
                    # keep track of stop time/frame for later
                    gazeCursor_Feedback.tStop = t  # not accounting for scr refresh
                    gazeCursor_Feedback.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'gazeCursor_Feedback.stopped')
                    # update status
                    gazeCursor_Feedback.status = FINISHED
                    gazeCursor_Feedback.setAutoDraw(False)
            # start/stop soundFeedback
            
            # if soundFeedback is starting this frame...
            if soundFeedback.status == NOT_STARTED and t >= 0-frameTolerance:
                # keep track of start time/frame for later
                soundFeedback.frameNStart = frameN  # exact frame index
                soundFeedback.tStart = t  # local t and not account for scr refresh
                soundFeedback.tStartRefresh = tThisFlipGlobal  # on global time
                # add timestamp to datafile
                thisExp.addData('soundFeedback.started', t)
                # update status
                soundFeedback.status = STARTED
                soundFeedback.play()  # start the sound (it finishes automatically)
            
            # if soundFeedback is stopping this frame...
            if soundFeedback.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > soundFeedback.tStartRefresh + 1-frameTolerance:
                    # keep track of stop time/frame for later
                    soundFeedback.tStop = t  # not accounting for scr refresh
                    soundFeedback.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.addData('soundFeedback.stopped', t)
                    # update status
                    soundFeedback.status = FINISHED
                    soundFeedback.stop()
            # *portShock1* updates
            
            # if portShock1 is starting this frame...
            if portShock1.status == NOT_STARTED and t >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                portShock1.frameNStart = frameN  # exact frame index
                portShock1.tStart = t  # local t and not account for scr refresh
                portShock1.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(portShock1, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.addData('portShock1.started', t)
                # update status
                portShock1.status = STARTED
                portShock1.status = STARTED
                win.callOnFlip(portShock1.setData, int(port_msg))
            
            # if portShock1 is stopping this frame...
            if portShock1.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > portShock1.tStartRefresh + 0.005-frameTolerance:
                    # keep track of stop time/frame for later
                    portShock1.tStop = t  # not accounting for scr refresh
                    portShock1.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.addData('portShock1.stopped', t)
                    # update status
                    portShock1.status = FINISHED
                    win.callOnFlip(portShock1.setData, int(0))
            # *portShock2* updates
            
            # if portShock2 is starting this frame...
            if portShock2.status == NOT_STARTED and t >= 0.05-frameTolerance:
                # keep track of start time/frame for later
                portShock2.frameNStart = frameN  # exact frame index
                portShock2.tStart = t  # local t and not account for scr refresh
                portShock2.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(portShock2, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.addData('portShock2.started', t)
                # update status
                portShock2.status = STARTED
                portShock2.status = STARTED
                win.callOnFlip(portShock2.setData, int(128 if port_msg == 132 else 0))
            
            # if portShock2 is stopping this frame...
            if portShock2.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > portShock2.tStartRefresh + 0.005-frameTolerance:
                    # keep track of stop time/frame for later
                    portShock2.tStop = t  # not accounting for scr refresh
                    portShock2.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.addData('portShock2.stopped', t)
                    # update status
                    portShock2.status = FINISHED
                    win.callOnFlip(portShock2.setData, int(0))
            # *portShock3* updates
            
            # if portShock3 is starting this frame...
            if portShock3.status == NOT_STARTED and t >= 0.1-frameTolerance:
                # keep track of start time/frame for later
                portShock3.frameNStart = frameN  # exact frame index
                portShock3.tStart = t  # local t and not account for scr refresh
                portShock3.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(portShock3, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.addData('portShock3.started', t)
                # update status
                portShock3.status = STARTED
                portShock3.status = STARTED
                win.callOnFlip(portShock3.setData, int(128 if port_msg == 132 else 0))
            
            # if portShock3 is stopping this frame...
            if portShock3.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > portShock3.tStartRefresh + 0.005-frameTolerance:
                    # keep track of stop time/frame for later
                    portShock3.tStop = t  # not accounting for scr refresh
                    portShock3.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.addData('portShock3.stopped', t)
                    # update status
                    portShock3.status = FINISHED
                    win.callOnFlip(portShock3.setData, int(0))
            
            # check for quit (typically the Esc key)
            if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
                core.quit()
                if eyetracker:
                    eyetracker.setConnectionState(False)
            
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
        routineForceEnded = not continueRoutine
        while continueRoutine:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # *fixcrossE* updates
            
            # if fixcrossE is starting this frame...
            if fixcrossE.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                fixcrossE.frameNStart = frameN  # exact frame index
                fixcrossE.tStart = t  # local t and not account for scr refresh
                fixcrossE.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(fixcrossE, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'fixcrossE.started')
                # update status
                fixcrossE.status = STARTED
                fixcrossE.setAutoDraw(True)
            
            # if fixcrossE is active this frame...
            if fixcrossE.status == STARTED:
                # update params
                pass
            
            # if fixcrossE is stopping this frame...
            if fixcrossE.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > fixcrossE.tStartRefresh +  3 + random()-frameTolerance:
                    # keep track of stop time/frame for later
                    fixcrossE.tStop = t  # not accounting for scr refresh
                    fixcrossE.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'fixcrossE.stopped')
                    # update status
                    fixcrossE.status = FINISHED
                    fixcrossE.setAutoDraw(False)
            
            # check for quit (typically the Esc key)
            if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
                core.quit()
                if eyetracker:
                    eyetracker.setConnectionState(False)
            
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
    # update component parameters for each repeat
    sliderPain.reset()
    spacePain.keys = []
    spacePain.rt = []
    _spacePain_allKeys = []
    # Run 'Begin Routine' code from codePainRating
    logging.log(level=logging.INFO, msg=f'PainRating')
    print("Pain Rating. Press Space to continue.")
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
    routineForceEnded = not continueRoutine
    while continueRoutine:
        # get current time
        t = routineTimer.getTime()
        tThisFlip = win.getFutureFlipTime(clock=routineTimer)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        
        # *textRate* updates
        
        # if textRate is starting this frame...
        if textRate.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            textRate.frameNStart = frameN  # exact frame index
            textRate.tStart = t  # local t and not account for scr refresh
            textRate.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(textRate, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'textRate.started')
            # update status
            textRate.status = STARTED
            textRate.setAutoDraw(True)
        
        # if textRate is active this frame...
        if textRate.status == STARTED:
            # update params
            pass
        
        # *text_noPain* updates
        
        # if text_noPain is starting this frame...
        if text_noPain.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            text_noPain.frameNStart = frameN  # exact frame index
            text_noPain.tStart = t  # local t and not account for scr refresh
            text_noPain.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(text_noPain, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'text_noPain.started')
            # update status
            text_noPain.status = STARTED
            text_noPain.setAutoDraw(True)
        
        # if text_noPain is active this frame...
        if text_noPain.status == STARTED:
            # update params
            pass
        
        # *text_slightPain* updates
        
        # if text_slightPain is starting this frame...
        if text_slightPain.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            text_slightPain.frameNStart = frameN  # exact frame index
            text_slightPain.tStart = t  # local t and not account for scr refresh
            text_slightPain.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(text_slightPain, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'text_slightPain.started')
            # update status
            text_slightPain.status = STARTED
            text_slightPain.setAutoDraw(True)
        
        # if text_slightPain is active this frame...
        if text_slightPain.status == STARTED:
            # update params
            pass
        
        # *text_highPain* updates
        
        # if text_highPain is starting this frame...
        if text_highPain.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            text_highPain.frameNStart = frameN  # exact frame index
            text_highPain.tStart = t  # local t and not account for scr refresh
            text_highPain.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(text_highPain, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'text_highPain.started')
            # update status
            text_highPain.status = STARTED
            text_highPain.setAutoDraw(True)
        
        # if text_highPain is active this frame...
        if text_highPain.status == STARTED:
            # update params
            pass
        
        # *sliderPain* updates
        
        # if sliderPain is starting this frame...
        if sliderPain.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            sliderPain.frameNStart = frameN  # exact frame index
            sliderPain.tStart = t  # local t and not account for scr refresh
            sliderPain.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(sliderPain, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'sliderPain.started')
            # update status
            sliderPain.status = STARTED
            sliderPain.setAutoDraw(True)
        
        # if sliderPain is active this frame...
        if sliderPain.status == STARTED:
            # update params
            pass
        
        # *textSpacePain* updates
        
        # if textSpacePain is starting this frame...
        if textSpacePain.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            textSpacePain.frameNStart = frameN  # exact frame index
            textSpacePain.tStart = t  # local t and not account for scr refresh
            textSpacePain.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(textSpacePain, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'textSpacePain.started')
            # update status
            textSpacePain.status = STARTED
            textSpacePain.setAutoDraw(True)
        
        # if textSpacePain is active this frame...
        if textSpacePain.status == STARTED:
            # update params
            pass
        
        # *spacePain* updates
        waitOnFlip = False
        
        # if spacePain is starting this frame...
        if spacePain.status == NOT_STARTED and sliderPain.rating:
            # keep track of start time/frame for later
            spacePain.frameNStart = frameN  # exact frame index
            spacePain.tStart = t  # local t and not account for scr refresh
            spacePain.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(spacePain, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'spacePain.started')
            # update status
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
        
        # check for quit (typically the Esc key)
        if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
            core.quit()
            if eyetracker:
                eyetracker.setConnectionState(False)
        
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
    ratingstrials = data.TrialHandler(nReps=1.0, method='random', 
        extraInfo=expInfo, originPath=-1,
        trialList=data.importConditions(stimFile),
        seed=None, name='ratingstrials')
    thisExp.addLoop(ratingstrials)  # add the loop to the experiment
    thisRatingstrial = ratingstrials.trialList[0]  # so we can initialise stimuli with some values
    # abbreviate parameter names if possible (e.g. rgb = thisRatingstrial.rgb)
    if thisRatingstrial != None:
        for paramName in thisRatingstrial:
            exec('{} = thisRatingstrial[paramName]'.format(paramName))
    
    for thisRatingstrial in ratingstrials:
        currentLoop = ratingstrials
        # abbreviate parameter names if possible (e.g. rgb = thisRatingstrial.rgb)
        if thisRatingstrial != None:
            for paramName in thisRatingstrial:
                exec('{} = thisRatingstrial[paramName]'.format(paramName))
        
        # --- Prepare to start Routine "stimRating" ---
        continueRoutine = True
        # update component parameters for each repeat
        imageRating.setPos((0, 0.1))
        imageRating.setImage(eval(stimtype))
        sliderStim.reset()
        spaceStim.keys = []
        spaceStim.rt = []
        _spaceStim_allKeys = []
        # Run 'Begin Routine' code from codeStimRating
        logging.log(level=logging.INFO, msg=f'StimRating_{stimtype}')
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
        routineForceEnded = not continueRoutine
        while continueRoutine:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # *textRateStim* updates
            
            # if textRateStim is starting this frame...
            if textRateStim.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                textRateStim.frameNStart = frameN  # exact frame index
                textRateStim.tStart = t  # local t and not account for scr refresh
                textRateStim.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(textRateStim, 'tStartRefresh')  # time at next scr refresh
                # update status
                textRateStim.status = STARTED
                textRateStim.setAutoDraw(True)
            
            # if textRateStim is active this frame...
            if textRateStim.status == STARTED:
                # update params
                pass
            
            # *imageRating* updates
            
            # if imageRating is starting this frame...
            if imageRating.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                # keep track of start time/frame for later
                imageRating.frameNStart = frameN  # exact frame index
                imageRating.tStart = t  # local t and not account for scr refresh
                imageRating.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(imageRating, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'imageRating.started')
                # update status
                imageRating.status = STARTED
                imageRating.setAutoDraw(True)
            
            # if imageRating is active this frame...
            if imageRating.status == STARTED:
                # update params
                imageRating.setSize(imagesize_rating, log=False)
            
            # *textUnpleasant* updates
            
            # if textUnpleasant is starting this frame...
            if textUnpleasant.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                textUnpleasant.frameNStart = frameN  # exact frame index
                textUnpleasant.tStart = t  # local t and not account for scr refresh
                textUnpleasant.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(textUnpleasant, 'tStartRefresh')  # time at next scr refresh
                # update status
                textUnpleasant.status = STARTED
                textUnpleasant.setAutoDraw(True)
            
            # if textUnpleasant is active this frame...
            if textUnpleasant.status == STARTED:
                # update params
                pass
            
            # *textPleasant* updates
            
            # if textPleasant is starting this frame...
            if textPleasant.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                textPleasant.frameNStart = frameN  # exact frame index
                textPleasant.tStart = t  # local t and not account for scr refresh
                textPleasant.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(textPleasant, 'tStartRefresh')  # time at next scr refresh
                # update status
                textPleasant.status = STARTED
                textPleasant.setAutoDraw(True)
            
            # if textPleasant is active this frame...
            if textPleasant.status == STARTED:
                # update params
                pass
            
            # *sliderStim* updates
            
            # if sliderStim is starting this frame...
            if sliderStim.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                sliderStim.frameNStart = frameN  # exact frame index
                sliderStim.tStart = t  # local t and not account for scr refresh
                sliderStim.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sliderStim, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'sliderStim.started')
                # update status
                sliderStim.status = STARTED
                sliderStim.setAutoDraw(True)
            
            # if sliderStim is active this frame...
            if sliderStim.status == STARTED:
                # update params
                pass
            
            # *textSpaceStim* updates
            
            # if textSpaceStim is starting this frame...
            if textSpaceStim.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                textSpaceStim.frameNStart = frameN  # exact frame index
                textSpaceStim.tStart = t  # local t and not account for scr refresh
                textSpaceStim.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(textSpaceStim, 'tStartRefresh')  # time at next scr refresh
                # update status
                textSpaceStim.status = STARTED
                textSpaceStim.setAutoDraw(True)
            
            # if textSpaceStim is active this frame...
            if textSpaceStim.status == STARTED:
                # update params
                pass
            
            # *spaceStim* updates
            waitOnFlip = False
            
            # if spaceStim is starting this frame...
            if spaceStim.status == NOT_STARTED and sliderStim.rating:
                # keep track of start time/frame for later
                spaceStim.frameNStart = frameN  # exact frame index
                spaceStim.tStart = t  # local t and not account for scr refresh
                spaceStim.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(spaceStim, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'spaceStim.started')
                # update status
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
                if eyetracker:
                    eyetracker.setConnectionState(False)
            
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
        ratingstrials.addData('sliderStim.response', sliderStim.getRating())
        ratingstrials.addData('sliderStim.rt', sliderStim.getRT())
        # check responses
        if spaceStim.keys in ['', [], None]:  # No response was made
            spaceStim.keys = None
        ratingstrials.addData('spaceStim.keys',spaceStim.keys)
        if spaceStim.keys != None:  # we had a response
            ratingstrials.addData('spaceStim.rt', spaceStim.rt)
        # the Routine "stimRating" was not non-slip safe, so reset the non-slip timer
        routineTimer.reset()
        thisExp.nextEntry()
        
    # completed 1.0 repeats of 'ratingstrials'
    
    
    # --- Prepare to start Routine "testInstr" ---
    continueRoutine = True
    # update component parameters for each repeat
    spaceTestInstr.keys = []
    spaceTestInstr.rt = []
    _spaceTestInstr_allKeys = []
    # Run 'Begin Routine' code from codeTestInstruction
    print("Start Test Phase. Press Space to continue.")
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
    routineForceEnded = not continueRoutine
    while continueRoutine:
        # get current time
        t = routineTimer.getTime()
        tThisFlip = win.getFutureFlipTime(clock=routineTimer)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        
        # *textTestInstr* updates
        
        # if textTestInstr is starting this frame...
        if textTestInstr.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            textTestInstr.frameNStart = frameN  # exact frame index
            textTestInstr.tStart = t  # local t and not account for scr refresh
            textTestInstr.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(textTestInstr, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'textTestInstr.started')
            # update status
            textTestInstr.status = STARTED
            textTestInstr.setAutoDraw(True)
        
        # if textTestInstr is active this frame...
        if textTestInstr.status == STARTED:
            # update params
            pass
        
        # *spaceTestInstr* updates
        waitOnFlip = False
        
        # if spaceTestInstr is starting this frame...
        if spaceTestInstr.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            spaceTestInstr.frameNStart = frameN  # exact frame index
            spaceTestInstr.tStart = t  # local t and not account for scr refresh
            spaceTestInstr.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(spaceTestInstr, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'spaceTestInstr.started')
            # update status
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
            if eyetracker:
                eyetracker.setConnectionState(False)
        
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
        # update component parameters for each repeat
        # Run 'Begin Routine' code from codeFixate
        logging.log(level=logging.INFO, msg=f'FixationCross')
        eyetracker.setRecordingState(True)
        
        if currentLoop.name == "trials":
            eyetracker.sendMessage("VP %s TRIALID %d CONDITION %s"%(expInfo['participant'], trials.thisN+1, trialtype))
            print("VP %s TRIALID %d CONDITION %s"%(expInfo['participant'], trials.thisN+1, trialtype))
            
        if currentLoop.name == "testtrials":
            eyetracker.sendMessage("VP %s TESTTRIALID %d CONDITION %s"%(expInfo['participant'], testtrials.thisN+1, stimtype))
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
        routineForceEnded = not continueRoutine
        while continueRoutine and routineTimer.getTime() < 1.0:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # *fixcross* updates
            
            # if fixcross is starting this frame...
            if fixcross.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                fixcross.frameNStart = frameN  # exact frame index
                fixcross.tStart = t  # local t and not account for scr refresh
                fixcross.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(fixcross, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'fixcross.started')
                # update status
                fixcross.status = STARTED
                fixcross.setAutoDraw(True)
            
            # if fixcross is active this frame...
            if fixcross.status == STARTED:
                # update params
                pass
            
            # if fixcross is stopping this frame...
            if fixcross.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > fixcross.tStartRefresh + 1-frameTolerance:
                    # keep track of stop time/frame for later
                    fixcross.tStop = t  # not accounting for scr refresh
                    fixcross.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'fixcross.stopped')
                    # update status
                    fixcross.status = FINISHED
                    fixcross.setAutoDraw(False)
            
            # check for quit (typically the Esc key)
            if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
                core.quit()
                if eyetracker:
                    eyetracker.setConnectionState(False)
            
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
        # update component parameters for each repeat
        imageTest.setPos((0, 0))
        imageTest.setImage(eval(stimtype))
        # Run 'Begin Routine' code from codeTesttrial
        from psychopy import logging
        logging.log(level=logging.INFO, msg=f'TestImageOnset_{stimtype}')
        # keep track of which components have finished
        testtrialComponents = [imageTest, portTestImage]
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
        routineForceEnded = not continueRoutine
        while continueRoutine and routineTimer.getTime() < 10.0:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # *imageTest* updates
            
            # if imageTest is starting this frame...
            if imageTest.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                # keep track of start time/frame for later
                imageTest.frameNStart = frameN  # exact frame index
                imageTest.tStart = t  # local t and not account for scr refresh
                imageTest.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(imageTest, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'imageTest.started')
                # update status
                imageTest.status = STARTED
                imageTest.setAutoDraw(True)
            
            # if imageTest is active this frame...
            if imageTest.status == STARTED:
                # update params
                imageTest.setSize(imagesize_test, log=False)
            
            # if imageTest is stopping this frame...
            if imageTest.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > imageTest.tStartRefresh + 10-frameTolerance:
                    # keep track of stop time/frame for later
                    imageTest.tStop = t  # not accounting for scr refresh
                    imageTest.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'imageTest.stopped')
                    # update status
                    imageTest.status = FINISHED
                    imageTest.setAutoDraw(False)
            # *portTestImage* updates
            
            # if portTestImage is starting this frame...
            if portTestImage.status == NOT_STARTED and t >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                portTestImage.frameNStart = frameN  # exact frame index
                portTestImage.tStart = t  # local t and not account for scr refresh
                portTestImage.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(portTestImage, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.addData('portTestImage.started', t)
                # update status
                portTestImage.status = STARTED
                portTestImage.status = STARTED
                win.callOnFlip(portTestImage.setData, int(2))
            
            # if portTestImage is stopping this frame...
            if portTestImage.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > portTestImage.tStartRefresh + 0.005-frameTolerance:
                    # keep track of stop time/frame for later
                    portTestImage.tStop = t  # not accounting for scr refresh
                    portTestImage.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.addData('portTestImage.stopped', t)
                    # update status
                    portTestImage.status = FINISHED
                    win.callOnFlip(portTestImage.setData, int(0))
            
            # check for quit (typically the Esc key)
            if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
                core.quit()
                if eyetracker:
                    eyetracker.setConnectionState(False)
            
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
        if portTestImage.status == STARTED:
            win.callOnFlip(portTestImage.setData, int(0))
        # using non-slip timing so subtract the expected duration of this Routine (unless ended on request)
        if routineForceEnded:
            routineTimer.reset()
        else:
            routineTimer.addTime(-10.000000)
        
        # --- Prepare to start Routine "crossEnd" ---
        continueRoutine = True
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
        routineForceEnded = not continueRoutine
        while continueRoutine:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # *fixcrossE* updates
            
            # if fixcrossE is starting this frame...
            if fixcrossE.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                fixcrossE.frameNStart = frameN  # exact frame index
                fixcrossE.tStart = t  # local t and not account for scr refresh
                fixcrossE.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(fixcrossE, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'fixcrossE.started')
                # update status
                fixcrossE.status = STARTED
                fixcrossE.setAutoDraw(True)
            
            # if fixcrossE is active this frame...
            if fixcrossE.status == STARTED:
                # update params
                pass
            
            # if fixcrossE is stopping this frame...
            if fixcrossE.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > fixcrossE.tStartRefresh +  3 + random()-frameTolerance:
                    # keep track of stop time/frame for later
                    fixcrossE.tStop = t  # not accounting for scr refresh
                    fixcrossE.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'fixcrossE.stopped')
                    # update status
                    fixcrossE.status = FINISHED
                    fixcrossE.setAutoDraw(False)
            
            # check for quit (typically the Esc key)
            if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
                core.quit()
                if eyetracker:
                    eyetracker.setConnectionState(False)
            
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
    ratingstrials2 = data.TrialHandler(nReps=1.0, method='random', 
        extraInfo=expInfo, originPath=-1,
        trialList=data.importConditions(stimFile),
        seed=None, name='ratingstrials2')
    thisExp.addLoop(ratingstrials2)  # add the loop to the experiment
    thisRatingstrials2 = ratingstrials2.trialList[0]  # so we can initialise stimuli with some values
    # abbreviate parameter names if possible (e.g. rgb = thisRatingstrials2.rgb)
    if thisRatingstrials2 != None:
        for paramName in thisRatingstrials2:
            exec('{} = thisRatingstrials2[paramName]'.format(paramName))
    
    for thisRatingstrials2 in ratingstrials2:
        currentLoop = ratingstrials2
        # abbreviate parameter names if possible (e.g. rgb = thisRatingstrials2.rgb)
        if thisRatingstrials2 != None:
            for paramName in thisRatingstrials2:
                exec('{} = thisRatingstrials2[paramName]'.format(paramName))
        
        # --- Prepare to start Routine "stimRating" ---
        continueRoutine = True
        # update component parameters for each repeat
        imageRating.setPos((0, 0.1))
        imageRating.setImage(eval(stimtype))
        sliderStim.reset()
        spaceStim.keys = []
        spaceStim.rt = []
        _spaceStim_allKeys = []
        # Run 'Begin Routine' code from codeStimRating
        logging.log(level=logging.INFO, msg=f'StimRating_{stimtype}')
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
        routineForceEnded = not continueRoutine
        while continueRoutine:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # *textRateStim* updates
            
            # if textRateStim is starting this frame...
            if textRateStim.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                textRateStim.frameNStart = frameN  # exact frame index
                textRateStim.tStart = t  # local t and not account for scr refresh
                textRateStim.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(textRateStim, 'tStartRefresh')  # time at next scr refresh
                # update status
                textRateStim.status = STARTED
                textRateStim.setAutoDraw(True)
            
            # if textRateStim is active this frame...
            if textRateStim.status == STARTED:
                # update params
                pass
            
            # *imageRating* updates
            
            # if imageRating is starting this frame...
            if imageRating.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                # keep track of start time/frame for later
                imageRating.frameNStart = frameN  # exact frame index
                imageRating.tStart = t  # local t and not account for scr refresh
                imageRating.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(imageRating, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'imageRating.started')
                # update status
                imageRating.status = STARTED
                imageRating.setAutoDraw(True)
            
            # if imageRating is active this frame...
            if imageRating.status == STARTED:
                # update params
                imageRating.setSize(imagesize_rating, log=False)
            
            # *textUnpleasant* updates
            
            # if textUnpleasant is starting this frame...
            if textUnpleasant.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                textUnpleasant.frameNStart = frameN  # exact frame index
                textUnpleasant.tStart = t  # local t and not account for scr refresh
                textUnpleasant.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(textUnpleasant, 'tStartRefresh')  # time at next scr refresh
                # update status
                textUnpleasant.status = STARTED
                textUnpleasant.setAutoDraw(True)
            
            # if textUnpleasant is active this frame...
            if textUnpleasant.status == STARTED:
                # update params
                pass
            
            # *textPleasant* updates
            
            # if textPleasant is starting this frame...
            if textPleasant.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                textPleasant.frameNStart = frameN  # exact frame index
                textPleasant.tStart = t  # local t and not account for scr refresh
                textPleasant.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(textPleasant, 'tStartRefresh')  # time at next scr refresh
                # update status
                textPleasant.status = STARTED
                textPleasant.setAutoDraw(True)
            
            # if textPleasant is active this frame...
            if textPleasant.status == STARTED:
                # update params
                pass
            
            # *sliderStim* updates
            
            # if sliderStim is starting this frame...
            if sliderStim.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                sliderStim.frameNStart = frameN  # exact frame index
                sliderStim.tStart = t  # local t and not account for scr refresh
                sliderStim.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sliderStim, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'sliderStim.started')
                # update status
                sliderStim.status = STARTED
                sliderStim.setAutoDraw(True)
            
            # if sliderStim is active this frame...
            if sliderStim.status == STARTED:
                # update params
                pass
            
            # *textSpaceStim* updates
            
            # if textSpaceStim is starting this frame...
            if textSpaceStim.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                textSpaceStim.frameNStart = frameN  # exact frame index
                textSpaceStim.tStart = t  # local t and not account for scr refresh
                textSpaceStim.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(textSpaceStim, 'tStartRefresh')  # time at next scr refresh
                # update status
                textSpaceStim.status = STARTED
                textSpaceStim.setAutoDraw(True)
            
            # if textSpaceStim is active this frame...
            if textSpaceStim.status == STARTED:
                # update params
                pass
            
            # *spaceStim* updates
            waitOnFlip = False
            
            # if spaceStim is starting this frame...
            if spaceStim.status == NOT_STARTED and sliderStim.rating:
                # keep track of start time/frame for later
                spaceStim.frameNStart = frameN  # exact frame index
                spaceStim.tStart = t  # local t and not account for scr refresh
                spaceStim.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(spaceStim, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'spaceStim.started')
                # update status
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
                if eyetracker:
                    eyetracker.setConnectionState(False)
            
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
        ratingstrials2.addData('sliderStim.response', sliderStim.getRating())
        ratingstrials2.addData('sliderStim.rt', sliderStim.getRT())
        # check responses
        if spaceStim.keys in ['', [], None]:  # No response was made
            spaceStim.keys = None
        ratingstrials2.addData('spaceStim.keys',spaceStim.keys)
        if spaceStim.keys != None:  # we had a response
            ratingstrials2.addData('spaceStim.rt', spaceStim.rt)
        # the Routine "stimRating" was not non-slip safe, so reset the non-slip timer
        routineTimer.reset()
        thisExp.nextEntry()
        
    # completed 1.0 repeats of 'ratingstrials2'
    
# completed 1.0 repeats of 'blocks'


# --- Prepare to start Routine "end" ---
continueRoutine = True
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
routineForceEnded = not continueRoutine
while continueRoutine and routineTimer.getTime() < 10.0:
    # get current time
    t = routineTimer.getTime()
    tThisFlip = win.getFutureFlipTime(clock=routineTimer)
    tThisFlipGlobal = win.getFutureFlipTime(clock=None)
    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
    # update/draw components on each frame
    
    # *textEnd* updates
    
    # if textEnd is starting this frame...
    if textEnd.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
        # keep track of start time/frame for later
        textEnd.frameNStart = frameN  # exact frame index
        textEnd.tStart = t  # local t and not account for scr refresh
        textEnd.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(textEnd, 'tStartRefresh')  # time at next scr refresh
        # update status
        textEnd.status = STARTED
        textEnd.setAutoDraw(True)
    
    # if textEnd is active this frame...
    if textEnd.status == STARTED:
        # update params
        pass
    
    # if textEnd is stopping this frame...
    if textEnd.status == STARTED:
        # is it time to stop? (based on global clock, using actual start)
        if tThisFlipGlobal > textEnd.tStartRefresh + 10-frameTolerance:
            # keep track of stop time/frame for later
            textEnd.tStop = t  # not accounting for scr refresh
            textEnd.frameNStop = frameN  # exact frame index
            # update status
            textEnd.status = FINISHED
            textEnd.setAutoDraw(False)
    
    # check for quit (typically the Esc key)
    if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
        core.quit()
        if eyetracker:
            eyetracker.setConnectionState(False)
    
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
# using non-slip timing so subtract the expected duration of this Routine (unless ended on request)
if routineForceEnded:
    routineTimer.reset()
else:
    routineTimer.addTime(-10.000000)

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
