#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This experiment was created using PsychoPy3 Experiment Builder (v2022.2.4),
    on Mai 16, 2024, at 15:33
If you publish work using this script the most relevant publication is:

    Peirce J, Gray JR, Simpson S, MacAskill M, Höchenberger R, Sogo H, Kastman E, Lindeløv JK. (2019) 
        PsychoPy2: Experiments in behavior made easy Behav Res 51: 195. 
        https://doi.org/10.3758/s13428-018-01193-y

"""

import psychopy
psychopy.useVersion('2022.2.4')


# --- Import packages ---
from psychopy import locale_setup
from psychopy import prefs
prefs.hardware['audioLib'] = 'pygame'
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

# Run 'Before Experiment' code from code_end
jittered_duration_cross = 3.5
# Run 'Before Experiment' code from codeBlank
jittered_duration_blank = 1.5
# Run 'Before Experiment' code from codeBlank
jittered_duration_blank = 1.5
# Run 'Before Experiment' code from codeTrial
feedback_color = (0,0,0)
feedback_opacity = 0
feedback_audio = "audio/silence.wav"
cursorcolor="white"
log_msg = ""

rectsize = (0.2, 0.2)
imagesize_test = [1, 1]

score = 0
points = 0

looked_at = False
paused = False
outcome = ""
# Run 'Before Experiment' code from codeFeedback
soundFeedback = sound.Sound('A', secs=1, stereo=True, hamming=True, name='soundFeedback')
# Run 'Before Experiment' code from code_end
jittered_duration_cross = 3.5
# Run 'Before Experiment' code from codeBlank
jittered_duration_blank = 1.5
# Run 'Before Experiment' code from codeBlank
jittered_duration_blank = 1.5
# Run 'Before Experiment' code from codeTrial
feedback_color = (0,0,0)
feedback_opacity = 0
feedback_audio = "audio/silence.wav"
cursorcolor="white"
log_msg = ""

rectsize = (0.2, 0.2)
imagesize_test = [1, 1]

score = 0
points = 0

looked_at = False
paused = False
outcome = ""
# Run 'Before Experiment' code from codeFeedback
soundFeedback = sound.Sound('A', secs=1, stereo=True, hamming=True, name='soundFeedback')
# Run 'Before Experiment' code from code_end
jittered_duration_cross = 3.5
# Run 'Before Experiment' code from codeBlankTest
jittered_duration_blank = 1.5
# Run 'Before Experiment' code from codeTestTrial
feedback_color = (0,0,0)
feedback_opacity = 0
feedback_audio = "audio/silence.wav"
cursorcolor="white"
log_msg = ""

rectsize = (0.2, 0.2)
imagesize_test = [1, 1]
feedback_position = (0, 0)

score = 0
points = 0

looked_at_roi1 = False
looked_at_roi2 = False
looked_at_roi3 = False
looked_at_roi4 = False
paused = False
outcome = ""
# Run 'Before Experiment' code from codeFeedbackTest
soundFeedback = sound.Sound('A', secs=1, stereo=True, hamming=True, name='soundFeedback')
# Run 'Before Experiment' code from code_end
jittered_duration_cross = 3.5
# Run 'Before Experiment' code from codeBlankTest
jittered_duration_blank = 1.5
# Run 'Before Experiment' code from codeTestTrial
feedback_color = (0,0,0)
feedback_opacity = 0
feedback_audio = "audio/silence.wav"
cursorcolor="white"
log_msg = ""

rectsize = (0.2, 0.2)
imagesize_test = [1, 1]
feedback_position = (0, 0)

score = 0
points = 0

looked_at_roi1 = False
looked_at_roi2 = False
looked_at_roi3 = False
looked_at_roi4 = False
paused = False
outcome = ""
# Run 'Before Experiment' code from codeFeedbackTest
soundFeedback = sound.Sound('A', secs=1, stereo=True, hamming=True, name='soundFeedback')
# Run 'Before Experiment' code from code_end
jittered_duration_cross = 3.5


# Ensure that relative paths start from the same directory as this script
_thisDir = os.path.dirname(os.path.abspath(__file__))
os.chdir(_thisDir)
# Store info about the experiment session
psychopyVersion = '2022.2.4'
expName = 'attentional_competition_task'  # from the Builder filename that created this script
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
    originPath='C:\\Users\\sag22id\\Documents\\Projects\\GCA\\gca_avoidance\\Study 2\\Experiment\\attentional_competition_task\\attentional_competition_task_lastrun.py',
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
ioServer = io.launchHubServer(window=win, experiment_code='attentional_competition_task', session_code=ioSession, datastore_name=filename, **ioConfig)
eyetracker = ioServer.getDevice('tracker')

# create a default keyboard (e.g. to check for escape)
defaultKeyboard = keyboard.Keyboard(backend='iohub')

# --- Initialize components for Routine "welcome" ---
textStart = visual.TextStim(win=win, name='textStart',
    text='Sehr geehrter Teilnehmer, sehr geehrte Teilnehmerin,\nvielen Dank für Ihre Teilnahme an unserem Experiment.\nJetzt beginnt der zweite Teil.\n\nBitte drücken Sie die Leertaste, um zu starten.',
    font='Open Sans',
    pos=(0, 0), height=0.06, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=0.0);
spaceStart = keyboard.Keyboard()

# --- Initialize components for Routine "startETCalibration" ---
text_ETCalibration = visual.TextStim(win=win, name='text_ETCalibration',
    text='Wir starten mit der Kalibrierung des Eye-Trackers.\n\nBitte drücken Sie wieder die Leertaste, um diese zu beginnen.',
    font='Open Sans',
    pos=(0, 0), height=0.06, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=0.0);
spaceETCalibration = keyboard.Keyboard()

# --- Initialize components for Routine "startExp" ---
textStartExp = visual.TextStim(win=win, name='textStartExp',
    text='Die Kalibrierung ist abgeschlossen.\n\nWährend des Experiments werden Sie mehrere Bilder sehen.\nBitte geben Sie nun an, wie Sie sich bei der Betrachtung der Bilder fühlen.\n\nDrücken Sie die Leertaste, um zu starten.',
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

# --- Initialize components for Routine "startSearchTask" ---
textStartSearch = visual.TextStim(win=win, name='textStartSearch',
    text='Zuerst möchten wir Ihnen das grundsätzliche Prinzip unserer Studie verdeutlichen: Die präsentierten Bilder reagieren auf Ihr Blickverhalten.\n\nDafür werden Sie nun mehrere Kreise sehen, welche ebenso auf Ihre Blicke reagieren. Probieren Sie dies nun in Ruhe aus.\n\nDrücken Sie die Leertaste, um zu starten.',
    font='Open Sans',
    pos=(0, 0), height=0.06, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=0.0);
spaceStartSearchTask = keyboard.Keyboard()

# --- Initialize components for Routine "cross" ---
fixcross = visual.TextStim(win=win, name='fixcross',
    text='+',
    font='Open Sans',
    pos=(0, 0), height=0.1, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-1.0);

# --- Initialize components for Routine "search_task" ---
# Run 'Begin Experiment' code from codeSearchTask
x_min = -.9
x_max = .9
y_min = -.8
y_max = .7
x_n = 9
y_n = 4
positions = []

for idx_x in range(x_n):
    for idx_y in range(y_n):
        positions.append([idx_x / (x_n-1) * (abs(x_min) + x_max) + x_min, idx_y/ (y_n-1) * (abs(y_min) + y_max) + y_min])

jitter = .06
elements_per_group = 15
elements = []

width = win.size[0]
height = win.size[1]
width_circle = 0.07
height_circle = 0.07 * (width/height)
gazeCursorSearch = visual.ShapeStim(
    win=win, name='gazeCursorSearch',
    size=(0.01, 0.02), vertices='circle',
    ori=0.0, pos=[0,0], anchor='center',
    lineWidth=1.0,     colorSpace='rgb',  lineColor=None, fillColor=[-1.0000, -1.0000, 1.0000],
    opacity=1.0, depth=-1.0, interpolate=True)

# --- Initialize components for Routine "crossEnd" ---
fixcrossE = visual.TextStim(win=win, name='fixcrossE',
    text='+',
    font='Open Sans',
    pos=(0, 0), height=0.1, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-1.0);

# --- Initialize components for Routine "startTask" ---
textStartTask = visual.TextStim(win=win, name='textStartTask',
    text='Wir starten nun mit dem Experiment.\n\nSie werden in jedem Durchgang in einer der vier Ecken des Bildschirms ein Bild präsentiert bekommen. Über Ihr Blickverhalten können Sie eine Belohnung in Form von Punkten erhalten. Allerdings lauert auch die Gefahr eines Verlustes von Punkten. Ihr Ziel ist es Ihren Belohnungsscore zu maximieren.\n\nBitte fixieren Sie am Anfang jedes Durchgangs das Fixationskreuz.\n\nDrücken Sie die Leertaste, um zu starten.\n',
    font='Open Sans',
    pos=(0, 0), height=0.06, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-1.0);
spaceStartTask = keyboard.Keyboard()

# --- Initialize components for Routine "blank" ---
blankScreen = visual.TextStim(win=win, name='blankScreen',
    text=None,
    font='Open Sans',
    pos=(0, 0), height=0.05, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-1.0);

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
    depth=-1.0);

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
    win=win, name='gazeCursor',
    size=(0.01, 0.02), vertices='circle',
    ori=0.0, pos=[0,0], anchor='center',
    lineWidth=1.0,     colorSpace='rgb',  lineColor=None, fillColor='white',
    opacity=1.0, depth=-3.0, interpolate=True)

# --- Initialize components for Routine "feedback" ---
textPoints = visual.TextStim(win=win, name='textPoints',
    text='',
    font='Open Sans',
    pos=(0, 0.15), height=0.1, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-1.0);
textScore = visual.TextStim(win=win, name='textScore',
    text='',
    font='Open Sans',
    pos=(0, 0), height=0.06, wrapWidth=None, ori=0.0, 
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
    win=win, name='gazeCursor_Feedback',
    size=(0.01, 0.02), vertices='circle',
    ori=0.0, pos=[0,0], anchor='center',
    lineWidth=1.0,     colorSpace='rgb',  lineColor=None, fillColor='white',
    opacity=1.0, depth=-5.0, interpolate=True)
soundFeedback = sound.Sound('A', secs=1, stereo=True, hamming=True,
    name='soundFeedback')
soundFeedback.setVolume(1.0)

# --- Initialize components for Routine "crossEnd" ---
fixcrossE = visual.TextStim(win=win, name='fixcrossE',
    text='+',
    font='Open Sans',
    pos=(0, 0), height=0.1, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-1.0);

# --- Initialize components for Routine "repeatTask" ---
textRepeatTask = visual.TextStim(win=win, name='textRepeatTask',
    text='--- Kurze Pause ---\n\nIm nächsten Block ist Ihre Aufgabe weiterhin Ihren Belohnungsscore zu maximieren.',
    font='Open Sans',
    pos=(0, 0.3), height=0.06, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-1.0);
textboxRepeatInstruction = visual.TextBox2(
     win, text='Sie können den Verlust von Punkten vermeiden, in dem Sie das entsprechende Bild nicht anschauen.', font='Open Sans',
     pos=(0, 0),     letterHeight=0.06,
     size=(None, None), borderWidth=2.0,
     color='white', colorSpace='rgb',
     opacity=None,
     bold=True, italic=False,
     lineSpacing=1.0,
     padding=0.0, alignment='center',
     anchor='center',
     fillColor=None, borderColor=None,
     flipHoriz=False, flipVert=False, languageStyle='LTR',
     editable=False,
     name='textboxRepeatInstruction',
     autoLog=True,
)
textRepeatSpace = visual.TextStim(win=win, name='textRepeatSpace',
    text='Bitte fixieren Sie am Anfang jedes Durchgangs das Fixationskreuz.\n\nDrücken Sie die Leertaste, um weiterzumachen.',
    font='Open Sans',
    pos=(0, -0.3), height=0.06, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-3.0);
spaceRepeatTask = keyboard.Keyboard()

# --- Initialize components for Routine "blank" ---
blankScreen = visual.TextStim(win=win, name='blankScreen',
    text=None,
    font='Open Sans',
    pos=(0, 0), height=0.05, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-1.0);

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
    depth=-1.0);

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
    win=win, name='gazeCursor',
    size=(0.01, 0.02), vertices='circle',
    ori=0.0, pos=[0,0], anchor='center',
    lineWidth=1.0,     colorSpace='rgb',  lineColor=None, fillColor='white',
    opacity=1.0, depth=-3.0, interpolate=True)

# --- Initialize components for Routine "feedback" ---
textPoints = visual.TextStim(win=win, name='textPoints',
    text='',
    font='Open Sans',
    pos=(0, 0.15), height=0.1, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-1.0);
textScore = visual.TextStim(win=win, name='textScore',
    text='',
    font='Open Sans',
    pos=(0, 0), height=0.06, wrapWidth=None, ori=0.0, 
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
    win=win, name='gazeCursor_Feedback',
    size=(0.01, 0.02), vertices='circle',
    ori=0.0, pos=[0,0], anchor='center',
    lineWidth=1.0,     colorSpace='rgb',  lineColor=None, fillColor='white',
    opacity=1.0, depth=-5.0, interpolate=True)
soundFeedback = sound.Sound('A', secs=1, stereo=True, hamming=True,
    name='soundFeedback')
soundFeedback.setVolume(1.0)

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

# --- Initialize components for Routine "blank_test" ---
blankScreenTest = visual.TextStim(win=win, name='blankScreenTest',
    text=None,
    font='Open Sans',
    pos=(0, 0), height=0.05, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-1.0);

# --- Initialize components for Routine "test_trial" ---
test_image1 = visual.ImageStim(
    win=win,
    name='test_image1', 
    image='sin', mask=None, anchor='center',
    ori=0.0, pos=[0,0], size=None,
    color=[1,1,1], colorSpace='rgb', opacity=None,
    flipHoriz=False, flipVert=False,
    texRes=128.0, interpolate=True, depth=-1.0)
test_roi1 = visual.ROI(win, name='test_roi1', device=eyetracker,
    debug=False,
    shape='rectangle',
    pos=[0,0], size=1.0, anchor='center', ori=0.0)
test_image2 = visual.ImageStim(
    win=win,
    name='test_image2', 
    image='sin', mask=None, anchor='center',
    ori=0.0, pos=[0,0], size=None,
    color=[1,1,1], colorSpace='rgb', opacity=None,
    flipHoriz=False, flipVert=False,
    texRes=128.0, interpolate=True, depth=-3.0)
test_roi2 = visual.ROI(win, name='test_roi2', device=eyetracker,
    debug=False,
    shape='rectangle',
    pos=[0,0], size=1.0, anchor='center', ori=0.0)
test_image3 = visual.ImageStim(
    win=win,
    name='test_image3', 
    image='sin', mask=None, anchor='center',
    ori=0.0, pos=[0,0], size=None,
    color=[1,1,1], colorSpace='rgb', opacity=None,
    flipHoriz=False, flipVert=False,
    texRes=128.0, interpolate=True, depth=-5.0)
test_roi3 = visual.ROI(win, name='test_roi3', device=eyetracker,
    debug=False,
    shape='rectangle',
    pos=[0,0], size=1.0, anchor='center', ori=0.0)
test_image4 = visual.ImageStim(
    win=win,
    name='test_image4', 
    image='sin', mask=None, anchor='center',
    ori=0.0, pos=[0,0], size=None,
    color=[1,1,1], colorSpace='rgb', opacity=None,
    flipHoriz=False, flipVert=False,
    texRes=128.0, interpolate=True, depth=-7.0)
test_roi4 = visual.ROI(win, name='test_roi4', device=eyetracker,
    debug=False,
    shape='rectangle',
    pos=[0,0], size=1.0, anchor='center', ori=0.0)
gazeCursorTest = visual.ShapeStim(
    win=win, name='gazeCursorTest',
    size=(0.01, 0.02), vertices='circle',
    ori=0.0, pos=[0,0], anchor='center',
    lineWidth=1.0,     colorSpace='rgb',  lineColor=None, fillColor='white',
    opacity=1.0, depth=-9.0, interpolate=True)

# --- Initialize components for Routine "test_feedback" ---
textScoreTest = visual.TextStim(win=win, name='textScoreTest',
    text='',
    font='Open Sans',
    pos=(0, 0), height=0.06, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-1.0);
polygonFeedbackTest = visual.Rect(
    win=win, name='polygonFeedbackTest',
    width=[1.0, 1.0][0], height=[1.0, 1.0][1],
    ori=0.0, pos=[0,0], anchor='center',
    lineWidth=100.0,     colorSpace='rgb',  lineColor='white', fillColor=[0.0000, 0.0000, 0.0000],
    opacity=1.0, depth=-2.0, interpolate=True)
test_image_feedback1 = visual.ImageStim(
    win=win,
    name='test_image_feedback1', 
    image='sin', mask=None, anchor='center',
    ori=0.0, pos=[0,0], size=None,
    color=[1,1,1], colorSpace='rgb', opacity=None,
    flipHoriz=False, flipVert=False,
    texRes=128.0, interpolate=True, depth=-3.0)
test_image_feedback2 = visual.ImageStim(
    win=win,
    name='test_image_feedback2', 
    image='sin', mask=None, anchor='center',
    ori=0.0, pos=[0,0], size=None,
    color=[1,1,1], colorSpace='rgb', opacity=None,
    flipHoriz=False, flipVert=False,
    texRes=128.0, interpolate=True, depth=-4.0)
test_iamge_feedback3 = visual.ImageStim(
    win=win,
    name='test_iamge_feedback3', 
    image='sin', mask=None, anchor='center',
    ori=0.0, pos=[0,0], size=None,
    color=[1,1,1], colorSpace='rgb', opacity=None,
    flipHoriz=False, flipVert=False,
    texRes=128.0, interpolate=True, depth=-5.0)
test_iamge_feedback4 = visual.ImageStim(
    win=win,
    name='test_iamge_feedback4', 
    image='sin', mask=None, anchor='center',
    ori=0.0, pos=[0,0], size=None,
    color=[1,1,1], colorSpace='rgb', opacity=None,
    flipHoriz=False, flipVert=False,
    texRes=128.0, interpolate=True, depth=-6.0)
gazeCursorTestFeedback = visual.ShapeStim(
    win=win, name='gazeCursorTestFeedback',
    size=(0.01, 0.02), vertices='circle',
    ori=0.0, pos=[0,0], anchor='center',
    lineWidth=1.0,     colorSpace='rgb',  lineColor=None, fillColor='white',
    opacity=1.0, depth=-7.0, interpolate=True)
soundFeedbackTest = sound.Sound('A', secs=0.5, stereo=True, hamming=True,
    name='soundFeedbackTest')
soundFeedbackTest.setVolume(1.0)

# --- Initialize components for Routine "crossEnd" ---
fixcrossE = visual.TextStim(win=win, name='fixcrossE',
    text='+',
    font='Open Sans',
    pos=(0, 0), height=0.1, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-1.0);

# --- Initialize components for Routine "cross" ---
fixcross = visual.TextStim(win=win, name='fixcross',
    text='+',
    font='Open Sans',
    pos=(0, 0), height=0.1, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-1.0);

# --- Initialize components for Routine "blank_test" ---
blankScreenTest = visual.TextStim(win=win, name='blankScreenTest',
    text=None,
    font='Open Sans',
    pos=(0, 0), height=0.05, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-1.0);

# --- Initialize components for Routine "test_trial" ---
test_image1 = visual.ImageStim(
    win=win,
    name='test_image1', 
    image='sin', mask=None, anchor='center',
    ori=0.0, pos=[0,0], size=None,
    color=[1,1,1], colorSpace='rgb', opacity=None,
    flipHoriz=False, flipVert=False,
    texRes=128.0, interpolate=True, depth=-1.0)
test_roi1 = visual.ROI(win, name='test_roi1', device=eyetracker,
    debug=False,
    shape='rectangle',
    pos=[0,0], size=1.0, anchor='center', ori=0.0)
test_image2 = visual.ImageStim(
    win=win,
    name='test_image2', 
    image='sin', mask=None, anchor='center',
    ori=0.0, pos=[0,0], size=None,
    color=[1,1,1], colorSpace='rgb', opacity=None,
    flipHoriz=False, flipVert=False,
    texRes=128.0, interpolate=True, depth=-3.0)
test_roi2 = visual.ROI(win, name='test_roi2', device=eyetracker,
    debug=False,
    shape='rectangle',
    pos=[0,0], size=1.0, anchor='center', ori=0.0)
test_image3 = visual.ImageStim(
    win=win,
    name='test_image3', 
    image='sin', mask=None, anchor='center',
    ori=0.0, pos=[0,0], size=None,
    color=[1,1,1], colorSpace='rgb', opacity=None,
    flipHoriz=False, flipVert=False,
    texRes=128.0, interpolate=True, depth=-5.0)
test_roi3 = visual.ROI(win, name='test_roi3', device=eyetracker,
    debug=False,
    shape='rectangle',
    pos=[0,0], size=1.0, anchor='center', ori=0.0)
test_image4 = visual.ImageStim(
    win=win,
    name='test_image4', 
    image='sin', mask=None, anchor='center',
    ori=0.0, pos=[0,0], size=None,
    color=[1,1,1], colorSpace='rgb', opacity=None,
    flipHoriz=False, flipVert=False,
    texRes=128.0, interpolate=True, depth=-7.0)
test_roi4 = visual.ROI(win, name='test_roi4', device=eyetracker,
    debug=False,
    shape='rectangle',
    pos=[0,0], size=1.0, anchor='center', ori=0.0)
gazeCursorTest = visual.ShapeStim(
    win=win, name='gazeCursorTest',
    size=(0.01, 0.02), vertices='circle',
    ori=0.0, pos=[0,0], anchor='center',
    lineWidth=1.0,     colorSpace='rgb',  lineColor=None, fillColor='white',
    opacity=1.0, depth=-9.0, interpolate=True)

# --- Initialize components for Routine "test_feedback" ---
textScoreTest = visual.TextStim(win=win, name='textScoreTest',
    text='',
    font='Open Sans',
    pos=(0, 0), height=0.06, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-1.0);
polygonFeedbackTest = visual.Rect(
    win=win, name='polygonFeedbackTest',
    width=[1.0, 1.0][0], height=[1.0, 1.0][1],
    ori=0.0, pos=[0,0], anchor='center',
    lineWidth=100.0,     colorSpace='rgb',  lineColor='white', fillColor=[0.0000, 0.0000, 0.0000],
    opacity=1.0, depth=-2.0, interpolate=True)
test_image_feedback1 = visual.ImageStim(
    win=win,
    name='test_image_feedback1', 
    image='sin', mask=None, anchor='center',
    ori=0.0, pos=[0,0], size=None,
    color=[1,1,1], colorSpace='rgb', opacity=None,
    flipHoriz=False, flipVert=False,
    texRes=128.0, interpolate=True, depth=-3.0)
test_image_feedback2 = visual.ImageStim(
    win=win,
    name='test_image_feedback2', 
    image='sin', mask=None, anchor='center',
    ori=0.0, pos=[0,0], size=None,
    color=[1,1,1], colorSpace='rgb', opacity=None,
    flipHoriz=False, flipVert=False,
    texRes=128.0, interpolate=True, depth=-4.0)
test_iamge_feedback3 = visual.ImageStim(
    win=win,
    name='test_iamge_feedback3', 
    image='sin', mask=None, anchor='center',
    ori=0.0, pos=[0,0], size=None,
    color=[1,1,1], colorSpace='rgb', opacity=None,
    flipHoriz=False, flipVert=False,
    texRes=128.0, interpolate=True, depth=-5.0)
test_iamge_feedback4 = visual.ImageStim(
    win=win,
    name='test_iamge_feedback4', 
    image='sin', mask=None, anchor='center',
    ori=0.0, pos=[0,0], size=None,
    color=[1,1,1], colorSpace='rgb', opacity=None,
    flipHoriz=False, flipVert=False,
    texRes=128.0, interpolate=True, depth=-6.0)
gazeCursorTestFeedback = visual.ShapeStim(
    win=win, name='gazeCursorTestFeedback',
    size=(0.01, 0.02), vertices='circle',
    ori=0.0, pos=[0,0], anchor='center',
    lineWidth=1.0,     colorSpace='rgb',  lineColor=None, fillColor='white',
    opacity=1.0, depth=-7.0, interpolate=True)
soundFeedbackTest = sound.Sound('A', secs=0.5, stereo=True, hamming=True,
    name='soundFeedbackTest')
soundFeedbackTest.setVolume(1.0)

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
    text='Sehr geehrter Teilnehmer, sehr geehrte Teilnehmerin,\ndieser Teil des Experiments ist beendet. \n\nBitte melden Sie sich bei der Versuchsleitung.',
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
# the Routine "welcome" was not non-slip safe, so reset the non-slip timer
routineTimer.reset()

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
        text_ETCalibration.setAutoDraw(True)
    
    # *spaceETCalibration* updates
    waitOnFlip = False
    if spaceETCalibration.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        # keep track of start time/frame for later
        spaceETCalibration.frameNStart = frameN  # exact frame index
        spaceETCalibration.tStart = t  # local t and not account for scr refresh
        spaceETCalibration.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(spaceETCalibration, 'tStartRefresh')  # time at next scr refresh
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
        logging.log(level=logging.INFO, msg=f'StimRating_{stim}')
        print("Rating of Stimulus: %s"%(stim))
        win.mouseVisible = True
        imageRating.setPos((0, 0.2))
        imageRating.setImage(eval(stim))
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
        # the Routine "stimRating" was not non-slip safe, so reset the non-slip timer
        routineTimer.reset()
        thisExp.nextEntry()
        
    # completed 1.0 repeats of 'ratingtrials1'
    
    
    # --- Prepare to start Routine "startSearchTask" ---
    continueRoutine = True
    routineForceEnded = False
    # update component parameters for each repeat
    spaceStartSearchTask.keys = []
    spaceStartSearchTask.rt = []
    _spaceStartSearchTask_allKeys = []
    # Run 'Begin Routine' code from codeStartSearchTask
    print("Start Experiment Screen. Press Space to continue.")
    win.mouseVisible = False
    # keep track of which components have finished
    startSearchTaskComponents = [textStartSearch, spaceStartSearchTask]
    for thisComponent in startSearchTaskComponents:
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
    
    # --- Run Routine "startSearchTask" ---
    while continueRoutine:
        # get current time
        t = routineTimer.getTime()
        tThisFlip = win.getFutureFlipTime(clock=routineTimer)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        
        # *textStartSearch* updates
        if textStartSearch.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            textStartSearch.frameNStart = frameN  # exact frame index
            textStartSearch.tStart = t  # local t and not account for scr refresh
            textStartSearch.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(textStartSearch, 'tStartRefresh')  # time at next scr refresh
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'textStartSearch.started')
            textStartSearch.setAutoDraw(True)
        
        # *spaceStartSearchTask* updates
        waitOnFlip = False
        if spaceStartSearchTask.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            spaceStartSearchTask.frameNStart = frameN  # exact frame index
            spaceStartSearchTask.tStart = t  # local t and not account for scr refresh
            spaceStartSearchTask.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(spaceStartSearchTask, 'tStartRefresh')  # time at next scr refresh
            spaceStartSearchTask.status = STARTED
            # keyboard checking is just starting
            waitOnFlip = True
            win.callOnFlip(spaceStartSearchTask.clock.reset)  # t=0 on next screen flip
            win.callOnFlip(spaceStartSearchTask.clearEvents, eventType='keyboard')  # clear events on next screen flip
        if spaceStartSearchTask.status == STARTED and not waitOnFlip:
            theseKeys = spaceStartSearchTask.getKeys(keyList=['space', 'enter'], waitRelease=False)
            _spaceStartSearchTask_allKeys.extend(theseKeys)
            if len(_spaceStartSearchTask_allKeys):
                spaceStartSearchTask.keys = _spaceStartSearchTask_allKeys[-1].name  # just the last key pressed
                spaceStartSearchTask.rt = _spaceStartSearchTask_allKeys[-1].rt
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
        for thisComponent in startSearchTaskComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # --- Ending Routine "startSearchTask" ---
    for thisComponent in startSearchTaskComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    # the Routine "startSearchTask" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset()
    
    # set up handler to look after randomisation of conditions etc
    search_task_trials = data.TrialHandler(nReps=2.0, method='random', 
        extraInfo=expInfo, originPath=-1,
        trialList=[None],
        seed=None, name='search_task_trials')
    thisExp.addLoop(search_task_trials)  # add the loop to the experiment
    thisSearch_task_trial = search_task_trials.trialList[0]  # so we can initialise stimuli with some values
    # abbreviate parameter names if possible (e.g. rgb = thisSearch_task_trial.rgb)
    if thisSearch_task_trial != None:
        for paramName in thisSearch_task_trial:
            exec('{} = thisSearch_task_trial[paramName]'.format(paramName))
    
    for thisSearch_task_trial in search_task_trials:
        currentLoop = search_task_trials
        # abbreviate parameter names if possible (e.g. rgb = thisSearch_task_trial.rgb)
        if thisSearch_task_trial != None:
            for paramName in thisSearch_task_trial:
                exec('{} = thisSearch_task_trial[paramName]'.format(paramName))
        
        # --- Prepare to start Routine "cross" ---
        continueRoutine = True
        routineForceEnded = False
        # update component parameters for each repeat
        # Run 'Begin Routine' code from codeFixate
        logging.log(level=logging.INFO, msg=f'FixationCross')
        eyetracker.setRecordingState(True)
         
        if currentLoop.name == "learning_trials":
            eyetracker.sendMessage("Trial " + str(learning_trials.thisN+1) + ", Condition " + stim)
            print("VP %s TRIALID %d CONDITION %s"%(expInfo['participant'], learning_trials.thisN+1, stim))
        
        if currentLoop.name == "trials":
            eyetracker.sendMessage("Trial " + str(trials.thisN+learning_trials.thisN+1) + ", Condition " + stim)
            print("VP %s TRIALID %d CONDITION %s"%(expInfo['participant'], trials.thisN+learning_trials.thisN+1, stim))
        
        if currentLoop.name == "testtrials":
            eyetracker.sendMessage("Test-Trial " + str(testtrials.thisN+1))
            print("VP %s TESTTRIALID %d"%(expInfo['participant'], testtrials.thisN+1))
            
        if currentLoop.name == "testtrials_novelty":
            eyetracker.sendMessage("Test-Trial " + str(testtrials_novelty.thisN+testtrials.thisN+1))
            print("VP %s TESTTRIALID %d"%(expInfo['participant'], testtrials_novelty.thisN+testtrials.thisN+1))
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
                fixcross.setAutoDraw(True)
            if fixcross.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > fixcross.tStartRefresh + 1-frameTolerance:
                    # keep track of stop time/frame for later
                    fixcross.tStop = t  # not accounting for scr refresh
                    fixcross.frameNStop = frameN  # exact frame index
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
        
        # --- Prepare to start Routine "search_task" ---
        continueRoutine = True
        routineForceEnded = False
        # update component parameters for each repeat
        # Run 'Begin Routine' code from codeSearchTask
        logging.log(level=logging.INFO, msg=f'VisualSearchTaskOnset')
        ioServer.sendMessageEvent(text='VisualSearchTaskOnset')
        eyetracker.sendMessage('VisualSearchTaskOnset')
        
        random_number = random()
        
        shuffle(positions)
        target_positions = positions[0:elements_per_group]
        distractor_positions = positions[elements_per_group:2 * elements_per_group]
        
        thisExp.addData('target_positions', target_positions)
        thisExp.addData('distractor_positions', distractor_positions)
        
        myTargets = []
        myTargetROIs = []
        for i in range(elements_per_group):
            position_x = target_positions[i][0]+np.random.uniform(-1, 1)*jitter
            position_y = target_positions[i][1]+np.random.uniform(-1, 1)*jitter
            roi = visual.ROI(win, name='roi', device=eyetracker, debug=False, shape='circle',
                pos=(position_x, position_y),
                size=(width_circle, height_circle), anchor='center', ori=0.0)
            circle = visual.ShapeStim(
                win=win, name='circle',
                size=(width_circle, height_circle), vertices='circle',ori=0.0, anchor='center',
                pos=(position_x, position_y),
                lineWidth=20.0, colorSpace='rgb',lineColor="white", fillColor="white",
                opacity=None, depth=-3, interpolate=True)
            myTargets.append(circle)
            myTargetROIs.append(roi)
        
        for target in myTargets:
            target.setAutoDraw(True)
        
        myDistractors = []
        myDistractorROIs = []
        myDistractorNoises = []
        for i in range(elements_per_group):
            position_x = distractor_positions[i][0]+np.random.uniform(-1, 1)*jitter
            position_y = distractor_positions[i][1]+np.random.uniform(-1, 1)*jitter
            roi = visual.ROI(win, name='roi', device=eyetracker, debug=False, shape='circle',
                pos=(position_x, position_y),
                size=(width_circle, height_circle), anchor='center', ori=0.0)
            circle = visual.ShapeStim(
                win=win, name='circle',
                size=(width_circle, height_circle), vertices='circle',ori=0.0, anchor='center',
                pos=(position_x, position_y),
                lineWidth=20.0, colorSpace='rgb', lineColor="white", fillColor="white",
                opacity=None, depth=-3, interpolate=True)
            noise = sound.Sound('audio/white.wav', secs=0.05, stereo=True, hamming=True, name='noise')
            noise.setVolume(1.0)
            myDistractors.append(circle)
            myDistractorROIs.append(roi)
            myDistractorNoises.append(noise)
        
        for distractor in myDistractors:
            distractor.setAutoDraw(True)
            
        
        # keep track of which components have finished
        search_taskComponents = [gazeCursorSearch]
        for thisComponent in search_taskComponents:
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
        
        # --- Run Routine "search_task" ---
        while continueRoutine and routineTimer.getTime() < 30.0:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            # Run 'Each Frame' code from codeSearchTask
            for i, roiTarget in enumerate(myTargetROIs):
                if roiTarget.isLookedIn:
                    myTargets[i].fillColor = "green"
                    myTargets[i].lineColor = "green"
                else:
                    myTargets[i].fillColor = "white"
            
            for i, roiDistractor in enumerate(myDistractorROIs):
                if roiDistractor.isLookedIn:
                    myDistractors[i].fillColor = "red"
                    myDistractors[i].lineColor = "red"
                    myDistractorNoises[i].play()
                else:
                    myDistractors[i].fillColor = "white"
                    myDistractorNoises[i].stop()
            
            # Use "p" to pause experiment
            if event.getKeys(keyList=["p"]) and not paused:
                print("Experiment Paused - Press 'p' to continue.")
                paused = True
                event.waitKeys(keyList=["p"], clearEvents=True)
                
                print("Experiment Continued")
                paused = False
                event.clearEvents(eventType="keyboard")
                
            
            # *gazeCursorSearch* updates
            if gazeCursorSearch.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                gazeCursorSearch.frameNStart = frameN  # exact frame index
                gazeCursorSearch.tStart = t  # local t and not account for scr refresh
                gazeCursorSearch.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(gazeCursorSearch, 'tStartRefresh')  # time at next scr refresh
                gazeCursorSearch.setAutoDraw(True)
            if gazeCursorSearch.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > gazeCursorSearch.tStartRefresh + 30-frameTolerance:
                    # keep track of stop time/frame for later
                    gazeCursorSearch.tStop = t  # not accounting for scr refresh
                    gazeCursorSearch.frameNStop = frameN  # exact frame index
                    gazeCursorSearch.setAutoDraw(False)
            if gazeCursorSearch.status == STARTED:  # only update if drawing
                gazeCursorSearch.setOpacity(1.0, log=False)
                gazeCursorSearch.setPos([eyetracker.getPos()], log=False)
            
            # check for quit (typically the Esc key)
            if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
                core.quit()
            
            # check if all components have finished
            if not continueRoutine:  # a component has requested a forced-end of Routine
                routineForceEnded = True
                break
            continueRoutine = False  # will revert to True if at least one component still running
            for thisComponent in search_taskComponents:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # --- Ending Routine "search_task" ---
        for thisComponent in search_taskComponents:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        # Run 'End Routine' code from codeSearchTask
        for target in myTargets:
            target.setAutoDraw(False)
        
        for distractor in myDistractors:
            distractor.setAutoDraw(False)
        # using non-slip timing so subtract the expected duration of this Routine (unless ended on request)
        if routineForceEnded:
            routineTimer.reset()
        else:
            routineTimer.addTime(-30.000000)
        
        # --- Prepare to start Routine "crossEnd" ---
        continueRoutine = True
        routineForceEnded = False
        # update component parameters for each repeat
        # Run 'Begin Routine' code from code_end
        eyetracker.setRecordingState(False)
        jittered_duration_cross = 3.0 + random()
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
                fixcrossE.setAutoDraw(True)
            if fixcrossE.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > fixcrossE.tStartRefresh + jittered_duration_cross-frameTolerance:
                    # keep track of stop time/frame for later
                    fixcrossE.tStop = t  # not accounting for scr refresh
                    fixcrossE.frameNStop = frameN  # exact frame index
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
        
    # completed 2.0 repeats of 'search_task_trials'
    
    
    # --- Prepare to start Routine "startTask" ---
    continueRoutine = True
    routineForceEnded = False
    # update component parameters for each repeat
    # Run 'Begin Routine' code from codeStartTask
    print("Start Task Screen. Press Space to continue.")
    win.mouseVisible = False
    spaceStartTask.keys = []
    spaceStartTask.rt = []
    _spaceStartTask_allKeys = []
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
    # the Routine "startTask" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset()
    
    # --- Prepare to start Routine "blank" ---
    continueRoutine = True
    routineForceEnded = False
    # update component parameters for each repeat
    # Run 'Begin Routine' code from codeBlank
    jittered_duration_blank = 1.0 + random()
    blankScreen.setText('')
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
            blankScreen.setAutoDraw(True)
        if blankScreen.status == STARTED:
            # is it time to stop? (based on global clock, using actual start)
            if tThisFlipGlobal > blankScreen.tStartRefresh + jittered_duration_blank-frameTolerance:
                # keep track of stop time/frame for later
                blankScreen.tStop = t  # not accounting for scr refresh
                blankScreen.frameNStop = frameN  # exact frame index
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
    # Run 'End Routine' code from codeBlank
    # print(f"Start: {round(blankScreen.tStart, 3)}, End: {round(blankScreen.tStop, 3)}, Duration: {round(blankScreen.tStop - blankScreen.tStart, 3)}")
    # the Routine "blank" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset()
    
    # set up handler to look after randomisation of conditions etc
    learning_trials = data.TrialHandler(nReps=2.0, method='random', 
        extraInfo=expInfo, originPath=-1,
        trialList=data.importConditions(posFile),
        seed=None, name='learning_trials')
    thisExp.addLoop(learning_trials)  # add the loop to the experiment
    thisLearning_trial = learning_trials.trialList[0]  # so we can initialise stimuli with some values
    # abbreviate parameter names if possible (e.g. rgb = thisLearning_trial.rgb)
    if thisLearning_trial != None:
        for paramName in thisLearning_trial:
            exec('{} = thisLearning_trial[paramName]'.format(paramName))
    
    for thisLearning_trial in learning_trials:
        currentLoop = learning_trials
        # abbreviate parameter names if possible (e.g. rgb = thisLearning_trial.rgb)
        if thisLearning_trial != None:
            for paramName in thisLearning_trial:
                exec('{} = thisLearning_trial[paramName]'.format(paramName))
        
        # --- Prepare to start Routine "cross" ---
        continueRoutine = True
        routineForceEnded = False
        # update component parameters for each repeat
        # Run 'Begin Routine' code from codeFixate
        logging.log(level=logging.INFO, msg=f'FixationCross')
        eyetracker.setRecordingState(True)
         
        if currentLoop.name == "learning_trials":
            eyetracker.sendMessage("Trial " + str(learning_trials.thisN+1) + ", Condition " + stim)
            print("VP %s TRIALID %d CONDITION %s"%(expInfo['participant'], learning_trials.thisN+1, stim))
        
        if currentLoop.name == "trials":
            eyetracker.sendMessage("Trial " + str(trials.thisN+learning_trials.thisN+1) + ", Condition " + stim)
            print("VP %s TRIALID %d CONDITION %s"%(expInfo['participant'], trials.thisN+learning_trials.thisN+1, stim))
        
        if currentLoop.name == "testtrials":
            eyetracker.sendMessage("Test-Trial " + str(testtrials.thisN+1))
            print("VP %s TESTTRIALID %d"%(expInfo['participant'], testtrials.thisN+1))
            
        if currentLoop.name == "testtrials_novelty":
            eyetracker.sendMessage("Test-Trial " + str(testtrials_novelty.thisN+testtrials.thisN+1))
            print("VP %s TESTTRIALID %d"%(expInfo['participant'], testtrials_novelty.thisN+testtrials.thisN+1))
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
                fixcross.setAutoDraw(True)
            if fixcross.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > fixcross.tStartRefresh + 1-frameTolerance:
                    # keep track of stop time/frame for later
                    fixcross.tStop = t  # not accounting for scr refresh
                    fixcross.frameNStop = frameN  # exact frame index
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
        # Run 'Begin Routine' code from codeBlank
        jittered_duration_blank = 1.0 + random()
        blankScreen.setText('')
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
                blankScreen.setAutoDraw(True)
            if blankScreen.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > blankScreen.tStartRefresh + jittered_duration_blank-frameTolerance:
                    # keep track of stop time/frame for later
                    blankScreen.tStop = t  # not accounting for scr refresh
                    blankScreen.frameNStop = frameN  # exact frame index
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
        # Run 'End Routine' code from codeBlank
        # print(f"Start: {round(blankScreen.tStart, 3)}, End: {round(blankScreen.tStop, 3)}, Duration: {round(blankScreen.tStop - blankScreen.tStart, 3)}")
        # the Routine "blank" was not non-slip safe, so reset the non-slip timer
        routineTimer.reset()
        
        # --- Prepare to start Routine "trial" ---
        continueRoutine = True
        routineForceEnded = False
        # update component parameters for each repeat
        # Run 'Begin Routine' code from codeTrial
        looked_at = False
        cursorcolor="white"
        
        logging.log(level=logging.INFO, msg=f'ImageOnset_{stim}')
        ioServer.sendMessageEvent(text='ImageOnset')
        eyetracker.sendMessage('ImageOnset')
        image.setPos(position)
        image.setImage(eval(stim))
        roi.setPos(position)
        roi.setSize([image.size])
        # clear any previous roi data
        roi.reset()
        # keep track of which components have finished
        trialComponents = [image, roi, gazeCursor]
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
        while continueRoutine and routineTimer.getTime() < 3.0:
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
                if tThisFlipGlobal > image.tStartRefresh + 3-frameTolerance:
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
                if tThisFlipGlobal > roi.tStartRefresh + 3-frameTolerance:
                    # keep track of stop time/frame for later
                    roi.tStop = t  # not accounting for scr refresh
                    roi.frameNStop = frameN  # exact frame index
                    roi.status = FINISHED
            
            # *gazeCursor* updates
            if gazeCursor.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                gazeCursor.frameNStart = frameN  # exact frame index
                gazeCursor.tStart = t  # local t and not account for scr refresh
                gazeCursor.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(gazeCursor, 'tStartRefresh')  # time at next scr refresh
                gazeCursor.setAutoDraw(True)
            if gazeCursor.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > gazeCursor.tStartRefresh + 3-frameTolerance:
                    # keep track of stop time/frame for later
                    gazeCursor.tStop = t  # not accounting for scr refresh
                    gazeCursor.frameNStop = frameN  # exact frame index
                    gazeCursor.setAutoDraw(False)
            if gazeCursor.status == STARTED:  # only update if drawing
                gazeCursor.setFillColor(cursorcolor, log=False)
                gazeCursor.setOpacity(0.0, log=False)
                gazeCursor.setPos([eyetracker.getPos()], log=False)
            
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
        if looked_at & ("neg" in stim):
            feedback_color = "red"
            feedback_opacity = 1
            cursorcolor="red"
            feedback_audio = "audio/error.wav"
            log_msg = "lookedat_cs_plus"
            points = -3
            score += points
            feedback_points = f"{points}"
            feedback_score = f"Score: {score}"
            outcome = "loss"
        elif looked_at & ("pos" in stim):
            feedback_color = "green"
            feedback_opacity = 1
            cursorcolor="green"
            feedback_audio = "audio/win.wav"
            log_msg = "lookedat_cs_minus"
            points = 3
            score += points
            feedback_points = f"+{points}"
            feedback_score = f"Score: {score}"
            outcome = "reward"
        else:
            feedback_color = (0,0,0)
            feedback_opacity = 0
            feedback_audio = "audio/silence.wav"
            log_msg = "no_look"
            points = 0
            feedback_points = ""
            feedback_score = ""
            outcome = "none"
        
        rectsize = [item * 1.02 for item in image.size]
        
        image_w = image.size[0]
        image_h = image.size[1]
        imagesize_test = [image_w * 1.2, image_h * 1.2]
        
        thisExp.addData('score', score)
        thisExp.addData('outcome', outcome)
        
        learning_trials.addData('roi.numLooks', roi.numLooks)
        if roi.numLooks:
           learning_trials.addData('roi.timesOn', roi.timesOn[0])
           learning_trials.addData('roi.timesOff', roi.timesOff[0])
        else:
           learning_trials.addData('roi.timesOn', "")
           learning_trials.addData('roi.timesOff', "")
        # using non-slip timing so subtract the expected duration of this Routine (unless ended on request)
        if routineForceEnded:
            routineTimer.reset()
        else:
            routineTimer.addTime(-3.000000)
        
        # --- Prepare to start Routine "feedback" ---
        continueRoutine = True
        routineForceEnded = False
        # update component parameters for each repeat
        # Run 'Begin Routine' code from codeFeedback
        logging.log(level=logging.INFO, msg=f'FeedbackOnset_{log_msg}')
        ioServer.sendMessageEvent(text='FeedbackOnset')
        eyetracker.sendMessage('FeedbackOnset')
        eyetracker.sendMessage(outcome)
        
        textPoints.setColor(feedback_color, colorSpace='rgb')
        polygon.setOpacity(feedback_opacity)
        polygon.setPos(position)
        polygon.setLineColor(feedback_color)
        imageFeedback.setPos(position)
        imageFeedback.setImage(eval(stim))
        soundFeedback.setSound(feedback_audio, secs=1, hamming=True)
        soundFeedback.setVolume(1.0, log=False)
        # keep track of which components have finished
        feedbackComponents = [textPoints, textScore, polygon, imageFeedback, gazeCursor_Feedback, soundFeedback]
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
        while continueRoutine and routineTimer.getTime() < 1.0:
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
                if tThisFlipGlobal > textPoints.tStartRefresh + 1-frameTolerance:
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
                if tThisFlipGlobal > textScore.tStartRefresh + 1-frameTolerance:
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
                if tThisFlipGlobal > polygon.tStartRefresh + 1-frameTolerance:
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
                if tThisFlipGlobal > imageFeedback.tStartRefresh + 1-frameTolerance:
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
                gazeCursor_Feedback.setAutoDraw(True)
            if gazeCursor_Feedback.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > gazeCursor_Feedback.tStartRefresh + 1-frameTolerance:
                    # keep track of stop time/frame for later
                    gazeCursor_Feedback.tStop = t  # not accounting for scr refresh
                    gazeCursor_Feedback.frameNStop = frameN  # exact frame index
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
                soundFeedback.play()  # start the sound (it finishes automatically)
            if soundFeedback.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > soundFeedback.tStartRefresh + 1-frameTolerance:
                    # keep track of stop time/frame for later
                    soundFeedback.tStop = t  # not accounting for scr refresh
                    soundFeedback.frameNStop = frameN  # exact frame index
                    soundFeedback.stop()
            
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
        # using non-slip timing so subtract the expected duration of this Routine (unless ended on request)
        if routineForceEnded:
            routineTimer.reset()
        else:
            routineTimer.addTime(-1.000000)
        
        # --- Prepare to start Routine "crossEnd" ---
        continueRoutine = True
        routineForceEnded = False
        # update component parameters for each repeat
        # Run 'Begin Routine' code from code_end
        eyetracker.setRecordingState(False)
        jittered_duration_cross = 3.0 + random()
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
                fixcrossE.setAutoDraw(True)
            if fixcrossE.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > fixcrossE.tStartRefresh + jittered_duration_cross-frameTolerance:
                    # keep track of stop time/frame for later
                    fixcrossE.tStop = t  # not accounting for scr refresh
                    fixcrossE.frameNStop = frameN  # exact frame index
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
        
    # completed 2.0 repeats of 'learning_trials'
    
    
    # --- Prepare to start Routine "repeatTask" ---
    continueRoutine = True
    routineForceEnded = False
    # update component parameters for each repeat
    # Run 'Begin Routine' code from codeRepeatTask
    print("Second Instruction Screen. Press Space to continue.")
    win.mouseVisible = False
    textboxRepeatInstruction.reset()
    spaceRepeatTask.keys = []
    spaceRepeatTask.rt = []
    _spaceRepeatTask_allKeys = []
    # keep track of which components have finished
    repeatTaskComponents = [textRepeatTask, textboxRepeatInstruction, textRepeatSpace, spaceRepeatTask]
    for thisComponent in repeatTaskComponents:
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
    
    # --- Run Routine "repeatTask" ---
    while continueRoutine:
        # get current time
        t = routineTimer.getTime()
        tThisFlip = win.getFutureFlipTime(clock=routineTimer)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        
        # *textRepeatTask* updates
        if textRepeatTask.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            textRepeatTask.frameNStart = frameN  # exact frame index
            textRepeatTask.tStart = t  # local t and not account for scr refresh
            textRepeatTask.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(textRepeatTask, 'tStartRefresh')  # time at next scr refresh
            textRepeatTask.setAutoDraw(True)
        
        # *textboxRepeatInstruction* updates
        if textboxRepeatInstruction.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            textboxRepeatInstruction.frameNStart = frameN  # exact frame index
            textboxRepeatInstruction.tStart = t  # local t and not account for scr refresh
            textboxRepeatInstruction.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(textboxRepeatInstruction, 'tStartRefresh')  # time at next scr refresh
            textboxRepeatInstruction.setAutoDraw(True)
        
        # *textRepeatSpace* updates
        if textRepeatSpace.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            textRepeatSpace.frameNStart = frameN  # exact frame index
            textRepeatSpace.tStart = t  # local t and not account for scr refresh
            textRepeatSpace.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(textRepeatSpace, 'tStartRefresh')  # time at next scr refresh
            textRepeatSpace.setAutoDraw(True)
        
        # *spaceRepeatTask* updates
        waitOnFlip = False
        if spaceRepeatTask.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            spaceRepeatTask.frameNStart = frameN  # exact frame index
            spaceRepeatTask.tStart = t  # local t and not account for scr refresh
            spaceRepeatTask.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(spaceRepeatTask, 'tStartRefresh')  # time at next scr refresh
            spaceRepeatTask.status = STARTED
            # keyboard checking is just starting
            waitOnFlip = True
            win.callOnFlip(spaceRepeatTask.clock.reset)  # t=0 on next screen flip
            win.callOnFlip(spaceRepeatTask.clearEvents, eventType='keyboard')  # clear events on next screen flip
        if spaceRepeatTask.status == STARTED and not waitOnFlip:
            theseKeys = spaceRepeatTask.getKeys(keyList=['space', 'enter'], waitRelease=False)
            _spaceRepeatTask_allKeys.extend(theseKeys)
            if len(_spaceRepeatTask_allKeys):
                spaceRepeatTask.keys = _spaceRepeatTask_allKeys[-1].name  # just the last key pressed
                spaceRepeatTask.rt = _spaceRepeatTask_allKeys[-1].rt
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
        for thisComponent in repeatTaskComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # --- Ending Routine "repeatTask" ---
    for thisComponent in repeatTaskComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    # the Routine "repeatTask" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset()
    
    # --- Prepare to start Routine "blank" ---
    continueRoutine = True
    routineForceEnded = False
    # update component parameters for each repeat
    # Run 'Begin Routine' code from codeBlank
    jittered_duration_blank = 1.0 + random()
    blankScreen.setText('')
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
            blankScreen.setAutoDraw(True)
        if blankScreen.status == STARTED:
            # is it time to stop? (based on global clock, using actual start)
            if tThisFlipGlobal > blankScreen.tStartRefresh + jittered_duration_blank-frameTolerance:
                # keep track of stop time/frame for later
                blankScreen.tStop = t  # not accounting for scr refresh
                blankScreen.frameNStop = frameN  # exact frame index
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
    # Run 'End Routine' code from codeBlank
    # print(f"Start: {round(blankScreen.tStart, 3)}, End: {round(blankScreen.tStop, 3)}, Duration: {round(blankScreen.tStop - blankScreen.tStart, 3)}")
    # the Routine "blank" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset()
    
    # set up handler to look after randomisation of conditions etc
    trials = data.TrialHandler(nReps=2.0, method='random', 
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
         
        if currentLoop.name == "learning_trials":
            eyetracker.sendMessage("Trial " + str(learning_trials.thisN+1) + ", Condition " + stim)
            print("VP %s TRIALID %d CONDITION %s"%(expInfo['participant'], learning_trials.thisN+1, stim))
        
        if currentLoop.name == "trials":
            eyetracker.sendMessage("Trial " + str(trials.thisN+learning_trials.thisN+1) + ", Condition " + stim)
            print("VP %s TRIALID %d CONDITION %s"%(expInfo['participant'], trials.thisN+learning_trials.thisN+1, stim))
        
        if currentLoop.name == "testtrials":
            eyetracker.sendMessage("Test-Trial " + str(testtrials.thisN+1))
            print("VP %s TESTTRIALID %d"%(expInfo['participant'], testtrials.thisN+1))
            
        if currentLoop.name == "testtrials_novelty":
            eyetracker.sendMessage("Test-Trial " + str(testtrials_novelty.thisN+testtrials.thisN+1))
            print("VP %s TESTTRIALID %d"%(expInfo['participant'], testtrials_novelty.thisN+testtrials.thisN+1))
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
                fixcross.setAutoDraw(True)
            if fixcross.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > fixcross.tStartRefresh + 1-frameTolerance:
                    # keep track of stop time/frame for later
                    fixcross.tStop = t  # not accounting for scr refresh
                    fixcross.frameNStop = frameN  # exact frame index
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
        # Run 'Begin Routine' code from codeBlank
        jittered_duration_blank = 1.0 + random()
        blankScreen.setText('')
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
                blankScreen.setAutoDraw(True)
            if blankScreen.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > blankScreen.tStartRefresh + jittered_duration_blank-frameTolerance:
                    # keep track of stop time/frame for later
                    blankScreen.tStop = t  # not accounting for scr refresh
                    blankScreen.frameNStop = frameN  # exact frame index
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
        # Run 'End Routine' code from codeBlank
        # print(f"Start: {round(blankScreen.tStart, 3)}, End: {round(blankScreen.tStop, 3)}, Duration: {round(blankScreen.tStop - blankScreen.tStart, 3)}")
        # the Routine "blank" was not non-slip safe, so reset the non-slip timer
        routineTimer.reset()
        
        # --- Prepare to start Routine "trial" ---
        continueRoutine = True
        routineForceEnded = False
        # update component parameters for each repeat
        # Run 'Begin Routine' code from codeTrial
        looked_at = False
        cursorcolor="white"
        
        logging.log(level=logging.INFO, msg=f'ImageOnset_{stim}')
        ioServer.sendMessageEvent(text='ImageOnset')
        eyetracker.sendMessage('ImageOnset')
        image.setPos(position)
        image.setImage(eval(stim))
        roi.setPos(position)
        roi.setSize([image.size])
        # clear any previous roi data
        roi.reset()
        # keep track of which components have finished
        trialComponents = [image, roi, gazeCursor]
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
        while continueRoutine and routineTimer.getTime() < 3.0:
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
                if tThisFlipGlobal > image.tStartRefresh + 3-frameTolerance:
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
                if tThisFlipGlobal > roi.tStartRefresh + 3-frameTolerance:
                    # keep track of stop time/frame for later
                    roi.tStop = t  # not accounting for scr refresh
                    roi.frameNStop = frameN  # exact frame index
                    roi.status = FINISHED
            
            # *gazeCursor* updates
            if gazeCursor.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                gazeCursor.frameNStart = frameN  # exact frame index
                gazeCursor.tStart = t  # local t and not account for scr refresh
                gazeCursor.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(gazeCursor, 'tStartRefresh')  # time at next scr refresh
                gazeCursor.setAutoDraw(True)
            if gazeCursor.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > gazeCursor.tStartRefresh + 3-frameTolerance:
                    # keep track of stop time/frame for later
                    gazeCursor.tStop = t  # not accounting for scr refresh
                    gazeCursor.frameNStop = frameN  # exact frame index
                    gazeCursor.setAutoDraw(False)
            if gazeCursor.status == STARTED:  # only update if drawing
                gazeCursor.setFillColor(cursorcolor, log=False)
                gazeCursor.setOpacity(0.0, log=False)
                gazeCursor.setPos([eyetracker.getPos()], log=False)
            
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
        if looked_at & ("neg" in stim):
            feedback_color = "red"
            feedback_opacity = 1
            cursorcolor="red"
            feedback_audio = "audio/error.wav"
            log_msg = "lookedat_cs_plus"
            points = -3
            score += points
            feedback_points = f"{points}"
            feedback_score = f"Score: {score}"
            outcome = "loss"
        elif looked_at & ("pos" in stim):
            feedback_color = "green"
            feedback_opacity = 1
            cursorcolor="green"
            feedback_audio = "audio/win.wav"
            log_msg = "lookedat_cs_minus"
            points = 3
            score += points
            feedback_points = f"+{points}"
            feedback_score = f"Score: {score}"
            outcome = "reward"
        else:
            feedback_color = (0,0,0)
            feedback_opacity = 0
            feedback_audio = "audio/silence.wav"
            log_msg = "no_look"
            points = 0
            feedback_points = ""
            feedback_score = ""
            outcome = "none"
        
        rectsize = [item * 1.02 for item in image.size]
        
        image_w = image.size[0]
        image_h = image.size[1]
        imagesize_test = [image_w * 1.2, image_h * 1.2]
        
        thisExp.addData('score', score)
        thisExp.addData('outcome', outcome)
        
        trials.addData('roi.numLooks', roi.numLooks)
        if roi.numLooks:
           trials.addData('roi.timesOn', roi.timesOn[0])
           trials.addData('roi.timesOff', roi.timesOff[0])
        else:
           trials.addData('roi.timesOn', "")
           trials.addData('roi.timesOff', "")
        # using non-slip timing so subtract the expected duration of this Routine (unless ended on request)
        if routineForceEnded:
            routineTimer.reset()
        else:
            routineTimer.addTime(-3.000000)
        
        # --- Prepare to start Routine "feedback" ---
        continueRoutine = True
        routineForceEnded = False
        # update component parameters for each repeat
        # Run 'Begin Routine' code from codeFeedback
        logging.log(level=logging.INFO, msg=f'FeedbackOnset_{log_msg}')
        ioServer.sendMessageEvent(text='FeedbackOnset')
        eyetracker.sendMessage('FeedbackOnset')
        eyetracker.sendMessage(outcome)
        
        textPoints.setColor(feedback_color, colorSpace='rgb')
        polygon.setOpacity(feedback_opacity)
        polygon.setPos(position)
        polygon.setLineColor(feedback_color)
        imageFeedback.setPos(position)
        imageFeedback.setImage(eval(stim))
        soundFeedback.setSound(feedback_audio, secs=1, hamming=True)
        soundFeedback.setVolume(1.0, log=False)
        # keep track of which components have finished
        feedbackComponents = [textPoints, textScore, polygon, imageFeedback, gazeCursor_Feedback, soundFeedback]
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
        while continueRoutine and routineTimer.getTime() < 1.0:
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
                if tThisFlipGlobal > textPoints.tStartRefresh + 1-frameTolerance:
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
                if tThisFlipGlobal > textScore.tStartRefresh + 1-frameTolerance:
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
                if tThisFlipGlobal > polygon.tStartRefresh + 1-frameTolerance:
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
                if tThisFlipGlobal > imageFeedback.tStartRefresh + 1-frameTolerance:
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
                gazeCursor_Feedback.setAutoDraw(True)
            if gazeCursor_Feedback.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > gazeCursor_Feedback.tStartRefresh + 1-frameTolerance:
                    # keep track of stop time/frame for later
                    gazeCursor_Feedback.tStop = t  # not accounting for scr refresh
                    gazeCursor_Feedback.frameNStop = frameN  # exact frame index
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
                soundFeedback.play()  # start the sound (it finishes automatically)
            if soundFeedback.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > soundFeedback.tStartRefresh + 1-frameTolerance:
                    # keep track of stop time/frame for later
                    soundFeedback.tStop = t  # not accounting for scr refresh
                    soundFeedback.frameNStop = frameN  # exact frame index
                    soundFeedback.stop()
            
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
        # using non-slip timing so subtract the expected duration of this Routine (unless ended on request)
        if routineForceEnded:
            routineTimer.reset()
        else:
            routineTimer.addTime(-1.000000)
        
        # --- Prepare to start Routine "crossEnd" ---
        continueRoutine = True
        routineForceEnded = False
        # update component parameters for each repeat
        # Run 'Begin Routine' code from code_end
        eyetracker.setRecordingState(False)
        jittered_duration_cross = 3.0 + random()
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
                fixcrossE.setAutoDraw(True)
            if fixcrossE.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > fixcrossE.tStartRefresh + jittered_duration_cross-frameTolerance:
                    # keep track of stop time/frame for later
                    fixcrossE.tStop = t  # not accounting for scr refresh
                    fixcrossE.frameNStop = frameN  # exact frame index
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
        
    # completed 2.0 repeats of 'trials'
    
    
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
        logging.log(level=logging.INFO, msg=f'StimRating_{stim}')
        print("Rating of Stimulus: %s"%(stim))
        win.mouseVisible = True
        imageRating.setPos((0, 0.2))
        imageRating.setImage(eval(stim))
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
    # the Routine "testInstr" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset()
    
    # set up handler to look after randomisation of conditions etc
    testtrials = data.TrialHandler(nReps=2.0, method='random', 
        extraInfo=expInfo, originPath=-1,
        trialList=data.importConditions(posFileTest),
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
         
        if currentLoop.name == "learning_trials":
            eyetracker.sendMessage("Trial " + str(learning_trials.thisN+1) + ", Condition " + stim)
            print("VP %s TRIALID %d CONDITION %s"%(expInfo['participant'], learning_trials.thisN+1, stim))
        
        if currentLoop.name == "trials":
            eyetracker.sendMessage("Trial " + str(trials.thisN+learning_trials.thisN+1) + ", Condition " + stim)
            print("VP %s TRIALID %d CONDITION %s"%(expInfo['participant'], trials.thisN+learning_trials.thisN+1, stim))
        
        if currentLoop.name == "testtrials":
            eyetracker.sendMessage("Test-Trial " + str(testtrials.thisN+1))
            print("VP %s TESTTRIALID %d"%(expInfo['participant'], testtrials.thisN+1))
            
        if currentLoop.name == "testtrials_novelty":
            eyetracker.sendMessage("Test-Trial " + str(testtrials_novelty.thisN+testtrials.thisN+1))
            print("VP %s TESTTRIALID %d"%(expInfo['participant'], testtrials_novelty.thisN+testtrials.thisN+1))
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
                fixcross.setAutoDraw(True)
            if fixcross.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > fixcross.tStartRefresh + 1-frameTolerance:
                    # keep track of stop time/frame for later
                    fixcross.tStop = t  # not accounting for scr refresh
                    fixcross.frameNStop = frameN  # exact frame index
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
        
        # --- Prepare to start Routine "blank_test" ---
        continueRoutine = True
        routineForceEnded = False
        # update component parameters for each repeat
        # Run 'Begin Routine' code from codeBlankTest
        jittered_duration_blank = 1.0 + random()
        looked_at_roi1_trial = False
        looked_at_roi2_trial = False
        looked_at_roi3_trial = False
        looked_at_roi4_trial = False
        blankScreenTest.setText('')
        # keep track of which components have finished
        blank_testComponents = [blankScreenTest]
        for thisComponent in blank_testComponents:
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
        
        # --- Run Routine "blank_test" ---
        while continueRoutine:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # *blankScreenTest* updates
            if blankScreenTest.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                blankScreenTest.frameNStart = frameN  # exact frame index
                blankScreenTest.tStart = t  # local t and not account for scr refresh
                blankScreenTest.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(blankScreenTest, 'tStartRefresh')  # time at next scr refresh
                blankScreenTest.setAutoDraw(True)
            if blankScreenTest.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > blankScreenTest.tStartRefresh + jittered_duration_blank-frameTolerance:
                    # keep track of stop time/frame for later
                    blankScreenTest.tStop = t  # not accounting for scr refresh
                    blankScreenTest.frameNStop = frameN  # exact frame index
                    blankScreenTest.setAutoDraw(False)
            
            # check for quit (typically the Esc key)
            if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
                core.quit()
            
            # check if all components have finished
            if not continueRoutine:  # a component has requested a forced-end of Routine
                routineForceEnded = True
                break
            continueRoutine = False  # will revert to True if at least one component still running
            for thisComponent in blank_testComponents:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # --- Ending Routine "blank_test" ---
        for thisComponent in blank_testComponents:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        # Run 'End Routine' code from codeBlankTest
        # print(f"Start: {round(blankScreen.tStart, 3)}, End: {round(blankScreen.tStop, 3)}, Duration: {round(blankScreen.tStop - blankScreen.tStart, 3)}")
        loop_timer = core.Clock()
        # the Routine "blank_test" was not non-slip safe, so reset the non-slip timer
        routineTimer.reset()
        
        # set up handler to look after randomisation of conditions etc
        loop_test = data.TrialHandler(nReps=100.0, method='sequential', 
            extraInfo=expInfo, originPath=-1,
            trialList=[None],
            seed=None, name='loop_test')
        thisExp.addLoop(loop_test)  # add the loop to the experiment
        thisLoop_test = loop_test.trialList[0]  # so we can initialise stimuli with some values
        # abbreviate parameter names if possible (e.g. rgb = thisLoop_test.rgb)
        if thisLoop_test != None:
            for paramName in thisLoop_test:
                exec('{} = thisLoop_test[paramName]'.format(paramName))
        
        for thisLoop_test in loop_test:
            currentLoop = loop_test
            # abbreviate parameter names if possible (e.g. rgb = thisLoop_test.rgb)
            if thisLoop_test != None:
                for paramName in thisLoop_test:
                    exec('{} = thisLoop_test[paramName]'.format(paramName))
            
            # --- Prepare to start Routine "test_trial" ---
            continueRoutine = True
            routineForceEnded = False
            # update component parameters for each repeat
            # Run 'Begin Routine' code from codeTestTrial
            looked_at_roi1 = False
            looked_at_roi2 = False
            looked_at_roi3 = False
            looked_at_roi4 = False
            cursorcolor="white"
            
            logging.log(level=logging.INFO, msg=f'ImageTestOnset')
            ioServer.sendMessageEvent(text='ImageTestOnset')
            eyetracker.sendMessage('ImageTestOnset')
            test_image1.setPos(position1)
            test_image1.setImage(eval(stim1))
            test_roi1.setPos(position1)
            test_roi1.setSize([test_image1.size])
            # clear any previous roi data
            test_roi1.reset()
            test_image2.setPos(position2)
            test_image2.setImage(eval(stim2))
            test_roi2.setPos(position2)
            test_roi2.setSize([test_image2.size])
            # clear any previous roi data
            test_roi2.reset()
            test_image3.setPos(position3)
            test_image3.setImage(eval(stim3))
            test_roi3.setPos(position3)
            test_roi3.setSize([test_image3.size])
            # clear any previous roi data
            test_roi3.reset()
            test_image4.setPos(position4)
            test_image4.setImage(eval(stim4))
            test_roi4.setPos(position4)
            test_roi4.setSize([test_image4.size])
            # clear any previous roi data
            test_roi4.reset()
            # keep track of which components have finished
            test_trialComponents = [test_image1, test_roi1, test_image2, test_roi2, test_image3, test_roi3, test_image4, test_roi4, gazeCursorTest]
            for thisComponent in test_trialComponents:
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
            
            # --- Run Routine "test_trial" ---
            while continueRoutine and routineTimer.getTime() < 2.0:
                # get current time
                t = routineTimer.getTime()
                tThisFlip = win.getFutureFlipTime(clock=routineTimer)
                tThisFlipGlobal = win.getFutureFlipTime(clock=None)
                frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
                # update/draw components on each frame
                # Run 'Each Frame' code from codeTestTrial
                if loop_timer.getTime() >= 2.0:
                    continueRoutine = False
                    currentLoop.finished = True
                
                test_roi1.size = test_image1.size
                test_roi1.pos = position1
                
                if test_roi1.isLookedIn and not looked_at_roi1_trial:
                    looked_at_roi1 = True
                    looked_at_roi1_trial = True
                    continueRoutine = False
                
                test_roi2.size = test_image2.size
                test_roi2.pos = position2
                
                if test_roi2.isLookedIn and not looked_at_roi2_trial:
                    looked_at_roi2 = True
                    looked_at_roi2_trial = True
                    continueRoutine = False
                    
                test_roi3.size = test_image3.size
                test_roi3.pos = position3
                
                if test_roi3.isLookedIn and not looked_at_roi3_trial:
                    looked_at_roi3 = True
                    looked_at_roi3_trial = True
                    continueRoutine = False
                   
                test_roi4.size = test_image4.size
                test_roi4.pos = position4
                
                if test_roi4.isLookedIn and not looked_at_roi4_trial:
                    looked_at_roi4 = True
                    looked_at_roi4_trial = True
                    continueRoutine = False
                    
                if event.getKeys(keyList=["p"]) and not paused:
                    print("Experiment Paused - Press 'p' to continue.")
                    paused = True
                    event.waitKeys(keyList=["p"], clearEvents=True)
                    
                    print("Experiment Continued")
                    paused = False
                    event.clearEvents(eventType="keyboard")
                    
                
                # *test_image1* updates
                if test_image1.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                    # keep track of start time/frame for later
                    test_image1.frameNStart = frameN  # exact frame index
                    test_image1.tStart = t  # local t and not account for scr refresh
                    test_image1.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(test_image1, 'tStartRefresh')  # time at next scr refresh
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'test_image1.started')
                    test_image1.setAutoDraw(True)
                if test_image1.status == STARTED:
                    # is it time to stop? (based on global clock, using actual start)
                    if tThisFlipGlobal > test_image1.tStartRefresh + 2-frameTolerance:
                        # keep track of stop time/frame for later
                        test_image1.tStop = t  # not accounting for scr refresh
                        test_image1.frameNStop = frameN  # exact frame index
                        # add timestamp to datafile
                        thisExp.timestampOnFlip(win, 'test_image1.stopped')
                        test_image1.setAutoDraw(False)
                if test_roi1.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                    # keep track of start time/frame for later
                    test_roi1.frameNStart = frameN  # exact frame index
                    test_roi1.tStart = t  # local t and not account for scr refresh
                    test_roi1.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(test_roi1, 'tStartRefresh')  # time at next scr refresh
                    test_roi1.status = STARTED
                if test_roi1.status == STARTED:
                    # check whether test_roi1 has been looked in
                    if test_roi1.isLookedIn:
                        if not test_roi1.wasLookedIn:
                            test_roi1.timesOn.append(routineTimer.getTime()) # store time of first look
                            test_roi1.timesOff.append(routineTimer.getTime()) # store time looked until
                        else:
                            test_roi1.timesOff[-1] = routineTimer.getTime() # update time looked until
                        test_roi1.wasLookedIn = True  # if test_roi1 is still looked at next frame, it is not a new look
                    else:
                        if test_roi1.wasLookedIn:
                            test_roi1.timesOff[-1] = routineTimer.getTime() # update time looked until
                        test_roi1.wasLookedIn = False  # if test_roi1 is looked at next frame, it is a new look
                else:
                    test_roi1.clock.reset() # keep clock at 0 if roi hasn't started / has finished
                    test_roi1.wasLookedIn = False  # if test_roi1 is looked at next frame, it is a new look
                if test_roi1.status == STARTED:
                    # is it time to stop? (based on global clock, using actual start)
                    if tThisFlipGlobal > test_roi1.tStartRefresh + 2-frameTolerance:
                        # keep track of stop time/frame for later
                        test_roi1.tStop = t  # not accounting for scr refresh
                        test_roi1.frameNStop = frameN  # exact frame index
                        test_roi1.status = FINISHED
                
                # *test_image2* updates
                if test_image2.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                    # keep track of start time/frame for later
                    test_image2.frameNStart = frameN  # exact frame index
                    test_image2.tStart = t  # local t and not account for scr refresh
                    test_image2.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(test_image2, 'tStartRefresh')  # time at next scr refresh
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'test_image2.started')
                    test_image2.setAutoDraw(True)
                if test_image2.status == STARTED:
                    # is it time to stop? (based on global clock, using actual start)
                    if tThisFlipGlobal > test_image2.tStartRefresh + 2-frameTolerance:
                        # keep track of stop time/frame for later
                        test_image2.tStop = t  # not accounting for scr refresh
                        test_image2.frameNStop = frameN  # exact frame index
                        # add timestamp to datafile
                        thisExp.timestampOnFlip(win, 'test_image2.stopped')
                        test_image2.setAutoDraw(False)
                if test_roi2.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                    # keep track of start time/frame for later
                    test_roi2.frameNStart = frameN  # exact frame index
                    test_roi2.tStart = t  # local t and not account for scr refresh
                    test_roi2.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(test_roi2, 'tStartRefresh')  # time at next scr refresh
                    test_roi2.status = STARTED
                if test_roi2.status == STARTED:
                    # check whether test_roi2 has been looked in
                    if test_roi2.isLookedIn:
                        if not test_roi2.wasLookedIn:
                            test_roi2.timesOn.append(routineTimer.getTime()) # store time of first look
                            test_roi2.timesOff.append(routineTimer.getTime()) # store time looked until
                        else:
                            test_roi2.timesOff[-1] = routineTimer.getTime() # update time looked until
                        test_roi2.wasLookedIn = True  # if test_roi2 is still looked at next frame, it is not a new look
                    else:
                        if test_roi2.wasLookedIn:
                            test_roi2.timesOff[-1] = routineTimer.getTime() # update time looked until
                        test_roi2.wasLookedIn = False  # if test_roi2 is looked at next frame, it is a new look
                else:
                    test_roi2.clock.reset() # keep clock at 0 if roi hasn't started / has finished
                    test_roi2.wasLookedIn = False  # if test_roi2 is looked at next frame, it is a new look
                if test_roi2.status == STARTED:
                    # is it time to stop? (based on global clock, using actual start)
                    if tThisFlipGlobal > test_roi2.tStartRefresh + 2-frameTolerance:
                        # keep track of stop time/frame for later
                        test_roi2.tStop = t  # not accounting for scr refresh
                        test_roi2.frameNStop = frameN  # exact frame index
                        test_roi2.status = FINISHED
                
                # *test_image3* updates
                if test_image3.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                    # keep track of start time/frame for later
                    test_image3.frameNStart = frameN  # exact frame index
                    test_image3.tStart = t  # local t and not account for scr refresh
                    test_image3.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(test_image3, 'tStartRefresh')  # time at next scr refresh
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'test_image3.started')
                    test_image3.setAutoDraw(True)
                if test_image3.status == STARTED:
                    # is it time to stop? (based on global clock, using actual start)
                    if tThisFlipGlobal > test_image3.tStartRefresh + 2-frameTolerance:
                        # keep track of stop time/frame for later
                        test_image3.tStop = t  # not accounting for scr refresh
                        test_image3.frameNStop = frameN  # exact frame index
                        # add timestamp to datafile
                        thisExp.timestampOnFlip(win, 'test_image3.stopped')
                        test_image3.setAutoDraw(False)
                if test_roi3.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                    # keep track of start time/frame for later
                    test_roi3.frameNStart = frameN  # exact frame index
                    test_roi3.tStart = t  # local t and not account for scr refresh
                    test_roi3.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(test_roi3, 'tStartRefresh')  # time at next scr refresh
                    test_roi3.status = STARTED
                if test_roi3.status == STARTED:
                    # check whether test_roi3 has been looked in
                    if test_roi3.isLookedIn:
                        if not test_roi3.wasLookedIn:
                            test_roi3.timesOn.append(routineTimer.getTime()) # store time of first look
                            test_roi3.timesOff.append(routineTimer.getTime()) # store time looked until
                        else:
                            test_roi3.timesOff[-1] = routineTimer.getTime() # update time looked until
                        test_roi3.wasLookedIn = True  # if test_roi3 is still looked at next frame, it is not a new look
                    else:
                        if test_roi3.wasLookedIn:
                            test_roi3.timesOff[-1] = routineTimer.getTime() # update time looked until
                        test_roi3.wasLookedIn = False  # if test_roi3 is looked at next frame, it is a new look
                else:
                    test_roi3.clock.reset() # keep clock at 0 if roi hasn't started / has finished
                    test_roi3.wasLookedIn = False  # if test_roi3 is looked at next frame, it is a new look
                if test_roi3.status == STARTED:
                    # is it time to stop? (based on global clock, using actual start)
                    if tThisFlipGlobal > test_roi3.tStartRefresh + 2-frameTolerance:
                        # keep track of stop time/frame for later
                        test_roi3.tStop = t  # not accounting for scr refresh
                        test_roi3.frameNStop = frameN  # exact frame index
                        test_roi3.status = FINISHED
                
                # *test_image4* updates
                if test_image4.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                    # keep track of start time/frame for later
                    test_image4.frameNStart = frameN  # exact frame index
                    test_image4.tStart = t  # local t and not account for scr refresh
                    test_image4.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(test_image4, 'tStartRefresh')  # time at next scr refresh
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'test_image4.started')
                    test_image4.setAutoDraw(True)
                if test_image4.status == STARTED:
                    # is it time to stop? (based on global clock, using actual start)
                    if tThisFlipGlobal > test_image4.tStartRefresh + 2-frameTolerance:
                        # keep track of stop time/frame for later
                        test_image4.tStop = t  # not accounting for scr refresh
                        test_image4.frameNStop = frameN  # exact frame index
                        # add timestamp to datafile
                        thisExp.timestampOnFlip(win, 'test_image4.stopped')
                        test_image4.setAutoDraw(False)
                if test_roi4.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                    # keep track of start time/frame for later
                    test_roi4.frameNStart = frameN  # exact frame index
                    test_roi4.tStart = t  # local t and not account for scr refresh
                    test_roi4.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(test_roi4, 'tStartRefresh')  # time at next scr refresh
                    test_roi4.status = STARTED
                if test_roi4.status == STARTED:
                    # check whether test_roi4 has been looked in
                    if test_roi4.isLookedIn:
                        if not test_roi4.wasLookedIn:
                            test_roi4.timesOn.append(routineTimer.getTime()) # store time of first look
                            test_roi4.timesOff.append(routineTimer.getTime()) # store time looked until
                        else:
                            test_roi4.timesOff[-1] = routineTimer.getTime() # update time looked until
                        test_roi4.wasLookedIn = True  # if test_roi4 is still looked at next frame, it is not a new look
                    else:
                        if test_roi4.wasLookedIn:
                            test_roi4.timesOff[-1] = routineTimer.getTime() # update time looked until
                        test_roi4.wasLookedIn = False  # if test_roi4 is looked at next frame, it is a new look
                else:
                    test_roi4.clock.reset() # keep clock at 0 if roi hasn't started / has finished
                    test_roi4.wasLookedIn = False  # if test_roi4 is looked at next frame, it is a new look
                if test_roi4.status == STARTED:
                    # is it time to stop? (based on global clock, using actual start)
                    if tThisFlipGlobal > test_roi4.tStartRefresh + 2-frameTolerance:
                        # keep track of stop time/frame for later
                        test_roi4.tStop = t  # not accounting for scr refresh
                        test_roi4.frameNStop = frameN  # exact frame index
                        test_roi4.status = FINISHED
                
                # *gazeCursorTest* updates
                if gazeCursorTest.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                    # keep track of start time/frame for later
                    gazeCursorTest.frameNStart = frameN  # exact frame index
                    gazeCursorTest.tStart = t  # local t and not account for scr refresh
                    gazeCursorTest.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(gazeCursorTest, 'tStartRefresh')  # time at next scr refresh
                    gazeCursorTest.setAutoDraw(True)
                if gazeCursorTest.status == STARTED:
                    # is it time to stop? (based on global clock, using actual start)
                    if tThisFlipGlobal > gazeCursorTest.tStartRefresh + 2-frameTolerance:
                        # keep track of stop time/frame for later
                        gazeCursorTest.tStop = t  # not accounting for scr refresh
                        gazeCursorTest.frameNStop = frameN  # exact frame index
                        gazeCursorTest.setAutoDraw(False)
                if gazeCursorTest.status == STARTED:  # only update if drawing
                    gazeCursorTest.setFillColor(cursorcolor, log=False)
                    gazeCursorTest.setOpacity(0.0, log=False)
                    gazeCursorTest.setPos([eyetracker.getPos()], log=False)
                
                # check for quit (typically the Esc key)
                if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
                    core.quit()
                
                # check if all components have finished
                if not continueRoutine:  # a component has requested a forced-end of Routine
                    routineForceEnded = True
                    break
                continueRoutine = False  # will revert to True if at least one component still running
                for thisComponent in test_trialComponents:
                    if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                        continueRoutine = True
                        break  # at least one component has not yet finished
                
                # refresh the screen
                if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                    win.flip()
            
            # --- Ending Routine "test_trial" ---
            for thisComponent in test_trialComponents:
                if hasattr(thisComponent, "setAutoDraw"):
                    thisComponent.setAutoDraw(False)
            # Run 'End Routine' code from codeTestTrial
            if looked_at_roi1:
                feedback_position = position1
                if "neg" in stim1:
                    outcome = "loss"
                elif "pos" in stim1:
                    outcome = "reward"
                else:
                    outcome = "none"
            if looked_at_roi2:
                feedback_position = position2
                if "neg" in stim2:
                    outcome = "loss"
                elif "pos" in stim2:
                    outcome = "reward"
                else:
                    outcome = "none"
            if looked_at_roi3:
                feedback_position = position3
                if "neg" in stim3:
                    outcome = "loss"
                elif "pos" in stim3:
                    outcome = "reward"
                else:
                    outcome = "none"
            if looked_at_roi4:
                feedback_position = position4
                if "neg" in stim4:
                    outcome = "loss"
                elif "pos" in stim4:
                    outcome = "reward"
                else:
                    outcome = "none"
                    
            if outcome == "loss":
                feedback_color = "red"
                feedback_opacity = 1
                cursorcolor="red"
                feedback_audio = "audio/error.wav"
                log_msg = "lookedat_cs_plus"
                points = -3
                score += points
                feedback_points = f"{points}"
                feedback_score = f"Score: {score}"
            elif outcome == "reward":
                feedback_color = "green"
                feedback_opacity = 1
                cursorcolor="green"
                feedback_audio = "audio/win.wav"
                log_msg = "lookedat_cs_minus"
                points = 3
                score += points
                feedback_points = f"+{points}"
                feedback_score = f"Score: {score}"
            else:
                feedback_color = (0,0,0)
                feedback_opacity = 0
                feedback_audio = "audio/silence.wav"
                log_msg = "no_look"
                points = 0
                feedback_points = ""
                feedback_score = ""
            
            rectsize = [item * 1.02 for item in test_image1.size]
            
            thisExp.addData('score', score)
            thisExp.addData('outcome', outcome)
            
            loop_test.addData('test_roi1.numLooks', test_roi1.numLooks)
            if test_roi1.numLooks:
               loop_test.addData('test_roi1.timesOn', test_roi1.timesOn[0])
               loop_test.addData('test_roi1.timesOff', test_roi1.timesOff[0])
            else:
               loop_test.addData('test_roi1.timesOn', "")
               loop_test.addData('test_roi1.timesOff', "")
            loop_test.addData('test_roi2.numLooks', test_roi2.numLooks)
            if test_roi2.numLooks:
               loop_test.addData('test_roi2.timesOn', test_roi2.timesOn[0])
               loop_test.addData('test_roi2.timesOff', test_roi2.timesOff[0])
            else:
               loop_test.addData('test_roi2.timesOn', "")
               loop_test.addData('test_roi2.timesOff', "")
            loop_test.addData('test_roi3.numLooks', test_roi3.numLooks)
            if test_roi3.numLooks:
               loop_test.addData('test_roi3.timesOn', test_roi3.timesOn[0])
               loop_test.addData('test_roi3.timesOff', test_roi3.timesOff[0])
            else:
               loop_test.addData('test_roi3.timesOn', "")
               loop_test.addData('test_roi3.timesOff', "")
            loop_test.addData('test_roi4.numLooks', test_roi4.numLooks)
            if test_roi4.numLooks:
               loop_test.addData('test_roi4.timesOn', test_roi4.timesOn[0])
               loop_test.addData('test_roi4.timesOff', test_roi4.timesOff[0])
            else:
               loop_test.addData('test_roi4.timesOn', "")
               loop_test.addData('test_roi4.timesOff', "")
            # using non-slip timing so subtract the expected duration of this Routine (unless ended on request)
            if routineForceEnded:
                routineTimer.reset()
            else:
                routineTimer.addTime(-2.000000)
            
            # --- Prepare to start Routine "test_feedback" ---
            continueRoutine = True
            routineForceEnded = False
            # update component parameters for each repeat
            # Run 'Begin Routine' code from codeFeedbackTest
            logging.log(level=logging.INFO, msg=f'FeedbackTestOnset_{log_msg}')
            ioServer.sendMessageEvent(text='FeedbackTestOnset')
            eyetracker.sendMessage('FeedbackTestOnset')
            eyetracker.sendMessage(outcome)
            
            polygonFeedbackTest.setOpacity(feedback_opacity)
            polygonFeedbackTest.setPos(feedback_position)
            polygonFeedbackTest.setLineColor(feedback_color)
            test_image_feedback1.setPos(position1)
            test_image_feedback1.setImage(eval(stim1))
            test_image_feedback2.setPos(position2)
            test_image_feedback2.setImage(eval(stim2))
            test_iamge_feedback3.setPos(position3)
            test_iamge_feedback3.setImage(eval(stim3))
            test_iamge_feedback4.setPos(position4)
            test_iamge_feedback4.setImage(eval(stim4))
            soundFeedbackTest.setSound(feedback_audio, secs=0.5, hamming=True)
            soundFeedbackTest.setVolume(1.0, log=False)
            # keep track of which components have finished
            test_feedbackComponents = [textScoreTest, polygonFeedbackTest, test_image_feedback1, test_image_feedback2, test_iamge_feedback3, test_iamge_feedback4, gazeCursorTestFeedback, soundFeedbackTest]
            for thisComponent in test_feedbackComponents:
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
            
            # --- Run Routine "test_feedback" ---
            while continueRoutine and routineTimer.getTime() < 0.5:
                # get current time
                t = routineTimer.getTime()
                tThisFlip = win.getFutureFlipTime(clock=routineTimer)
                tThisFlipGlobal = win.getFutureFlipTime(clock=None)
                frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
                # update/draw components on each frame
                # Run 'Each Frame' code from codeFeedbackTest
                if loop_timer.getTime() >= 2.0:
                    continueRoutine = False
                    currentLoop.finished = True
                
                # *textScoreTest* updates
                if textScoreTest.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                    # keep track of start time/frame for later
                    textScoreTest.frameNStart = frameN  # exact frame index
                    textScoreTest.tStart = t  # local t and not account for scr refresh
                    textScoreTest.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(textScoreTest, 'tStartRefresh')  # time at next scr refresh
                    textScoreTest.setAutoDraw(True)
                if textScoreTest.status == STARTED:
                    # is it time to stop? (based on global clock, using actual start)
                    if tThisFlipGlobal > textScoreTest.tStartRefresh + 0.5-frameTolerance:
                        # keep track of stop time/frame for later
                        textScoreTest.tStop = t  # not accounting for scr refresh
                        textScoreTest.frameNStop = frameN  # exact frame index
                        textScoreTest.setAutoDraw(False)
                if textScoreTest.status == STARTED:  # only update if drawing
                    textScoreTest.setText(feedback_score, log=False)
                
                # *polygonFeedbackTest* updates
                if polygonFeedbackTest.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                    # keep track of start time/frame for later
                    polygonFeedbackTest.frameNStart = frameN  # exact frame index
                    polygonFeedbackTest.tStart = t  # local t and not account for scr refresh
                    polygonFeedbackTest.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(polygonFeedbackTest, 'tStartRefresh')  # time at next scr refresh
                    polygonFeedbackTest.setAutoDraw(True)
                if polygonFeedbackTest.status == STARTED:
                    # is it time to stop? (based on global clock, using actual start)
                    if tThisFlipGlobal > polygonFeedbackTest.tStartRefresh + 0.5-frameTolerance:
                        # keep track of stop time/frame for later
                        polygonFeedbackTest.tStop = t  # not accounting for scr refresh
                        polygonFeedbackTest.frameNStop = frameN  # exact frame index
                        polygonFeedbackTest.setAutoDraw(False)
                if polygonFeedbackTest.status == STARTED:  # only update if drawing
                    polygonFeedbackTest.setSize(rectsize, log=False)
                
                # *test_image_feedback1* updates
                if test_image_feedback1.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                    # keep track of start time/frame for later
                    test_image_feedback1.frameNStart = frameN  # exact frame index
                    test_image_feedback1.tStart = t  # local t and not account for scr refresh
                    test_image_feedback1.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(test_image_feedback1, 'tStartRefresh')  # time at next scr refresh
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'test_image_feedback1.started')
                    test_image_feedback1.setAutoDraw(True)
                if test_image_feedback1.status == STARTED:
                    # is it time to stop? (based on global clock, using actual start)
                    if tThisFlipGlobal > test_image_feedback1.tStartRefresh + 0.5-frameTolerance:
                        # keep track of stop time/frame for later
                        test_image_feedback1.tStop = t  # not accounting for scr refresh
                        test_image_feedback1.frameNStop = frameN  # exact frame index
                        # add timestamp to datafile
                        thisExp.timestampOnFlip(win, 'test_image_feedback1.stopped')
                        test_image_feedback1.setAutoDraw(False)
                
                # *test_image_feedback2* updates
                if test_image_feedback2.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                    # keep track of start time/frame for later
                    test_image_feedback2.frameNStart = frameN  # exact frame index
                    test_image_feedback2.tStart = t  # local t and not account for scr refresh
                    test_image_feedback2.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(test_image_feedback2, 'tStartRefresh')  # time at next scr refresh
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'test_image_feedback2.started')
                    test_image_feedback2.setAutoDraw(True)
                if test_image_feedback2.status == STARTED:
                    # is it time to stop? (based on global clock, using actual start)
                    if tThisFlipGlobal > test_image_feedback2.tStartRefresh + 0.5-frameTolerance:
                        # keep track of stop time/frame for later
                        test_image_feedback2.tStop = t  # not accounting for scr refresh
                        test_image_feedback2.frameNStop = frameN  # exact frame index
                        # add timestamp to datafile
                        thisExp.timestampOnFlip(win, 'test_image_feedback2.stopped')
                        test_image_feedback2.setAutoDraw(False)
                
                # *test_iamge_feedback3* updates
                if test_iamge_feedback3.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                    # keep track of start time/frame for later
                    test_iamge_feedback3.frameNStart = frameN  # exact frame index
                    test_iamge_feedback3.tStart = t  # local t and not account for scr refresh
                    test_iamge_feedback3.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(test_iamge_feedback3, 'tStartRefresh')  # time at next scr refresh
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'test_iamge_feedback3.started')
                    test_iamge_feedback3.setAutoDraw(True)
                if test_iamge_feedback3.status == STARTED:
                    # is it time to stop? (based on global clock, using actual start)
                    if tThisFlipGlobal > test_iamge_feedback3.tStartRefresh + 0.5-frameTolerance:
                        # keep track of stop time/frame for later
                        test_iamge_feedback3.tStop = t  # not accounting for scr refresh
                        test_iamge_feedback3.frameNStop = frameN  # exact frame index
                        # add timestamp to datafile
                        thisExp.timestampOnFlip(win, 'test_iamge_feedback3.stopped')
                        test_iamge_feedback3.setAutoDraw(False)
                
                # *test_iamge_feedback4* updates
                if test_iamge_feedback4.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                    # keep track of start time/frame for later
                    test_iamge_feedback4.frameNStart = frameN  # exact frame index
                    test_iamge_feedback4.tStart = t  # local t and not account for scr refresh
                    test_iamge_feedback4.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(test_iamge_feedback4, 'tStartRefresh')  # time at next scr refresh
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'test_iamge_feedback4.started')
                    test_iamge_feedback4.setAutoDraw(True)
                if test_iamge_feedback4.status == STARTED:
                    # is it time to stop? (based on global clock, using actual start)
                    if tThisFlipGlobal > test_iamge_feedback4.tStartRefresh + 0.5-frameTolerance:
                        # keep track of stop time/frame for later
                        test_iamge_feedback4.tStop = t  # not accounting for scr refresh
                        test_iamge_feedback4.frameNStop = frameN  # exact frame index
                        # add timestamp to datafile
                        thisExp.timestampOnFlip(win, 'test_iamge_feedback4.stopped')
                        test_iamge_feedback4.setAutoDraw(False)
                
                # *gazeCursorTestFeedback* updates
                if gazeCursorTestFeedback.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                    # keep track of start time/frame for later
                    gazeCursorTestFeedback.frameNStart = frameN  # exact frame index
                    gazeCursorTestFeedback.tStart = t  # local t and not account for scr refresh
                    gazeCursorTestFeedback.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(gazeCursorTestFeedback, 'tStartRefresh')  # time at next scr refresh
                    gazeCursorTestFeedback.setAutoDraw(True)
                if gazeCursorTestFeedback.status == STARTED:
                    # is it time to stop? (based on global clock, using actual start)
                    if tThisFlipGlobal > gazeCursorTestFeedback.tStartRefresh + 0.5-frameTolerance:
                        # keep track of stop time/frame for later
                        gazeCursorTestFeedback.tStop = t  # not accounting for scr refresh
                        gazeCursorTestFeedback.frameNStop = frameN  # exact frame index
                        gazeCursorTestFeedback.setAutoDraw(False)
                if gazeCursorTestFeedback.status == STARTED:  # only update if drawing
                    gazeCursorTestFeedback.setFillColor(cursorcolor, log=False)
                    gazeCursorTestFeedback.setOpacity(0.0, log=False)
                    gazeCursorTestFeedback.setPos([eyetracker.getPos()], log=False)
                # start/stop soundFeedbackTest
                if soundFeedbackTest.status == NOT_STARTED and t >= 0-frameTolerance:
                    # keep track of start time/frame for later
                    soundFeedbackTest.frameNStart = frameN  # exact frame index
                    soundFeedbackTest.tStart = t  # local t and not account for scr refresh
                    soundFeedbackTest.tStartRefresh = tThisFlipGlobal  # on global time
                    soundFeedbackTest.play()  # start the sound (it finishes automatically)
                if soundFeedbackTest.status == STARTED:
                    # is it time to stop? (based on global clock, using actual start)
                    if tThisFlipGlobal > soundFeedbackTest.tStartRefresh + 0.5-frameTolerance:
                        # keep track of stop time/frame for later
                        soundFeedbackTest.tStop = t  # not accounting for scr refresh
                        soundFeedbackTest.frameNStop = frameN  # exact frame index
                        soundFeedbackTest.stop()
                
                # check for quit (typically the Esc key)
                if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
                    core.quit()
                
                # check if all components have finished
                if not continueRoutine:  # a component has requested a forced-end of Routine
                    routineForceEnded = True
                    break
                continueRoutine = False  # will revert to True if at least one component still running
                for thisComponent in test_feedbackComponents:
                    if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                        continueRoutine = True
                        break  # at least one component has not yet finished
                
                # refresh the screen
                if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                    win.flip()
            
            # --- Ending Routine "test_feedback" ---
            for thisComponent in test_feedbackComponents:
                if hasattr(thisComponent, "setAutoDraw"):
                    thisComponent.setAutoDraw(False)
            soundFeedbackTest.stop()  # ensure sound has stopped at end of routine
            # using non-slip timing so subtract the expected duration of this Routine (unless ended on request)
            if routineForceEnded:
                routineTimer.reset()
            else:
                routineTimer.addTime(-0.500000)
        # completed 100.0 repeats of 'loop_test'
        
        
        # --- Prepare to start Routine "crossEnd" ---
        continueRoutine = True
        routineForceEnded = False
        # update component parameters for each repeat
        # Run 'Begin Routine' code from code_end
        eyetracker.setRecordingState(False)
        jittered_duration_cross = 3.0 + random()
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
                fixcrossE.setAutoDraw(True)
            if fixcrossE.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > fixcrossE.tStartRefresh + jittered_duration_cross-frameTolerance:
                    # keep track of stop time/frame for later
                    fixcrossE.tStop = t  # not accounting for scr refresh
                    fixcrossE.frameNStop = frameN  # exact frame index
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
        
    # completed 2.0 repeats of 'testtrials'
    
    
    # set up handler to look after randomisation of conditions etc
    testtrials_novelty = data.TrialHandler(nReps=1.0, method='random', 
        extraInfo=expInfo, originPath=-1,
        trialList=data.importConditions(posFileTestNew, selection=np.concatenate((np.random.permutation(np.arange(0, 24))[0:8], np.random.permutation(np.arange(24, 48))[0:10], np.random.permutation(np.arange(48, 72))[0:10], np.random.permutation(np.arange(72, 96))[0:10], np.random.permutation(np.arange(96, 120))[0:10]))),
        seed=None, name='testtrials_novelty')
    thisExp.addLoop(testtrials_novelty)  # add the loop to the experiment
    thisTesttrials_novelty = testtrials_novelty.trialList[0]  # so we can initialise stimuli with some values
    # abbreviate parameter names if possible (e.g. rgb = thisTesttrials_novelty.rgb)
    if thisTesttrials_novelty != None:
        for paramName in thisTesttrials_novelty:
            exec('{} = thisTesttrials_novelty[paramName]'.format(paramName))
    
    for thisTesttrials_novelty in testtrials_novelty:
        currentLoop = testtrials_novelty
        # abbreviate parameter names if possible (e.g. rgb = thisTesttrials_novelty.rgb)
        if thisTesttrials_novelty != None:
            for paramName in thisTesttrials_novelty:
                exec('{} = thisTesttrials_novelty[paramName]'.format(paramName))
        
        # --- Prepare to start Routine "cross" ---
        continueRoutine = True
        routineForceEnded = False
        # update component parameters for each repeat
        # Run 'Begin Routine' code from codeFixate
        logging.log(level=logging.INFO, msg=f'FixationCross')
        eyetracker.setRecordingState(True)
         
        if currentLoop.name == "learning_trials":
            eyetracker.sendMessage("Trial " + str(learning_trials.thisN+1) + ", Condition " + stim)
            print("VP %s TRIALID %d CONDITION %s"%(expInfo['participant'], learning_trials.thisN+1, stim))
        
        if currentLoop.name == "trials":
            eyetracker.sendMessage("Trial " + str(trials.thisN+learning_trials.thisN+1) + ", Condition " + stim)
            print("VP %s TRIALID %d CONDITION %s"%(expInfo['participant'], trials.thisN+learning_trials.thisN+1, stim))
        
        if currentLoop.name == "testtrials":
            eyetracker.sendMessage("Test-Trial " + str(testtrials.thisN+1))
            print("VP %s TESTTRIALID %d"%(expInfo['participant'], testtrials.thisN+1))
            
        if currentLoop.name == "testtrials_novelty":
            eyetracker.sendMessage("Test-Trial " + str(testtrials_novelty.thisN+testtrials.thisN+1))
            print("VP %s TESTTRIALID %d"%(expInfo['participant'], testtrials_novelty.thisN+testtrials.thisN+1))
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
                fixcross.setAutoDraw(True)
            if fixcross.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > fixcross.tStartRefresh + 1-frameTolerance:
                    # keep track of stop time/frame for later
                    fixcross.tStop = t  # not accounting for scr refresh
                    fixcross.frameNStop = frameN  # exact frame index
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
        
        # --- Prepare to start Routine "blank_test" ---
        continueRoutine = True
        routineForceEnded = False
        # update component parameters for each repeat
        # Run 'Begin Routine' code from codeBlankTest
        jittered_duration_blank = 1.0 + random()
        looked_at_roi1_trial = False
        looked_at_roi2_trial = False
        looked_at_roi3_trial = False
        looked_at_roi4_trial = False
        blankScreenTest.setText('')
        # keep track of which components have finished
        blank_testComponents = [blankScreenTest]
        for thisComponent in blank_testComponents:
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
        
        # --- Run Routine "blank_test" ---
        while continueRoutine:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # *blankScreenTest* updates
            if blankScreenTest.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                blankScreenTest.frameNStart = frameN  # exact frame index
                blankScreenTest.tStart = t  # local t and not account for scr refresh
                blankScreenTest.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(blankScreenTest, 'tStartRefresh')  # time at next scr refresh
                blankScreenTest.setAutoDraw(True)
            if blankScreenTest.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > blankScreenTest.tStartRefresh + jittered_duration_blank-frameTolerance:
                    # keep track of stop time/frame for later
                    blankScreenTest.tStop = t  # not accounting for scr refresh
                    blankScreenTest.frameNStop = frameN  # exact frame index
                    blankScreenTest.setAutoDraw(False)
            
            # check for quit (typically the Esc key)
            if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
                core.quit()
            
            # check if all components have finished
            if not continueRoutine:  # a component has requested a forced-end of Routine
                routineForceEnded = True
                break
            continueRoutine = False  # will revert to True if at least one component still running
            for thisComponent in blank_testComponents:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # --- Ending Routine "blank_test" ---
        for thisComponent in blank_testComponents:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        # Run 'End Routine' code from codeBlankTest
        # print(f"Start: {round(blankScreen.tStart, 3)}, End: {round(blankScreen.tStop, 3)}, Duration: {round(blankScreen.tStop - blankScreen.tStart, 3)}")
        loop_timer = core.Clock()
        # the Routine "blank_test" was not non-slip safe, so reset the non-slip timer
        routineTimer.reset()
        
        # set up handler to look after randomisation of conditions etc
        loop_test_novelty = data.TrialHandler(nReps=100.0, method='sequential', 
            extraInfo=expInfo, originPath=-1,
            trialList=[None],
            seed=None, name='loop_test_novelty')
        thisExp.addLoop(loop_test_novelty)  # add the loop to the experiment
        thisLoop_test_novelty = loop_test_novelty.trialList[0]  # so we can initialise stimuli with some values
        # abbreviate parameter names if possible (e.g. rgb = thisLoop_test_novelty.rgb)
        if thisLoop_test_novelty != None:
            for paramName in thisLoop_test_novelty:
                exec('{} = thisLoop_test_novelty[paramName]'.format(paramName))
        
        for thisLoop_test_novelty in loop_test_novelty:
            currentLoop = loop_test_novelty
            # abbreviate parameter names if possible (e.g. rgb = thisLoop_test_novelty.rgb)
            if thisLoop_test_novelty != None:
                for paramName in thisLoop_test_novelty:
                    exec('{} = thisLoop_test_novelty[paramName]'.format(paramName))
            
            # --- Prepare to start Routine "test_trial" ---
            continueRoutine = True
            routineForceEnded = False
            # update component parameters for each repeat
            # Run 'Begin Routine' code from codeTestTrial
            looked_at_roi1 = False
            looked_at_roi2 = False
            looked_at_roi3 = False
            looked_at_roi4 = False
            cursorcolor="white"
            
            logging.log(level=logging.INFO, msg=f'ImageTestOnset')
            ioServer.sendMessageEvent(text='ImageTestOnset')
            eyetracker.sendMessage('ImageTestOnset')
            test_image1.setPos(position1)
            test_image1.setImage(eval(stim1))
            test_roi1.setPos(position1)
            test_roi1.setSize([test_image1.size])
            # clear any previous roi data
            test_roi1.reset()
            test_image2.setPos(position2)
            test_image2.setImage(eval(stim2))
            test_roi2.setPos(position2)
            test_roi2.setSize([test_image2.size])
            # clear any previous roi data
            test_roi2.reset()
            test_image3.setPos(position3)
            test_image3.setImage(eval(stim3))
            test_roi3.setPos(position3)
            test_roi3.setSize([test_image3.size])
            # clear any previous roi data
            test_roi3.reset()
            test_image4.setPos(position4)
            test_image4.setImage(eval(stim4))
            test_roi4.setPos(position4)
            test_roi4.setSize([test_image4.size])
            # clear any previous roi data
            test_roi4.reset()
            # keep track of which components have finished
            test_trialComponents = [test_image1, test_roi1, test_image2, test_roi2, test_image3, test_roi3, test_image4, test_roi4, gazeCursorTest]
            for thisComponent in test_trialComponents:
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
            
            # --- Run Routine "test_trial" ---
            while continueRoutine and routineTimer.getTime() < 2.0:
                # get current time
                t = routineTimer.getTime()
                tThisFlip = win.getFutureFlipTime(clock=routineTimer)
                tThisFlipGlobal = win.getFutureFlipTime(clock=None)
                frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
                # update/draw components on each frame
                # Run 'Each Frame' code from codeTestTrial
                if loop_timer.getTime() >= 2.0:
                    continueRoutine = False
                    currentLoop.finished = True
                
                test_roi1.size = test_image1.size
                test_roi1.pos = position1
                
                if test_roi1.isLookedIn and not looked_at_roi1_trial:
                    looked_at_roi1 = True
                    looked_at_roi1_trial = True
                    continueRoutine = False
                
                test_roi2.size = test_image2.size
                test_roi2.pos = position2
                
                if test_roi2.isLookedIn and not looked_at_roi2_trial:
                    looked_at_roi2 = True
                    looked_at_roi2_trial = True
                    continueRoutine = False
                    
                test_roi3.size = test_image3.size
                test_roi3.pos = position3
                
                if test_roi3.isLookedIn and not looked_at_roi3_trial:
                    looked_at_roi3 = True
                    looked_at_roi3_trial = True
                    continueRoutine = False
                   
                test_roi4.size = test_image4.size
                test_roi4.pos = position4
                
                if test_roi4.isLookedIn and not looked_at_roi4_trial:
                    looked_at_roi4 = True
                    looked_at_roi4_trial = True
                    continueRoutine = False
                    
                if event.getKeys(keyList=["p"]) and not paused:
                    print("Experiment Paused - Press 'p' to continue.")
                    paused = True
                    event.waitKeys(keyList=["p"], clearEvents=True)
                    
                    print("Experiment Continued")
                    paused = False
                    event.clearEvents(eventType="keyboard")
                    
                
                # *test_image1* updates
                if test_image1.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                    # keep track of start time/frame for later
                    test_image1.frameNStart = frameN  # exact frame index
                    test_image1.tStart = t  # local t and not account for scr refresh
                    test_image1.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(test_image1, 'tStartRefresh')  # time at next scr refresh
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'test_image1.started')
                    test_image1.setAutoDraw(True)
                if test_image1.status == STARTED:
                    # is it time to stop? (based on global clock, using actual start)
                    if tThisFlipGlobal > test_image1.tStartRefresh + 2-frameTolerance:
                        # keep track of stop time/frame for later
                        test_image1.tStop = t  # not accounting for scr refresh
                        test_image1.frameNStop = frameN  # exact frame index
                        # add timestamp to datafile
                        thisExp.timestampOnFlip(win, 'test_image1.stopped')
                        test_image1.setAutoDraw(False)
                if test_roi1.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                    # keep track of start time/frame for later
                    test_roi1.frameNStart = frameN  # exact frame index
                    test_roi1.tStart = t  # local t and not account for scr refresh
                    test_roi1.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(test_roi1, 'tStartRefresh')  # time at next scr refresh
                    test_roi1.status = STARTED
                if test_roi1.status == STARTED:
                    # check whether test_roi1 has been looked in
                    if test_roi1.isLookedIn:
                        if not test_roi1.wasLookedIn:
                            test_roi1.timesOn.append(routineTimer.getTime()) # store time of first look
                            test_roi1.timesOff.append(routineTimer.getTime()) # store time looked until
                        else:
                            test_roi1.timesOff[-1] = routineTimer.getTime() # update time looked until
                        test_roi1.wasLookedIn = True  # if test_roi1 is still looked at next frame, it is not a new look
                    else:
                        if test_roi1.wasLookedIn:
                            test_roi1.timesOff[-1] = routineTimer.getTime() # update time looked until
                        test_roi1.wasLookedIn = False  # if test_roi1 is looked at next frame, it is a new look
                else:
                    test_roi1.clock.reset() # keep clock at 0 if roi hasn't started / has finished
                    test_roi1.wasLookedIn = False  # if test_roi1 is looked at next frame, it is a new look
                if test_roi1.status == STARTED:
                    # is it time to stop? (based on global clock, using actual start)
                    if tThisFlipGlobal > test_roi1.tStartRefresh + 2-frameTolerance:
                        # keep track of stop time/frame for later
                        test_roi1.tStop = t  # not accounting for scr refresh
                        test_roi1.frameNStop = frameN  # exact frame index
                        test_roi1.status = FINISHED
                
                # *test_image2* updates
                if test_image2.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                    # keep track of start time/frame for later
                    test_image2.frameNStart = frameN  # exact frame index
                    test_image2.tStart = t  # local t and not account for scr refresh
                    test_image2.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(test_image2, 'tStartRefresh')  # time at next scr refresh
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'test_image2.started')
                    test_image2.setAutoDraw(True)
                if test_image2.status == STARTED:
                    # is it time to stop? (based on global clock, using actual start)
                    if tThisFlipGlobal > test_image2.tStartRefresh + 2-frameTolerance:
                        # keep track of stop time/frame for later
                        test_image2.tStop = t  # not accounting for scr refresh
                        test_image2.frameNStop = frameN  # exact frame index
                        # add timestamp to datafile
                        thisExp.timestampOnFlip(win, 'test_image2.stopped')
                        test_image2.setAutoDraw(False)
                if test_roi2.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                    # keep track of start time/frame for later
                    test_roi2.frameNStart = frameN  # exact frame index
                    test_roi2.tStart = t  # local t and not account for scr refresh
                    test_roi2.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(test_roi2, 'tStartRefresh')  # time at next scr refresh
                    test_roi2.status = STARTED
                if test_roi2.status == STARTED:
                    # check whether test_roi2 has been looked in
                    if test_roi2.isLookedIn:
                        if not test_roi2.wasLookedIn:
                            test_roi2.timesOn.append(routineTimer.getTime()) # store time of first look
                            test_roi2.timesOff.append(routineTimer.getTime()) # store time looked until
                        else:
                            test_roi2.timesOff[-1] = routineTimer.getTime() # update time looked until
                        test_roi2.wasLookedIn = True  # if test_roi2 is still looked at next frame, it is not a new look
                    else:
                        if test_roi2.wasLookedIn:
                            test_roi2.timesOff[-1] = routineTimer.getTime() # update time looked until
                        test_roi2.wasLookedIn = False  # if test_roi2 is looked at next frame, it is a new look
                else:
                    test_roi2.clock.reset() # keep clock at 0 if roi hasn't started / has finished
                    test_roi2.wasLookedIn = False  # if test_roi2 is looked at next frame, it is a new look
                if test_roi2.status == STARTED:
                    # is it time to stop? (based on global clock, using actual start)
                    if tThisFlipGlobal > test_roi2.tStartRefresh + 2-frameTolerance:
                        # keep track of stop time/frame for later
                        test_roi2.tStop = t  # not accounting for scr refresh
                        test_roi2.frameNStop = frameN  # exact frame index
                        test_roi2.status = FINISHED
                
                # *test_image3* updates
                if test_image3.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                    # keep track of start time/frame for later
                    test_image3.frameNStart = frameN  # exact frame index
                    test_image3.tStart = t  # local t and not account for scr refresh
                    test_image3.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(test_image3, 'tStartRefresh')  # time at next scr refresh
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'test_image3.started')
                    test_image3.setAutoDraw(True)
                if test_image3.status == STARTED:
                    # is it time to stop? (based on global clock, using actual start)
                    if tThisFlipGlobal > test_image3.tStartRefresh + 2-frameTolerance:
                        # keep track of stop time/frame for later
                        test_image3.tStop = t  # not accounting for scr refresh
                        test_image3.frameNStop = frameN  # exact frame index
                        # add timestamp to datafile
                        thisExp.timestampOnFlip(win, 'test_image3.stopped')
                        test_image3.setAutoDraw(False)
                if test_roi3.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                    # keep track of start time/frame for later
                    test_roi3.frameNStart = frameN  # exact frame index
                    test_roi3.tStart = t  # local t and not account for scr refresh
                    test_roi3.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(test_roi3, 'tStartRefresh')  # time at next scr refresh
                    test_roi3.status = STARTED
                if test_roi3.status == STARTED:
                    # check whether test_roi3 has been looked in
                    if test_roi3.isLookedIn:
                        if not test_roi3.wasLookedIn:
                            test_roi3.timesOn.append(routineTimer.getTime()) # store time of first look
                            test_roi3.timesOff.append(routineTimer.getTime()) # store time looked until
                        else:
                            test_roi3.timesOff[-1] = routineTimer.getTime() # update time looked until
                        test_roi3.wasLookedIn = True  # if test_roi3 is still looked at next frame, it is not a new look
                    else:
                        if test_roi3.wasLookedIn:
                            test_roi3.timesOff[-1] = routineTimer.getTime() # update time looked until
                        test_roi3.wasLookedIn = False  # if test_roi3 is looked at next frame, it is a new look
                else:
                    test_roi3.clock.reset() # keep clock at 0 if roi hasn't started / has finished
                    test_roi3.wasLookedIn = False  # if test_roi3 is looked at next frame, it is a new look
                if test_roi3.status == STARTED:
                    # is it time to stop? (based on global clock, using actual start)
                    if tThisFlipGlobal > test_roi3.tStartRefresh + 2-frameTolerance:
                        # keep track of stop time/frame for later
                        test_roi3.tStop = t  # not accounting for scr refresh
                        test_roi3.frameNStop = frameN  # exact frame index
                        test_roi3.status = FINISHED
                
                # *test_image4* updates
                if test_image4.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                    # keep track of start time/frame for later
                    test_image4.frameNStart = frameN  # exact frame index
                    test_image4.tStart = t  # local t and not account for scr refresh
                    test_image4.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(test_image4, 'tStartRefresh')  # time at next scr refresh
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'test_image4.started')
                    test_image4.setAutoDraw(True)
                if test_image4.status == STARTED:
                    # is it time to stop? (based on global clock, using actual start)
                    if tThisFlipGlobal > test_image4.tStartRefresh + 2-frameTolerance:
                        # keep track of stop time/frame for later
                        test_image4.tStop = t  # not accounting for scr refresh
                        test_image4.frameNStop = frameN  # exact frame index
                        # add timestamp to datafile
                        thisExp.timestampOnFlip(win, 'test_image4.stopped')
                        test_image4.setAutoDraw(False)
                if test_roi4.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                    # keep track of start time/frame for later
                    test_roi4.frameNStart = frameN  # exact frame index
                    test_roi4.tStart = t  # local t and not account for scr refresh
                    test_roi4.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(test_roi4, 'tStartRefresh')  # time at next scr refresh
                    test_roi4.status = STARTED
                if test_roi4.status == STARTED:
                    # check whether test_roi4 has been looked in
                    if test_roi4.isLookedIn:
                        if not test_roi4.wasLookedIn:
                            test_roi4.timesOn.append(routineTimer.getTime()) # store time of first look
                            test_roi4.timesOff.append(routineTimer.getTime()) # store time looked until
                        else:
                            test_roi4.timesOff[-1] = routineTimer.getTime() # update time looked until
                        test_roi4.wasLookedIn = True  # if test_roi4 is still looked at next frame, it is not a new look
                    else:
                        if test_roi4.wasLookedIn:
                            test_roi4.timesOff[-1] = routineTimer.getTime() # update time looked until
                        test_roi4.wasLookedIn = False  # if test_roi4 is looked at next frame, it is a new look
                else:
                    test_roi4.clock.reset() # keep clock at 0 if roi hasn't started / has finished
                    test_roi4.wasLookedIn = False  # if test_roi4 is looked at next frame, it is a new look
                if test_roi4.status == STARTED:
                    # is it time to stop? (based on global clock, using actual start)
                    if tThisFlipGlobal > test_roi4.tStartRefresh + 2-frameTolerance:
                        # keep track of stop time/frame for later
                        test_roi4.tStop = t  # not accounting for scr refresh
                        test_roi4.frameNStop = frameN  # exact frame index
                        test_roi4.status = FINISHED
                
                # *gazeCursorTest* updates
                if gazeCursorTest.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                    # keep track of start time/frame for later
                    gazeCursorTest.frameNStart = frameN  # exact frame index
                    gazeCursorTest.tStart = t  # local t and not account for scr refresh
                    gazeCursorTest.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(gazeCursorTest, 'tStartRefresh')  # time at next scr refresh
                    gazeCursorTest.setAutoDraw(True)
                if gazeCursorTest.status == STARTED:
                    # is it time to stop? (based on global clock, using actual start)
                    if tThisFlipGlobal > gazeCursorTest.tStartRefresh + 2-frameTolerance:
                        # keep track of stop time/frame for later
                        gazeCursorTest.tStop = t  # not accounting for scr refresh
                        gazeCursorTest.frameNStop = frameN  # exact frame index
                        gazeCursorTest.setAutoDraw(False)
                if gazeCursorTest.status == STARTED:  # only update if drawing
                    gazeCursorTest.setFillColor(cursorcolor, log=False)
                    gazeCursorTest.setOpacity(0.0, log=False)
                    gazeCursorTest.setPos([eyetracker.getPos()], log=False)
                
                # check for quit (typically the Esc key)
                if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
                    core.quit()
                
                # check if all components have finished
                if not continueRoutine:  # a component has requested a forced-end of Routine
                    routineForceEnded = True
                    break
                continueRoutine = False  # will revert to True if at least one component still running
                for thisComponent in test_trialComponents:
                    if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                        continueRoutine = True
                        break  # at least one component has not yet finished
                
                # refresh the screen
                if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                    win.flip()
            
            # --- Ending Routine "test_trial" ---
            for thisComponent in test_trialComponents:
                if hasattr(thisComponent, "setAutoDraw"):
                    thisComponent.setAutoDraw(False)
            # Run 'End Routine' code from codeTestTrial
            if looked_at_roi1:
                feedback_position = position1
                if "neg" in stim1:
                    outcome = "loss"
                elif "pos" in stim1:
                    outcome = "reward"
                else:
                    outcome = "none"
            if looked_at_roi2:
                feedback_position = position2
                if "neg" in stim2:
                    outcome = "loss"
                elif "pos" in stim2:
                    outcome = "reward"
                else:
                    outcome = "none"
            if looked_at_roi3:
                feedback_position = position3
                if "neg" in stim3:
                    outcome = "loss"
                elif "pos" in stim3:
                    outcome = "reward"
                else:
                    outcome = "none"
            if looked_at_roi4:
                feedback_position = position4
                if "neg" in stim4:
                    outcome = "loss"
                elif "pos" in stim4:
                    outcome = "reward"
                else:
                    outcome = "none"
                    
            if outcome == "loss":
                feedback_color = "red"
                feedback_opacity = 1
                cursorcolor="red"
                feedback_audio = "audio/error.wav"
                log_msg = "lookedat_cs_plus"
                points = -3
                score += points
                feedback_points = f"{points}"
                feedback_score = f"Score: {score}"
            elif outcome == "reward":
                feedback_color = "green"
                feedback_opacity = 1
                cursorcolor="green"
                feedback_audio = "audio/win.wav"
                log_msg = "lookedat_cs_minus"
                points = 3
                score += points
                feedback_points = f"+{points}"
                feedback_score = f"Score: {score}"
            else:
                feedback_color = (0,0,0)
                feedback_opacity = 0
                feedback_audio = "audio/silence.wav"
                log_msg = "no_look"
                points = 0
                feedback_points = ""
                feedback_score = ""
            
            rectsize = [item * 1.02 for item in test_image1.size]
            
            thisExp.addData('score', score)
            thisExp.addData('outcome', outcome)
            
            loop_test_novelty.addData('test_roi1.numLooks', test_roi1.numLooks)
            if test_roi1.numLooks:
               loop_test_novelty.addData('test_roi1.timesOn', test_roi1.timesOn[0])
               loop_test_novelty.addData('test_roi1.timesOff', test_roi1.timesOff[0])
            else:
               loop_test_novelty.addData('test_roi1.timesOn', "")
               loop_test_novelty.addData('test_roi1.timesOff', "")
            loop_test_novelty.addData('test_roi2.numLooks', test_roi2.numLooks)
            if test_roi2.numLooks:
               loop_test_novelty.addData('test_roi2.timesOn', test_roi2.timesOn[0])
               loop_test_novelty.addData('test_roi2.timesOff', test_roi2.timesOff[0])
            else:
               loop_test_novelty.addData('test_roi2.timesOn', "")
               loop_test_novelty.addData('test_roi2.timesOff', "")
            loop_test_novelty.addData('test_roi3.numLooks', test_roi3.numLooks)
            if test_roi3.numLooks:
               loop_test_novelty.addData('test_roi3.timesOn', test_roi3.timesOn[0])
               loop_test_novelty.addData('test_roi3.timesOff', test_roi3.timesOff[0])
            else:
               loop_test_novelty.addData('test_roi3.timesOn', "")
               loop_test_novelty.addData('test_roi3.timesOff', "")
            loop_test_novelty.addData('test_roi4.numLooks', test_roi4.numLooks)
            if test_roi4.numLooks:
               loop_test_novelty.addData('test_roi4.timesOn', test_roi4.timesOn[0])
               loop_test_novelty.addData('test_roi4.timesOff', test_roi4.timesOff[0])
            else:
               loop_test_novelty.addData('test_roi4.timesOn', "")
               loop_test_novelty.addData('test_roi4.timesOff', "")
            # using non-slip timing so subtract the expected duration of this Routine (unless ended on request)
            if routineForceEnded:
                routineTimer.reset()
            else:
                routineTimer.addTime(-2.000000)
            
            # --- Prepare to start Routine "test_feedback" ---
            continueRoutine = True
            routineForceEnded = False
            # update component parameters for each repeat
            # Run 'Begin Routine' code from codeFeedbackTest
            logging.log(level=logging.INFO, msg=f'FeedbackTestOnset_{log_msg}')
            ioServer.sendMessageEvent(text='FeedbackTestOnset')
            eyetracker.sendMessage('FeedbackTestOnset')
            eyetracker.sendMessage(outcome)
            
            polygonFeedbackTest.setOpacity(feedback_opacity)
            polygonFeedbackTest.setPos(feedback_position)
            polygonFeedbackTest.setLineColor(feedback_color)
            test_image_feedback1.setPos(position1)
            test_image_feedback1.setImage(eval(stim1))
            test_image_feedback2.setPos(position2)
            test_image_feedback2.setImage(eval(stim2))
            test_iamge_feedback3.setPos(position3)
            test_iamge_feedback3.setImage(eval(stim3))
            test_iamge_feedback4.setPos(position4)
            test_iamge_feedback4.setImage(eval(stim4))
            soundFeedbackTest.setSound(feedback_audio, secs=0.5, hamming=True)
            soundFeedbackTest.setVolume(1.0, log=False)
            # keep track of which components have finished
            test_feedbackComponents = [textScoreTest, polygonFeedbackTest, test_image_feedback1, test_image_feedback2, test_iamge_feedback3, test_iamge_feedback4, gazeCursorTestFeedback, soundFeedbackTest]
            for thisComponent in test_feedbackComponents:
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
            
            # --- Run Routine "test_feedback" ---
            while continueRoutine and routineTimer.getTime() < 0.5:
                # get current time
                t = routineTimer.getTime()
                tThisFlip = win.getFutureFlipTime(clock=routineTimer)
                tThisFlipGlobal = win.getFutureFlipTime(clock=None)
                frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
                # update/draw components on each frame
                # Run 'Each Frame' code from codeFeedbackTest
                if loop_timer.getTime() >= 2.0:
                    continueRoutine = False
                    currentLoop.finished = True
                
                # *textScoreTest* updates
                if textScoreTest.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                    # keep track of start time/frame for later
                    textScoreTest.frameNStart = frameN  # exact frame index
                    textScoreTest.tStart = t  # local t and not account for scr refresh
                    textScoreTest.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(textScoreTest, 'tStartRefresh')  # time at next scr refresh
                    textScoreTest.setAutoDraw(True)
                if textScoreTest.status == STARTED:
                    # is it time to stop? (based on global clock, using actual start)
                    if tThisFlipGlobal > textScoreTest.tStartRefresh + 0.5-frameTolerance:
                        # keep track of stop time/frame for later
                        textScoreTest.tStop = t  # not accounting for scr refresh
                        textScoreTest.frameNStop = frameN  # exact frame index
                        textScoreTest.setAutoDraw(False)
                if textScoreTest.status == STARTED:  # only update if drawing
                    textScoreTest.setText(feedback_score, log=False)
                
                # *polygonFeedbackTest* updates
                if polygonFeedbackTest.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                    # keep track of start time/frame for later
                    polygonFeedbackTest.frameNStart = frameN  # exact frame index
                    polygonFeedbackTest.tStart = t  # local t and not account for scr refresh
                    polygonFeedbackTest.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(polygonFeedbackTest, 'tStartRefresh')  # time at next scr refresh
                    polygonFeedbackTest.setAutoDraw(True)
                if polygonFeedbackTest.status == STARTED:
                    # is it time to stop? (based on global clock, using actual start)
                    if tThisFlipGlobal > polygonFeedbackTest.tStartRefresh + 0.5-frameTolerance:
                        # keep track of stop time/frame for later
                        polygonFeedbackTest.tStop = t  # not accounting for scr refresh
                        polygonFeedbackTest.frameNStop = frameN  # exact frame index
                        polygonFeedbackTest.setAutoDraw(False)
                if polygonFeedbackTest.status == STARTED:  # only update if drawing
                    polygonFeedbackTest.setSize(rectsize, log=False)
                
                # *test_image_feedback1* updates
                if test_image_feedback1.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                    # keep track of start time/frame for later
                    test_image_feedback1.frameNStart = frameN  # exact frame index
                    test_image_feedback1.tStart = t  # local t and not account for scr refresh
                    test_image_feedback1.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(test_image_feedback1, 'tStartRefresh')  # time at next scr refresh
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'test_image_feedback1.started')
                    test_image_feedback1.setAutoDraw(True)
                if test_image_feedback1.status == STARTED:
                    # is it time to stop? (based on global clock, using actual start)
                    if tThisFlipGlobal > test_image_feedback1.tStartRefresh + 0.5-frameTolerance:
                        # keep track of stop time/frame for later
                        test_image_feedback1.tStop = t  # not accounting for scr refresh
                        test_image_feedback1.frameNStop = frameN  # exact frame index
                        # add timestamp to datafile
                        thisExp.timestampOnFlip(win, 'test_image_feedback1.stopped')
                        test_image_feedback1.setAutoDraw(False)
                
                # *test_image_feedback2* updates
                if test_image_feedback2.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                    # keep track of start time/frame for later
                    test_image_feedback2.frameNStart = frameN  # exact frame index
                    test_image_feedback2.tStart = t  # local t and not account for scr refresh
                    test_image_feedback2.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(test_image_feedback2, 'tStartRefresh')  # time at next scr refresh
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'test_image_feedback2.started')
                    test_image_feedback2.setAutoDraw(True)
                if test_image_feedback2.status == STARTED:
                    # is it time to stop? (based on global clock, using actual start)
                    if tThisFlipGlobal > test_image_feedback2.tStartRefresh + 0.5-frameTolerance:
                        # keep track of stop time/frame for later
                        test_image_feedback2.tStop = t  # not accounting for scr refresh
                        test_image_feedback2.frameNStop = frameN  # exact frame index
                        # add timestamp to datafile
                        thisExp.timestampOnFlip(win, 'test_image_feedback2.stopped')
                        test_image_feedback2.setAutoDraw(False)
                
                # *test_iamge_feedback3* updates
                if test_iamge_feedback3.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                    # keep track of start time/frame for later
                    test_iamge_feedback3.frameNStart = frameN  # exact frame index
                    test_iamge_feedback3.tStart = t  # local t and not account for scr refresh
                    test_iamge_feedback3.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(test_iamge_feedback3, 'tStartRefresh')  # time at next scr refresh
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'test_iamge_feedback3.started')
                    test_iamge_feedback3.setAutoDraw(True)
                if test_iamge_feedback3.status == STARTED:
                    # is it time to stop? (based on global clock, using actual start)
                    if tThisFlipGlobal > test_iamge_feedback3.tStartRefresh + 0.5-frameTolerance:
                        # keep track of stop time/frame for later
                        test_iamge_feedback3.tStop = t  # not accounting for scr refresh
                        test_iamge_feedback3.frameNStop = frameN  # exact frame index
                        # add timestamp to datafile
                        thisExp.timestampOnFlip(win, 'test_iamge_feedback3.stopped')
                        test_iamge_feedback3.setAutoDraw(False)
                
                # *test_iamge_feedback4* updates
                if test_iamge_feedback4.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                    # keep track of start time/frame for later
                    test_iamge_feedback4.frameNStart = frameN  # exact frame index
                    test_iamge_feedback4.tStart = t  # local t and not account for scr refresh
                    test_iamge_feedback4.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(test_iamge_feedback4, 'tStartRefresh')  # time at next scr refresh
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'test_iamge_feedback4.started')
                    test_iamge_feedback4.setAutoDraw(True)
                if test_iamge_feedback4.status == STARTED:
                    # is it time to stop? (based on global clock, using actual start)
                    if tThisFlipGlobal > test_iamge_feedback4.tStartRefresh + 0.5-frameTolerance:
                        # keep track of stop time/frame for later
                        test_iamge_feedback4.tStop = t  # not accounting for scr refresh
                        test_iamge_feedback4.frameNStop = frameN  # exact frame index
                        # add timestamp to datafile
                        thisExp.timestampOnFlip(win, 'test_iamge_feedback4.stopped')
                        test_iamge_feedback4.setAutoDraw(False)
                
                # *gazeCursorTestFeedback* updates
                if gazeCursorTestFeedback.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                    # keep track of start time/frame for later
                    gazeCursorTestFeedback.frameNStart = frameN  # exact frame index
                    gazeCursorTestFeedback.tStart = t  # local t and not account for scr refresh
                    gazeCursorTestFeedback.tStartRefresh = tThisFlipGlobal  # on global time
                    win.timeOnFlip(gazeCursorTestFeedback, 'tStartRefresh')  # time at next scr refresh
                    gazeCursorTestFeedback.setAutoDraw(True)
                if gazeCursorTestFeedback.status == STARTED:
                    # is it time to stop? (based on global clock, using actual start)
                    if tThisFlipGlobal > gazeCursorTestFeedback.tStartRefresh + 0.5-frameTolerance:
                        # keep track of stop time/frame for later
                        gazeCursorTestFeedback.tStop = t  # not accounting for scr refresh
                        gazeCursorTestFeedback.frameNStop = frameN  # exact frame index
                        gazeCursorTestFeedback.setAutoDraw(False)
                if gazeCursorTestFeedback.status == STARTED:  # only update if drawing
                    gazeCursorTestFeedback.setFillColor(cursorcolor, log=False)
                    gazeCursorTestFeedback.setOpacity(0.0, log=False)
                    gazeCursorTestFeedback.setPos([eyetracker.getPos()], log=False)
                # start/stop soundFeedbackTest
                if soundFeedbackTest.status == NOT_STARTED and t >= 0-frameTolerance:
                    # keep track of start time/frame for later
                    soundFeedbackTest.frameNStart = frameN  # exact frame index
                    soundFeedbackTest.tStart = t  # local t and not account for scr refresh
                    soundFeedbackTest.tStartRefresh = tThisFlipGlobal  # on global time
                    soundFeedbackTest.play()  # start the sound (it finishes automatically)
                if soundFeedbackTest.status == STARTED:
                    # is it time to stop? (based on global clock, using actual start)
                    if tThisFlipGlobal > soundFeedbackTest.tStartRefresh + 0.5-frameTolerance:
                        # keep track of stop time/frame for later
                        soundFeedbackTest.tStop = t  # not accounting for scr refresh
                        soundFeedbackTest.frameNStop = frameN  # exact frame index
                        soundFeedbackTest.stop()
                
                # check for quit (typically the Esc key)
                if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
                    core.quit()
                
                # check if all components have finished
                if not continueRoutine:  # a component has requested a forced-end of Routine
                    routineForceEnded = True
                    break
                continueRoutine = False  # will revert to True if at least one component still running
                for thisComponent in test_feedbackComponents:
                    if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                        continueRoutine = True
                        break  # at least one component has not yet finished
                
                # refresh the screen
                if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                    win.flip()
            
            # --- Ending Routine "test_feedback" ---
            for thisComponent in test_feedbackComponents:
                if hasattr(thisComponent, "setAutoDraw"):
                    thisComponent.setAutoDraw(False)
            soundFeedbackTest.stop()  # ensure sound has stopped at end of routine
            # using non-slip timing so subtract the expected duration of this Routine (unless ended on request)
            if routineForceEnded:
                routineTimer.reset()
            else:
                routineTimer.addTime(-0.500000)
            thisExp.nextEntry()
            
        # completed 100.0 repeats of 'loop_test_novelty'
        
        
        # --- Prepare to start Routine "crossEnd" ---
        continueRoutine = True
        routineForceEnded = False
        # update component parameters for each repeat
        # Run 'Begin Routine' code from code_end
        eyetracker.setRecordingState(False)
        jittered_duration_cross = 3.0 + random()
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
                fixcrossE.setAutoDraw(True)
            if fixcrossE.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > fixcrossE.tStartRefresh + jittered_duration_cross-frameTolerance:
                    # keep track of stop time/frame for later
                    fixcrossE.tStop = t  # not accounting for scr refresh
                    fixcrossE.frameNStop = frameN  # exact frame index
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
        
    # completed 1.0 repeats of 'testtrials_novelty'
    
    
    # set up handler to look after randomisation of conditions etc
    ratingtrials3 = data.TrialHandler(nReps=1.0, method='random', 
        extraInfo=expInfo, originPath=-1,
        trialList=data.importConditions(stimFileNew),
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
        logging.log(level=logging.INFO, msg=f'StimRating_{stim}')
        print("Rating of Stimulus: %s"%(stim))
        win.mouseVisible = True
        imageRating.setPos((0, 0.2))
        imageRating.setImage(eval(stim))
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
