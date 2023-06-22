#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This experiment was created using PsychoPy3 Experiment Builder (v2023.1.1),
    on Juni 22, 2023, at 11:57
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

# Run 'Before Experiment' code from code_trial
feedback_color = "grey"
feedback_opacity = 0
feedback_audio = "audio/silence.wav"
cursorcolor="white"
digitimer_msg = 0
biopac_feedback_msg = 0

rectsize = (0.2, 0.2)
imagesize_test = [0.2, 0.2]
imagesize_rating = [0.2, 0.2]

score = 0
points = 0


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
    originPath='C:\\Users\\sag22id\\Documents\\Projects\\GCA\\gaze_avoidance\\gca_lastrun.py',
    savePickle=True, saveWideText=True,
    dataFileName=filename)
# save a log file for detail verbose info
logFile = logging.LogFile(filename+'.log', level=logging.EXP)
logging.console.setLevel(logging.WARNING)  # this outputs to the screen, not a file

endExpNow = False  # flag for 'escape' or other condition => quit the exp
frameTolerance = 0.001  # how close to onset before 'same' frame

# Start Code - component code to be run after the window creation

# --- Setup the Window ---
win = visual.Window(
    size=[2194, 1234], fullscr=True, screen=0, 
    winType='pyglet', allowStencil=False,
    monitor='officeMonitor', color=[0,0,0], colorSpace='rgb',
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
ioServer = io.launchHubServer(window=win, **ioConfig)
eyetracker = ioServer.getDevice('tracker')

# create a default keyboard (e.g. to check for escape)
defaultKeyboard = keyboard.Keyboard(backend='iohub')

# --- Initialize components for Routine "welcome" ---
welcome_msg = visual.TextStim(win=win, name='welcome_msg',
    text='Sehr geehrter Teilnehmer, sehr geehrte Teilnehmerin,\nvielen Dank für Ihre Teilnahme an unserem Experiment.\n\nBitte drücken Sie die Leertaste um zu starten.',
    font='Open Sans',
    pos=(0, 0), height=0.03, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=0.0);
space_resp = keyboard.Keyboard()

# --- Initialize components for Routine "painRating" ---
textRate = visual.TextStim(win=win, name='textRate',
    text='Wie unangenehm finden Sie den präsentierten Reiz?',
    font='Open Sans',
    pos=(0, 0.3), height=0.05, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=0.0);
text_noPain = visual.TextStim(win=win, name='text_noPain',
    text='kein\nSchmerz',
    font='Open Sans',
    pos=(-0.5, 0.1), height=0.03, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-1.0);
text_highPain = visual.TextStim(win=win, name='text_highPain',
    text='sehr starker\nSchmerz',
    font='Open Sans',
    pos=(0.5, 0.1), height=0.03, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-2.0);
sliderPain = visual.Slider(win=win, name='sliderPain',
    startValue=None, size=(1.0, 0.05), pos=(0, 0), units=win.units,
    labels=(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), ticks=(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), granularity=1.0,
    style='rating', styleTweaks=(), opacity=None,
    labelColor='White', markerColor='Red', lineColor='White', colorSpace='rgb',
    font='Open Sans', labelHeight=0.02,
    flip=False, ori=0.0, depth=-3, readOnly=False)
PainRating_Shock1 = parallel.ParallelPort(address='0x0378')
PainRating_Shock2 = parallel.ParallelPort(address='0x0378')
PainRating_Shock3 = parallel.ParallelPort(address='0x0378')
textSpacePain = visual.TextStim(win=win, name='textSpacePain',
    text='Bitte drücken Sie die Leertaste, um die Antwort zu bestätigen.',
    font='Open Sans',
    pos=(0, -0.3), height=0.05, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-7.0);
spacePain = keyboard.Keyboard()

# --- Initialize components for Routine "startRecord" ---
startRecording = hardware.eyetracker.EyetrackerControl(
    tracker=eyetracker,
    actionType='Start Only'
)

# --- Initialize components for Routine "cross" ---
fixateStart = visual.ShapeStim(
    win=win, name='fixateStart', vertices='cross',
    size=(0.03, 0.03),
    ori=0.0, pos=(0, 0), anchor='center',
    lineWidth=1.0,     colorSpace='rgb',  lineColor='white', fillColor='white',
    opacity=None, depth=0.0, interpolate=True)

# --- Initialize components for Routine "blank" ---
blank_screen = visual.TextStim(win=win, name='blank_screen',
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
    opacity=0.0, depth=-3.0, interpolate=True)
Biopac_ImageOnset = parallel.ParallelPort(address='0x0378')

# --- Initialize components for Routine "feedback" ---
textPoints = visual.TextStim(win=win, name='textPoints',
    text='',
    font='Open Sans',
    pos=(0, 0.3), height=0.06, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=0.0);
textScore = visual.TextStim(win=win, name='textScore',
    text='',
    font='Open Sans',
    pos=(0, 0.24), height=0.06, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-1.0);
polygon = visual.Rect(
    win=win, name='polygon',
    width=[1.0, 1.0][0], height=[1.0, 1.0][1],
    ori=0.0, pos=[0,0], anchor='center',
    lineWidth=1.0,     colorSpace='rgb',  lineColor=feedback_color, fillColor='white',
    opacity=1.0, depth=-2.0, interpolate=True)
imageFeedback = visual.ImageStim(
    win=win,
    name='imageFeedback', 
    image='default.png', mask=None, anchor='center',
    ori=0.0, pos=[0,0], size=None,
    color=[1,1,1], colorSpace='rgb', opacity=None,
    flipHoriz=False, flipVert=False,
    texRes=128.0, interpolate=True, depth=-3.0)
gazeCursor_Feedback = visual.ShapeStim(
    win=win, name='gazeCursor_Feedback', vertices='star7',
    size=(0.02, 0.02),
    ori=0.0, pos=[0,0], anchor='center',
    lineWidth=1.0,     colorSpace='rgb',  lineColor=None, fillColor='white',
    opacity=1.0, depth=-4.0, interpolate=True)
soundFeedback = sound.Sound('A', secs=0.5, stereo=True, hamming=True,
    name='soundFeedback')
soundFeedback.setVolume(1.0)
Digitimer_Shock1 = parallel.ParallelPort(address='0x0378')
Digitimer_Shock2 = parallel.ParallelPort(address='0x0378')
Digitimer_Shock3 = parallel.ParallelPort(address='0x0378')
Biopac_Feedback = parallel.ParallelPort(address='0x0378')

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
text_unpleasant = visual.TextStim(win=win, name='text_unpleasant',
    text='sehr \nunwohl',
    font='Open Sans',
    pos=(-0.5, -0.1), height=0.03, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-2.0);
text_pleasant = visual.TextStim(win=win, name='text_pleasant',
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

# --- Initialize components for Routine "cross" ---
fixateStart = visual.ShapeStim(
    win=win, name='fixateStart', vertices='cross',
    size=(0.03, 0.03),
    ori=0.0, pos=(0, 0), anchor='center',
    lineWidth=1.0,     colorSpace='rgb',  lineColor='white', fillColor='white',
    opacity=None, depth=0.0, interpolate=True)

# --- Initialize components for Routine "testtrial" ---
imageTest = visual.ImageStim(
    win=win,
    name='imageTest', 
    image='default.png', mask=None, anchor='center',
    ori=0.0, pos=[0,0], size=1.0,
    color=[1,1,1], colorSpace='rgb', opacity=None,
    flipHoriz=False, flipVert=False,
    texRes=128.0, interpolate=True, depth=0.0)
Biopac_TestImageOnset = parallel.ParallelPort(address='0x0378')

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
text_unpleasant = visual.TextStim(win=win, name='text_unpleasant',
    text='sehr \nunwohl',
    font='Open Sans',
    pos=(-0.5, -0.1), height=0.03, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-2.0);
text_pleasant = visual.TextStim(win=win, name='text_pleasant',
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
end_msg = visual.TextStim(win=win, name='end_msg',
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
space_resp.keys = []
space_resp.rt = []
_space_resp_allKeys = []
# keep track of which components have finished
welcomeComponents = [welcome_msg, space_resp]
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
    
    # *welcome_msg* updates
    
    # if welcome_msg is starting this frame...
    if welcome_msg.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        # keep track of start time/frame for later
        welcome_msg.frameNStart = frameN  # exact frame index
        welcome_msg.tStart = t  # local t and not account for scr refresh
        welcome_msg.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(welcome_msg, 'tStartRefresh')  # time at next scr refresh
        # add timestamp to datafile
        thisExp.timestampOnFlip(win, 'welcome_msg.started')
        # update status
        welcome_msg.status = STARTED
        welcome_msg.setAutoDraw(True)
    
    # if welcome_msg is active this frame...
    if welcome_msg.status == STARTED:
        # update params
        pass
    
    # *space_resp* updates
    waitOnFlip = False
    
    # if space_resp is starting this frame...
    if space_resp.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        # keep track of start time/frame for later
        space_resp.frameNStart = frameN  # exact frame index
        space_resp.tStart = t  # local t and not account for scr refresh
        space_resp.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(space_resp, 'tStartRefresh')  # time at next scr refresh
        # add timestamp to datafile
        thisExp.timestampOnFlip(win, 'space_resp.started')
        # update status
        space_resp.status = STARTED
        # keyboard checking is just starting
        waitOnFlip = True
        win.callOnFlip(space_resp.clock.reset)  # t=0 on next screen flip
        win.callOnFlip(space_resp.clearEvents, eventType='keyboard')  # clear events on next screen flip
    if space_resp.status == STARTED and not waitOnFlip:
        theseKeys = space_resp.getKeys(keyList=['space', 'enter'], waitRelease=False)
        _space_resp_allKeys.extend(theseKeys)
        if len(_space_resp_allKeys):
            space_resp.keys = _space_resp_allKeys[-1].name  # just the last key pressed
            space_resp.rt = _space_resp_allKeys[-1].rt
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
if space_resp.keys in ['', [], None]:  # No response was made
    space_resp.keys = None
thisExp.addData('space_resp.keys',space_resp.keys)
if space_resp.keys != None:  # we had a response
    thisExp.addData('space_resp.rt', space_resp.rt)
thisExp.nextEntry()
# the Routine "welcome" was not non-slip safe, so reset the non-slip timer
routineTimer.reset()

# --- Prepare to start Routine "painRating" ---
continueRoutine = True
# update component parameters for each repeat
sliderPain.reset()
spacePain.keys = []
spacePain.rt = []
_spacePain_allKeys = []
# keep track of which components have finished
painRatingComponents = [textRate, text_noPain, text_highPain, sliderPain, PainRating_Shock1, PainRating_Shock2, PainRating_Shock3, textSpacePain, spacePain]
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
    # *PainRating_Shock1* updates
    
    # if PainRating_Shock1 is starting this frame...
    if PainRating_Shock1.status == NOT_STARTED and t >= 2-frameTolerance:
        # keep track of start time/frame for later
        PainRating_Shock1.frameNStart = frameN  # exact frame index
        PainRating_Shock1.tStart = t  # local t and not account for scr refresh
        PainRating_Shock1.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(PainRating_Shock1, 'tStartRefresh')  # time at next scr refresh
        # add timestamp to datafile
        thisExp.addData('PainRating_Shock1.started', t)
        # update status
        PainRating_Shock1.status = STARTED
        PainRating_Shock1.status = STARTED
        win.callOnFlip(PainRating_Shock1.setData, int(1))
    
    # if PainRating_Shock1 is stopping this frame...
    if PainRating_Shock1.status == STARTED:
        # is it time to stop? (based on global clock, using actual start)
        if tThisFlipGlobal > PainRating_Shock1.tStartRefresh + 0.002-frameTolerance:
            # keep track of stop time/frame for later
            PainRating_Shock1.tStop = t  # not accounting for scr refresh
            PainRating_Shock1.frameNStop = frameN  # exact frame index
            # add timestamp to datafile
            thisExp.addData('PainRating_Shock1.stopped', t)
            # update status
            PainRating_Shock1.status = FINISHED
            win.callOnFlip(PainRating_Shock1.setData, int(0))
    # *PainRating_Shock2* updates
    
    # if PainRating_Shock2 is starting this frame...
    if PainRating_Shock2.status == NOT_STARTED and t >= 2.05-frameTolerance:
        # keep track of start time/frame for later
        PainRating_Shock2.frameNStart = frameN  # exact frame index
        PainRating_Shock2.tStart = t  # local t and not account for scr refresh
        PainRating_Shock2.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(PainRating_Shock2, 'tStartRefresh')  # time at next scr refresh
        # add timestamp to datafile
        thisExp.addData('PainRating_Shock2.started', t)
        # update status
        PainRating_Shock2.status = STARTED
        PainRating_Shock2.status = STARTED
        win.callOnFlip(PainRating_Shock2.setData, int(1))
    
    # if PainRating_Shock2 is stopping this frame...
    if PainRating_Shock2.status == STARTED:
        # is it time to stop? (based on global clock, using actual start)
        if tThisFlipGlobal > PainRating_Shock2.tStartRefresh + 0.002-frameTolerance:
            # keep track of stop time/frame for later
            PainRating_Shock2.tStop = t  # not accounting for scr refresh
            PainRating_Shock2.frameNStop = frameN  # exact frame index
            # add timestamp to datafile
            thisExp.addData('PainRating_Shock2.stopped', t)
            # update status
            PainRating_Shock2.status = FINISHED
            win.callOnFlip(PainRating_Shock2.setData, int(0))
    # *PainRating_Shock3* updates
    
    # if PainRating_Shock3 is starting this frame...
    if PainRating_Shock3.status == NOT_STARTED and t >= 2.1-frameTolerance:
        # keep track of start time/frame for later
        PainRating_Shock3.frameNStart = frameN  # exact frame index
        PainRating_Shock3.tStart = t  # local t and not account for scr refresh
        PainRating_Shock3.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(PainRating_Shock3, 'tStartRefresh')  # time at next scr refresh
        # add timestamp to datafile
        thisExp.addData('PainRating_Shock3.started', t)
        # update status
        PainRating_Shock3.status = STARTED
        PainRating_Shock3.status = STARTED
        win.callOnFlip(PainRating_Shock3.setData, int(1))
    
    # if PainRating_Shock3 is stopping this frame...
    if PainRating_Shock3.status == STARTED:
        # is it time to stop? (based on global clock, using actual start)
        if tThisFlipGlobal > PainRating_Shock3.tStartRefresh + 0.002-frameTolerance:
            # keep track of stop time/frame for later
            PainRating_Shock3.tStop = t  # not accounting for scr refresh
            PainRating_Shock3.frameNStop = frameN  # exact frame index
            # add timestamp to datafile
            thisExp.addData('PainRating_Shock3.stopped', t)
            # update status
            PainRating_Shock3.status = FINISHED
            win.callOnFlip(PainRating_Shock3.setData, int(0))
    
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
thisExp.addData('sliderPain.response', sliderPain.getRating())
thisExp.addData('sliderPain.rt', sliderPain.getRT())
if PainRating_Shock1.status == STARTED:
    win.callOnFlip(PainRating_Shock1.setData, int(0))
if PainRating_Shock2.status == STARTED:
    win.callOnFlip(PainRating_Shock2.setData, int(0))
if PainRating_Shock3.status == STARTED:
    win.callOnFlip(PainRating_Shock3.setData, int(0))
# check responses
if spacePain.keys in ['', [], None]:  # No response was made
    spacePain.keys = None
thisExp.addData('spacePain.keys',spacePain.keys)
if spacePain.keys != None:  # we had a response
    thisExp.addData('spacePain.rt', spacePain.rt)
thisExp.nextEntry()
# the Routine "painRating" was not non-slip safe, so reset the non-slip timer
routineTimer.reset()

# --- Prepare to start Routine "startRecord" ---
continueRoutine = True
# update component parameters for each repeat
# keep track of which components have finished
startRecordComponents = [startRecording]
for thisComponent in startRecordComponents:
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

# --- Run Routine "startRecord" ---
routineForceEnded = not continueRoutine
while continueRoutine:
    # get current time
    t = routineTimer.getTime()
    tThisFlip = win.getFutureFlipTime(clock=routineTimer)
    tThisFlipGlobal = win.getFutureFlipTime(clock=None)
    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
    # update/draw components on each frame
    # *startRecording* updates
    
    # if startRecording is starting this frame...
    if startRecording.status == NOT_STARTED and t >= 0.0-frameTolerance:
        # keep track of start time/frame for later
        startRecording.frameNStart = frameN  # exact frame index
        startRecording.tStart = t  # local t and not account for scr refresh
        startRecording.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(startRecording, 'tStartRefresh')  # time at next scr refresh
        # add timestamp to datafile
        thisExp.addData('startRecording.started', t)
        # update status
        startRecording.status = STARTED
    
    # if startRecording is stopping this frame...
    if startRecording.status == STARTED:
        # is it time to stop? (based on global clock, using actual start)
        if tThisFlipGlobal > startRecording.tStartRefresh + 0-frameTolerance:
            # keep track of stop time/frame for later
            startRecording.tStop = t  # not accounting for scr refresh
            startRecording.frameNStop = frameN  # exact frame index
            # add timestamp to datafile
            thisExp.addData('startRecording.stopped', t)
            # update status
            startRecording.status = FINISHED
    
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
    for thisComponent in startRecordComponents:
        if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
            continueRoutine = True
            break  # at least one component has not yet finished
    
    # refresh the screen
    if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
        win.flip()

# --- Ending Routine "startRecord" ---
for thisComponent in startRecordComponents:
    if hasattr(thisComponent, "setAutoDraw"):
        thisComponent.setAutoDraw(False)
# make sure the eyetracker recording stops
if startRecording.status != FINISHED:
    startRecording.status = FINISHED
# the Routine "startRecord" was not non-slip safe, so reset the non-slip timer
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
    trials = data.TrialHandler(nReps=1.0, method='random', 
        extraInfo=expInfo, originPath=-1,
        trialList=data.importConditions(posFile, selection=random(4)*16),
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
        # keep track of which components have finished
        crossComponents = [fixateStart]
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
        while continueRoutine:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # *fixateStart* updates
            
            # if fixateStart is starting this frame...
            if fixateStart.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                fixateStart.frameNStart = frameN  # exact frame index
                fixateStart.tStart = t  # local t and not account for scr refresh
                fixateStart.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(fixateStart, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'fixateStart.started')
                # update status
                fixateStart.status = STARTED
                fixateStart.setAutoDraw(True)
            
            # if fixateStart is active this frame...
            if fixateStart.status == STARTED:
                # update params
                pass
            
            # if fixateStart is stopping this frame...
            if fixateStart.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > fixateStart.tStartRefresh + 4 + random()-frameTolerance:
                    # keep track of stop time/frame for later
                    fixateStart.tStop = t  # not accounting for scr refresh
                    fixateStart.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'fixateStart.stopped')
                    # update status
                    fixateStart.status = FINISHED
                    fixateStart.setAutoDraw(False)
            
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
        # the Routine "cross" was not non-slip safe, so reset the non-slip timer
        routineTimer.reset()
        
        # --- Prepare to start Routine "blank" ---
        continueRoutine = True
        # update component parameters for each repeat
        # keep track of which components have finished
        blankComponents = [blank_screen]
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
            
            # *blank_screen* updates
            
            # if blank_screen is starting this frame...
            if blank_screen.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                blank_screen.frameNStart = frameN  # exact frame index
                blank_screen.tStart = t  # local t and not account for scr refresh
                blank_screen.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(blank_screen, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'blank_screen.started')
                # update status
                blank_screen.status = STARTED
                blank_screen.setAutoDraw(True)
            
            # if blank_screen is active this frame...
            if blank_screen.status == STARTED:
                # update params
                pass
            
            # if blank_screen is stopping this frame...
            if blank_screen.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > blank_screen.tStartRefresh + 1.0 + random()-frameTolerance:
                    # keep track of stop time/frame for later
                    blank_screen.tStop = t  # not accounting for scr refresh
                    blank_screen.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'blank_screen.stopped')
                    # update status
                    blank_screen.status = FINISHED
                    blank_screen.setAutoDraw(False)
            
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
        # Run 'Begin Routine' code from code_trial
        looked_at = False
        cursorcolor="white"
        image.setPos(position)
        image.setImage(eval(trialtype))
        # clear any previous roi data
        roi.reset()
        # keep track of which components have finished
        trialComponents = [image, roi, gazeCursor, Biopac_ImageOnset]
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
            # Run 'Each Frame' code from code_trial
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
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'roi.started')
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
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'roi.stopped')
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
            # *Biopac_ImageOnset* updates
            
            # if Biopac_ImageOnset is starting this frame...
            if Biopac_ImageOnset.status == NOT_STARTED and t >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                Biopac_ImageOnset.frameNStart = frameN  # exact frame index
                Biopac_ImageOnset.tStart = t  # local t and not account for scr refresh
                Biopac_ImageOnset.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(Biopac_ImageOnset, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.addData('Biopac_ImageOnset.started', t)
                # update status
                Biopac_ImageOnset.status = STARTED
                Biopac_ImageOnset.status = STARTED
                win.callOnFlip(Biopac_ImageOnset.setData, int(1))
            
            # if Biopac_ImageOnset is stopping this frame...
            if Biopac_ImageOnset.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > Biopac_ImageOnset.tStartRefresh + 1.0-frameTolerance:
                    # keep track of stop time/frame for later
                    Biopac_ImageOnset.tStop = t  # not accounting for scr refresh
                    Biopac_ImageOnset.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.addData('Biopac_ImageOnset.stopped', t)
                    # update status
                    Biopac_ImageOnset.status = FINISHED
                    win.callOnFlip(Biopac_ImageOnset.setData, int(0))
            
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
        # Run 'End Routine' code from code_trial
        if looked_at & ("plus" in trialtype):
            feedback_color = "red"
            feedback_opacity = 1
            cursorcolor="red"
            feedback_audio = "audio/error.wav"
            digitimer_msg = 128
            biopac_feedback_msg = 2  # Error
            points = 0
            feedback_points = ""
            feedback_score = ""
        elif looked_at & ("minus" in trialtype):
            feedback_color = "green"
            feedback_opacity = 1
            cursorcolor="green"
            feedback_audio = "audio/win.wav"
            digitimer_msg = 0
            biopac_feedback_msg = 3  # Reward
            points = 5
            score += points
            feedback_points = f"+ {points}"
            feedback_score = f"Score: {score}"
        else:
            feedback_color = "grey"
            feedback_opacity = 0
            feedback_audio = "audio/silence.wav"
            digitimer_msg = 0
            biopac_feedback_msg = 4  # None
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
        if Biopac_ImageOnset.status == STARTED:
            win.callOnFlip(Biopac_ImageOnset.setData, int(0))
        # using non-slip timing so subtract the expected duration of this Routine (unless ended on request)
        if routineForceEnded:
            routineTimer.reset()
        else:
            routineTimer.addTime(-5.000000)
        
        # --- Prepare to start Routine "feedback" ---
        continueRoutine = True
        # update component parameters for each repeat
        polygon.setFillColor(feedback_color)
        polygon.setOpacity(feedback_opacity)
        polygon.setPos(position)
        imageFeedback.setPos(position)
        imageFeedback.setImage(eval(trialtype))
        gazeCursor_Feedback.setOpacity(0.0)
        soundFeedback.setSound(feedback_audio, secs=0.5, hamming=True)
        soundFeedback.setVolume(1.0, log=False)
        # keep track of which components have finished
        feedbackComponents = [textPoints, textScore, polygon, imageFeedback, gazeCursor_Feedback, soundFeedback, Digitimer_Shock1, Digitimer_Shock2, Digitimer_Shock3, Biopac_Feedback]
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
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'textPoints.started')
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
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'textPoints.stopped')
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
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'textScore.started')
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
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'textScore.stopped')
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
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'polygon.started')
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
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'polygon.stopped')
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
            if soundFeedback.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                # keep track of start time/frame for later
                soundFeedback.frameNStart = frameN  # exact frame index
                soundFeedback.tStart = t  # local t and not account for scr refresh
                soundFeedback.tStartRefresh = tThisFlipGlobal  # on global time
                # add timestamp to datafile
                thisExp.addData('soundFeedback.started', tThisFlipGlobal)
                # update status
                soundFeedback.status = STARTED
                soundFeedback.play(when=win)  # sync with win flip
            
            # if soundFeedback is stopping this frame...
            if soundFeedback.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > soundFeedback.tStartRefresh + 0.5-frameTolerance:
                    # keep track of stop time/frame for later
                    soundFeedback.tStop = t  # not accounting for scr refresh
                    soundFeedback.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'soundFeedback.stopped')
                    # update status
                    soundFeedback.status = FINISHED
                    soundFeedback.stop()
            # *Digitimer_Shock1* updates
            
            # if Digitimer_Shock1 is starting this frame...
            if Digitimer_Shock1.status == NOT_STARTED and t >= 0-frameTolerance:
                # keep track of start time/frame for later
                Digitimer_Shock1.frameNStart = frameN  # exact frame index
                Digitimer_Shock1.tStart = t  # local t and not account for scr refresh
                Digitimer_Shock1.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(Digitimer_Shock1, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.addData('Digitimer_Shock1.started', t)
                # update status
                Digitimer_Shock1.status = STARTED
                Digitimer_Shock1.status = STARTED
                win.callOnFlip(Digitimer_Shock1.setData, int(digitimer_msg))
            
            # if Digitimer_Shock1 is stopping this frame...
            if Digitimer_Shock1.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > Digitimer_Shock1.tStartRefresh + 0.002-frameTolerance:
                    # keep track of stop time/frame for later
                    Digitimer_Shock1.tStop = t  # not accounting for scr refresh
                    Digitimer_Shock1.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.addData('Digitimer_Shock1.stopped', t)
                    # update status
                    Digitimer_Shock1.status = FINISHED
                    win.callOnFlip(Digitimer_Shock1.setData, int(0))
            # *Digitimer_Shock2* updates
            
            # if Digitimer_Shock2 is starting this frame...
            if Digitimer_Shock2.status == NOT_STARTED and t >= 0.05-frameTolerance:
                # keep track of start time/frame for later
                Digitimer_Shock2.frameNStart = frameN  # exact frame index
                Digitimer_Shock2.tStart = t  # local t and not account for scr refresh
                Digitimer_Shock2.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(Digitimer_Shock2, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.addData('Digitimer_Shock2.started', t)
                # update status
                Digitimer_Shock2.status = STARTED
                Digitimer_Shock2.status = STARTED
                win.callOnFlip(Digitimer_Shock2.setData, int(digitimer_msg))
            
            # if Digitimer_Shock2 is stopping this frame...
            if Digitimer_Shock2.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > Digitimer_Shock2.tStartRefresh + 0.002-frameTolerance:
                    # keep track of stop time/frame for later
                    Digitimer_Shock2.tStop = t  # not accounting for scr refresh
                    Digitimer_Shock2.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.addData('Digitimer_Shock2.stopped', t)
                    # update status
                    Digitimer_Shock2.status = FINISHED
                    win.callOnFlip(Digitimer_Shock2.setData, int(0))
            # *Digitimer_Shock3* updates
            
            # if Digitimer_Shock3 is starting this frame...
            if Digitimer_Shock3.status == NOT_STARTED and t >= 0.1-frameTolerance:
                # keep track of start time/frame for later
                Digitimer_Shock3.frameNStart = frameN  # exact frame index
                Digitimer_Shock3.tStart = t  # local t and not account for scr refresh
                Digitimer_Shock3.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(Digitimer_Shock3, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.addData('Digitimer_Shock3.started', t)
                # update status
                Digitimer_Shock3.status = STARTED
                Digitimer_Shock3.status = STARTED
                win.callOnFlip(Digitimer_Shock3.setData, int(digitimer_msg))
            
            # if Digitimer_Shock3 is stopping this frame...
            if Digitimer_Shock3.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > Digitimer_Shock3.tStartRefresh + 0.002-frameTolerance:
                    # keep track of stop time/frame for later
                    Digitimer_Shock3.tStop = t  # not accounting for scr refresh
                    Digitimer_Shock3.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.addData('Digitimer_Shock3.stopped', t)
                    # update status
                    Digitimer_Shock3.status = FINISHED
                    win.callOnFlip(Digitimer_Shock3.setData, int(0))
            # *Biopac_Feedback* updates
            
            # if Biopac_Feedback is starting this frame...
            if Biopac_Feedback.status == NOT_STARTED and t >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                Biopac_Feedback.frameNStart = frameN  # exact frame index
                Biopac_Feedback.tStart = t  # local t and not account for scr refresh
                Biopac_Feedback.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(Biopac_Feedback, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.addData('Biopac_Feedback.started', t)
                # update status
                Biopac_Feedback.status = STARTED
                Biopac_Feedback.status = STARTED
                win.callOnFlip(Biopac_Feedback.setData, int(biopac_feedback_msg))
            
            # if Biopac_Feedback is stopping this frame...
            if Biopac_Feedback.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > Biopac_Feedback.tStartRefresh + 1.0-frameTolerance:
                    # keep track of stop time/frame for later
                    Biopac_Feedback.tStop = t  # not accounting for scr refresh
                    Biopac_Feedback.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.addData('Biopac_Feedback.stopped', t)
                    # update status
                    Biopac_Feedback.status = FINISHED
                    win.callOnFlip(Biopac_Feedback.setData, int(0))
            
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
        if Digitimer_Shock1.status == STARTED:
            win.callOnFlip(Digitimer_Shock1.setData, int(0))
        if Digitimer_Shock2.status == STARTED:
            win.callOnFlip(Digitimer_Shock2.setData, int(0))
        if Digitimer_Shock3.status == STARTED:
            win.callOnFlip(Digitimer_Shock3.setData, int(0))
        if Biopac_Feedback.status == STARTED:
            win.callOnFlip(Biopac_Feedback.setData, int(0))
        # using non-slip timing so subtract the expected duration of this Routine (unless ended on request)
        if routineForceEnded:
            routineTimer.reset()
        else:
            routineTimer.addTime(-3.000000)
        thisExp.nextEntry()
        
    # completed 1.0 repeats of 'trials'
    
    
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
        # keep track of which components have finished
        stimRatingComponents = [textRateStim, imageRating, text_unpleasant, text_pleasant, sliderStim, textSpaceStim, spaceStim]
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
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'textRateStim.started')
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
            
            # *text_unpleasant* updates
            
            # if text_unpleasant is starting this frame...
            if text_unpleasant.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                text_unpleasant.frameNStart = frameN  # exact frame index
                text_unpleasant.tStart = t  # local t and not account for scr refresh
                text_unpleasant.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(text_unpleasant, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'text_unpleasant.started')
                # update status
                text_unpleasant.status = STARTED
                text_unpleasant.setAutoDraw(True)
            
            # if text_unpleasant is active this frame...
            if text_unpleasant.status == STARTED:
                # update params
                pass
            
            # *text_pleasant* updates
            
            # if text_pleasant is starting this frame...
            if text_pleasant.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                text_pleasant.frameNStart = frameN  # exact frame index
                text_pleasant.tStart = t  # local t and not account for scr refresh
                text_pleasant.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(text_pleasant, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'text_pleasant.started')
                # update status
                text_pleasant.status = STARTED
                text_pleasant.setAutoDraw(True)
            
            # if text_pleasant is active this frame...
            if text_pleasant.status == STARTED:
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
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'textSpaceStim.started')
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
    
    
    # set up handler to look after randomisation of conditions etc
    testtrials = data.TrialHandler(nReps=2.0, method='random', 
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
        # keep track of which components have finished
        crossComponents = [fixateStart]
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
        while continueRoutine:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # *fixateStart* updates
            
            # if fixateStart is starting this frame...
            if fixateStart.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                fixateStart.frameNStart = frameN  # exact frame index
                fixateStart.tStart = t  # local t and not account for scr refresh
                fixateStart.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(fixateStart, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'fixateStart.started')
                # update status
                fixateStart.status = STARTED
                fixateStart.setAutoDraw(True)
            
            # if fixateStart is active this frame...
            if fixateStart.status == STARTED:
                # update params
                pass
            
            # if fixateStart is stopping this frame...
            if fixateStart.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > fixateStart.tStartRefresh + 4 + random()-frameTolerance:
                    # keep track of stop time/frame for later
                    fixateStart.tStop = t  # not accounting for scr refresh
                    fixateStart.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'fixateStart.stopped')
                    # update status
                    fixateStart.status = FINISHED
                    fixateStart.setAutoDraw(False)
            
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
        # the Routine "cross" was not non-slip safe, so reset the non-slip timer
        routineTimer.reset()
        
        # --- Prepare to start Routine "testtrial" ---
        continueRoutine = True
        # update component parameters for each repeat
        imageTest.setPos((0, 0))
        imageTest.setImage(eval(stimtype))
        # keep track of which components have finished
        testtrialComponents = [imageTest, Biopac_TestImageOnset]
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
        while continueRoutine and routineTimer.getTime() < 6.0:
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
                if tThisFlipGlobal > imageTest.tStartRefresh + 6-frameTolerance:
                    # keep track of stop time/frame for later
                    imageTest.tStop = t  # not accounting for scr refresh
                    imageTest.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'imageTest.stopped')
                    # update status
                    imageTest.status = FINISHED
                    imageTest.setAutoDraw(False)
            # *Biopac_TestImageOnset* updates
            
            # if Biopac_TestImageOnset is starting this frame...
            if Biopac_TestImageOnset.status == NOT_STARTED and t >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                Biopac_TestImageOnset.frameNStart = frameN  # exact frame index
                Biopac_TestImageOnset.tStart = t  # local t and not account for scr refresh
                Biopac_TestImageOnset.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(Biopac_TestImageOnset, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.addData('Biopac_TestImageOnset.started', t)
                # update status
                Biopac_TestImageOnset.status = STARTED
                Biopac_TestImageOnset.status = STARTED
                win.callOnFlip(Biopac_TestImageOnset.setData, int(5))
            
            # if Biopac_TestImageOnset is stopping this frame...
            if Biopac_TestImageOnset.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > Biopac_TestImageOnset.tStartRefresh + 1.0-frameTolerance:
                    # keep track of stop time/frame for later
                    Biopac_TestImageOnset.tStop = t  # not accounting for scr refresh
                    Biopac_TestImageOnset.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.addData('Biopac_TestImageOnset.stopped', t)
                    # update status
                    Biopac_TestImageOnset.status = FINISHED
                    win.callOnFlip(Biopac_TestImageOnset.setData, int(0))
            
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
        if Biopac_TestImageOnset.status == STARTED:
            win.callOnFlip(Biopac_TestImageOnset.setData, int(0))
        # using non-slip timing so subtract the expected duration of this Routine (unless ended on request)
        if routineForceEnded:
            routineTimer.reset()
        else:
            routineTimer.addTime(-6.000000)
        thisExp.nextEntry()
        
    # completed 2.0 repeats of 'testtrials'
    
    
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
        # keep track of which components have finished
        stimRatingComponents = [textRateStim, imageRating, text_unpleasant, text_pleasant, sliderStim, textSpaceStim, spaceStim]
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
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'textRateStim.started')
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
            
            # *text_unpleasant* updates
            
            # if text_unpleasant is starting this frame...
            if text_unpleasant.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                text_unpleasant.frameNStart = frameN  # exact frame index
                text_unpleasant.tStart = t  # local t and not account for scr refresh
                text_unpleasant.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(text_unpleasant, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'text_unpleasant.started')
                # update status
                text_unpleasant.status = STARTED
                text_unpleasant.setAutoDraw(True)
            
            # if text_unpleasant is active this frame...
            if text_unpleasant.status == STARTED:
                # update params
                pass
            
            # *text_pleasant* updates
            
            # if text_pleasant is starting this frame...
            if text_pleasant.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                text_pleasant.frameNStart = frameN  # exact frame index
                text_pleasant.tStart = t  # local t and not account for scr refresh
                text_pleasant.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(text_pleasant, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'text_pleasant.started')
                # update status
                text_pleasant.status = STARTED
                text_pleasant.setAutoDraw(True)
            
            # if text_pleasant is active this frame...
            if text_pleasant.status == STARTED:
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
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'textSpaceStim.started')
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
endComponents = [end_msg]
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
while continueRoutine:
    # get current time
    t = routineTimer.getTime()
    tThisFlip = win.getFutureFlipTime(clock=routineTimer)
    tThisFlipGlobal = win.getFutureFlipTime(clock=None)
    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
    # update/draw components on each frame
    
    # *end_msg* updates
    
    # if end_msg is starting this frame...
    if end_msg.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        # keep track of start time/frame for later
        end_msg.frameNStart = frameN  # exact frame index
        end_msg.tStart = t  # local t and not account for scr refresh
        end_msg.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(end_msg, 'tStartRefresh')  # time at next scr refresh
        # add timestamp to datafile
        thisExp.timestampOnFlip(win, 'end_msg.started')
        # update status
        end_msg.status = STARTED
        end_msg.setAutoDraw(True)
    
    # if end_msg is active this frame...
    if end_msg.status == STARTED:
        # update params
        pass
    
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
# the Routine "end" was not non-slip safe, so reset the non-slip timer
routineTimer.reset()

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
