#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This experiment was created using PsychoPy3 Experiment Builder (v2022.2.4),
    on November 21, 2023, at 11:49
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

# Run 'Before Experiment' code from codeBlank
jittered_duration_blank = 1.5
# Run 'Before Experiment' code from codeResponse
random_number = random()

if random_number > 0.5:
    same = 'ä'
    different = 'l'
else:
    same = 'l'
    different = 'ä'
# Run 'Before Experiment' code from code_end
jittered_duration_cross = 1
# Run 'Before Experiment' code from codeBlank
jittered_duration_blank = 1.5
# Run 'Before Experiment' code from codeResponse
random_number = random()

if random_number > 0.5:
    same = 'ä'
    different = 'l'
else:
    same = 'l'
    different = 'ä'
# Run 'Before Experiment' code from code_end
jittered_duration_cross = 1


# Ensure that relative paths start from the same directory as this script
_thisDir = os.path.dirname(os.path.abspath(__file__))
os.chdir(_thisDir)
# Store info about the experiment session
psychopyVersion = '2022.2.4'
expName = 'gca_discrimination'  # from the Builder filename that created this script
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
    originPath='C:\\Users\\sag22id\\Documents\\Projects\\GCA\\gaze_discrimination\\gca_discrimination.py',
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

# Setup iohub keyboard
ioConfig['Keyboard'] = dict(use_keymap='psychopy')

ioSession = '1'
if 'session' in expInfo:
    ioSession = str(expInfo['session'])
ioServer = io.launchHubServer(window=win, **ioConfig)
eyetracker = None

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

# --- Initialize components for Routine "startExp" ---
# Run 'Begin Experiment' code from codeStartExp
instruction_start = f"Bitte geben Sie an, ob die beiden Bilder, die Sie sehen gleich oder unterschiedlich sind.\nDrücken Sie die '{same}'-Taste, wenn die beiden Bilder gleich sind und die '{different}'-Taste, wenn die beiden Bilder unterschiedlich sind. \n\nFixieren Sie während des Experiments bitte stets das Fixationskreuz. \n\nDrücken Sie die Leertaste, um mit einer Übung zu starten."
textStartExp = visual.TextStim(win=win, name='textStartExp',
    text=instruction_start,
    font='Open Sans',
    pos=(0, 0), height=0.06, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-1.0);
spaceStartExp = keyboard.Keyboard()

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

# --- Initialize components for Routine "image1" ---
firstImage = visual.ImageStim(
    win=win,
    name='firstImage', 
    image='sin', mask=None, anchor='center',
    ori=0.0, pos=[0,0], size=None,
    color=[1,1,1], colorSpace='rgb', opacity=None,
    flipHoriz=False, flipVert=False,
    texRes=128.0, interpolate=True, depth=0.0)
fixcrossImage1 = visual.TextStim(win=win, name='fixcrossImage1',
    text='+',
    font='Open Sans',
    pos=(0, 0), height=0.1, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-1.0);

# --- Initialize components for Routine "interstimulus_cross" ---
fixcrossITI = visual.TextStim(win=win, name='fixcrossITI',
    text='+',
    font='Open Sans',
    pos=(0, 0), height=0.1, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-1.0);

# --- Initialize components for Routine "image2" ---
secondImage = visual.ImageStim(
    win=win,
    name='secondImage', 
    image='sin', mask=None, anchor='center',
    ori=0.0, pos=[0,0], size=None,
    color=[1,1,1], colorSpace='rgb', opacity=None,
    flipHoriz=False, flipVert=False,
    texRes=128.0, interpolate=True, depth=0.0)
fixcrossImage2 = visual.TextStim(win=win, name='fixcrossImage2',
    text='+',
    font='Open Sans',
    pos=(0, 0), height=0.1, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-1.0);

# --- Initialize components for Routine "response" ---
# Run 'Begin Experiment' code from codeResponse
print(f"Key Mappings: Same = {same}, Different = {different}")
sameKey = keyboard.Keyboard()
differentKey = keyboard.Keyboard()

# --- Initialize components for Routine "feedback" ---

# --- Initialize components for Routine "crossEnd" ---
fixcrossEnd = visual.TextStim(win=win, name='fixcrossEnd',
    text='+',
    font='Open Sans',
    pos=(0, 0), height=0.1, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-1.0);

# --- Initialize components for Routine "startTask" ---
# Run 'Begin Experiment' code from codeStartTask
instruction_task = f"Bitte geben Sie auch im nächsten Block an, ob die beiden Bilder, die Sie sehen gleich oder unterschiedlich sind.\nDrücken Sie die '{same}'-Taste, wenn die beiden Bilder gleich sind und die '{different}'-Taste, wenn die beiden Bilder unterschiedlich sind. \n\nFixieren Sie während des Experiments bitte stets das Fixationskreuz. \n\nDrücken Sie die Leertaste, um zu starten."
textStartTask = visual.TextStim(win=win, name='textStartTask',
    text=instruction_task,
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

# --- Initialize components for Routine "image1" ---
firstImage = visual.ImageStim(
    win=win,
    name='firstImage', 
    image='sin', mask=None, anchor='center',
    ori=0.0, pos=[0,0], size=None,
    color=[1,1,1], colorSpace='rgb', opacity=None,
    flipHoriz=False, flipVert=False,
    texRes=128.0, interpolate=True, depth=0.0)
fixcrossImage1 = visual.TextStim(win=win, name='fixcrossImage1',
    text='+',
    font='Open Sans',
    pos=(0, 0), height=0.1, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-1.0);

# --- Initialize components for Routine "interstimulus_cross" ---
fixcrossITI = visual.TextStim(win=win, name='fixcrossITI',
    text='+',
    font='Open Sans',
    pos=(0, 0), height=0.1, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-1.0);

# --- Initialize components for Routine "image2" ---
secondImage = visual.ImageStim(
    win=win,
    name='secondImage', 
    image='sin', mask=None, anchor='center',
    ori=0.0, pos=[0,0], size=None,
    color=[1,1,1], colorSpace='rgb', opacity=None,
    flipHoriz=False, flipVert=False,
    texRes=128.0, interpolate=True, depth=0.0)
fixcrossImage2 = visual.TextStim(win=win, name='fixcrossImage2',
    text='+',
    font='Open Sans',
    pos=(0, 0), height=0.1, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-1.0);

# --- Initialize components for Routine "response" ---
# Run 'Begin Experiment' code from codeResponse
print(f"Key Mappings: Same = {same}, Different = {different}")
sameKey = keyboard.Keyboard()
differentKey = keyboard.Keyboard()

# --- Initialize components for Routine "crossEnd" ---
fixcrossEnd = visual.TextStim(win=win, name='fixcrossEnd',
    text='+',
    font='Open Sans',
    pos=(0, 0), height=0.1, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-1.0);

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
    # Run 'Begin Routine' code from codeStartExp
    print("Start Experiment Screen. Press Space to continue.")
    
    
    spaceStartExp.keys = []
    spaceStartExp.rt = []
    _spaceStartExp_allKeys = []
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
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'blankScreen.started')
            blankScreen.setAutoDraw(True)
        if blankScreen.status == STARTED:
            # is it time to stop? (based on global clock, using actual start)
            if tThisFlipGlobal > blankScreen.tStartRefresh + jittered_duration_blank-frameTolerance:
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
    # Run 'End Routine' code from codeBlank
    # print(f"Start: {round(blankScreen.tStart, 3)}, End: {round(blankScreen.tStop, 3)}, Duration: {round(blankScreen.tStop - blankScreen.tStart, 3)}")
    # the Routine "blank" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset()
    
    # set up handler to look after randomisation of conditions etc
    learning_trials = data.TrialHandler(nReps=1.0, method='random', 
        extraInfo=expInfo, originPath=-1,
        trialList=data.importConditions(testFile, selection='0:10'),
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
        
        # --- Prepare to start Routine "image1" ---
        continueRoutine = True
        routineForceEnded = False
        # update component parameters for each repeat
        firstImage.setPos(position1)
        firstImage.setImage(eval(trialtype))
        # keep track of which components have finished
        image1Components = [firstImage, fixcrossImage1]
        for thisComponent in image1Components:
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
        
        # --- Run Routine "image1" ---
        while continueRoutine and routineTimer.getTime() < 0.15:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # *firstImage* updates
            if firstImage.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                # keep track of start time/frame for later
                firstImage.frameNStart = frameN  # exact frame index
                firstImage.tStart = t  # local t and not account for scr refresh
                firstImage.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(firstImage, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'firstImage.started')
                firstImage.setAutoDraw(True)
            if firstImage.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > firstImage.tStartRefresh + 0.15-frameTolerance:
                    # keep track of stop time/frame for later
                    firstImage.tStop = t  # not accounting for scr refresh
                    firstImage.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'firstImage.stopped')
                    firstImage.setAutoDraw(False)
            
            # *fixcrossImage1* updates
            if fixcrossImage1.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                fixcrossImage1.frameNStart = frameN  # exact frame index
                fixcrossImage1.tStart = t  # local t and not account for scr refresh
                fixcrossImage1.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(fixcrossImage1, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'fixcrossImage1.started')
                fixcrossImage1.setAutoDraw(True)
            if fixcrossImage1.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > fixcrossImage1.tStartRefresh + 0.15-frameTolerance:
                    # keep track of stop time/frame for later
                    fixcrossImage1.tStop = t  # not accounting for scr refresh
                    fixcrossImage1.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'fixcrossImage1.stopped')
                    fixcrossImage1.setAutoDraw(False)
            
            # check for quit (typically the Esc key)
            if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
                core.quit()
            
            # check if all components have finished
            if not continueRoutine:  # a component has requested a forced-end of Routine
                routineForceEnded = True
                break
            continueRoutine = False  # will revert to True if at least one component still running
            for thisComponent in image1Components:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # --- Ending Routine "image1" ---
        for thisComponent in image1Components:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        # using non-slip timing so subtract the expected duration of this Routine (unless ended on request)
        if routineForceEnded:
            routineTimer.reset()
        else:
            routineTimer.addTime(-0.150000)
        
        # --- Prepare to start Routine "interstimulus_cross" ---
        continueRoutine = True
        routineForceEnded = False
        # update component parameters for each repeat
        # keep track of which components have finished
        interstimulus_crossComponents = [fixcrossITI]
        for thisComponent in interstimulus_crossComponents:
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
        
        # --- Run Routine "interstimulus_cross" ---
        while continueRoutine and routineTimer.getTime() < 0.5:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # *fixcrossITI* updates
            if fixcrossITI.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                fixcrossITI.frameNStart = frameN  # exact frame index
                fixcrossITI.tStart = t  # local t and not account for scr refresh
                fixcrossITI.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(fixcrossITI, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'fixcrossITI.started')
                fixcrossITI.setAutoDraw(True)
            if fixcrossITI.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > fixcrossITI.tStartRefresh + 0.5-frameTolerance:
                    # keep track of stop time/frame for later
                    fixcrossITI.tStop = t  # not accounting for scr refresh
                    fixcrossITI.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'fixcrossITI.stopped')
                    fixcrossITI.setAutoDraw(False)
            
            # check for quit (typically the Esc key)
            if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
                core.quit()
            
            # check if all components have finished
            if not continueRoutine:  # a component has requested a forced-end of Routine
                routineForceEnded = True
                break
            continueRoutine = False  # will revert to True if at least one component still running
            for thisComponent in interstimulus_crossComponents:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # --- Ending Routine "interstimulus_cross" ---
        for thisComponent in interstimulus_crossComponents:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        # using non-slip timing so subtract the expected duration of this Routine (unless ended on request)
        if routineForceEnded:
            routineTimer.reset()
        else:
            routineTimer.addTime(-0.500000)
        
        # --- Prepare to start Routine "image2" ---
        continueRoutine = True
        routineForceEnded = False
        # update component parameters for each repeat
        secondImage.setPos(position2)
        secondImage.setImage(eval(comptype))
        # keep track of which components have finished
        image2Components = [secondImage, fixcrossImage2]
        for thisComponent in image2Components:
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
        
        # --- Run Routine "image2" ---
        while continueRoutine and routineTimer.getTime() < 0.15:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # *secondImage* updates
            if secondImage.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                # keep track of start time/frame for later
                secondImage.frameNStart = frameN  # exact frame index
                secondImage.tStart = t  # local t and not account for scr refresh
                secondImage.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(secondImage, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'secondImage.started')
                secondImage.setAutoDraw(True)
            if secondImage.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > secondImage.tStartRefresh + 0.15-frameTolerance:
                    # keep track of stop time/frame for later
                    secondImage.tStop = t  # not accounting for scr refresh
                    secondImage.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'secondImage.stopped')
                    secondImage.setAutoDraw(False)
            
            # *fixcrossImage2* updates
            if fixcrossImage2.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                fixcrossImage2.frameNStart = frameN  # exact frame index
                fixcrossImage2.tStart = t  # local t and not account for scr refresh
                fixcrossImage2.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(fixcrossImage2, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'fixcrossImage2.started')
                fixcrossImage2.setAutoDraw(True)
            if fixcrossImage2.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > fixcrossImage2.tStartRefresh + 0.15-frameTolerance:
                    # keep track of stop time/frame for later
                    fixcrossImage2.tStop = t  # not accounting for scr refresh
                    fixcrossImage2.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'fixcrossImage2.stopped')
                    fixcrossImage2.setAutoDraw(False)
            
            # check for quit (typically the Esc key)
            if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
                core.quit()
            
            # check if all components have finished
            if not continueRoutine:  # a component has requested a forced-end of Routine
                routineForceEnded = True
                break
            continueRoutine = False  # will revert to True if at least one component still running
            for thisComponent in image2Components:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # --- Ending Routine "image2" ---
        for thisComponent in image2Components:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        # using non-slip timing so subtract the expected duration of this Routine (unless ended on request)
        if routineForceEnded:
            routineTimer.reset()
        else:
            routineTimer.addTime(-0.150000)
        
        # --- Prepare to start Routine "response" ---
        continueRoutine = True
        routineForceEnded = False
        # update component parameters for each repeat
        sameKey.keys = []
        sameKey.rt = []
        _sameKey_allKeys = []
        differentKey.keys = []
        differentKey.rt = []
        _differentKey_allKeys = []
        # keep track of which components have finished
        responseComponents = [sameKey, differentKey]
        for thisComponent in responseComponents:
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
        
        # --- Run Routine "response" ---
        while continueRoutine and routineTimer.getTime() < 2.0:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            # Run 'Each Frame' code from codeResponse
            if sameKey.
            
            # *sameKey* updates
            waitOnFlip = False
            if sameKey.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                sameKey.frameNStart = frameN  # exact frame index
                sameKey.tStart = t  # local t and not account for scr refresh
                sameKey.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sameKey, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'sameKey.started')
                sameKey.status = STARTED
                # AllowedKeys looks like a variable named `same`
                if not type(same) in [list, tuple, np.ndarray]:
                    if not isinstance(same, str):
                        logging.error('AllowedKeys variable `same` is not string- or list-like.')
                        core.quit()
                    elif not ',' in same:
                        same = (same,)
                    else:
                        same = eval(same)
                # keyboard checking is just starting
                waitOnFlip = True
                win.callOnFlip(sameKey.clock.reset)  # t=0 on next screen flip
                win.callOnFlip(sameKey.clearEvents, eventType='keyboard')  # clear events on next screen flip
            if sameKey.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > sameKey.tStartRefresh + 2-frameTolerance:
                    # keep track of stop time/frame for later
                    sameKey.tStop = t  # not accounting for scr refresh
                    sameKey.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'sameKey.stopped')
                    sameKey.status = FINISHED
            if sameKey.status == STARTED and not waitOnFlip:
                theseKeys = sameKey.getKeys(keyList=list(same), waitRelease=False)
                _sameKey_allKeys.extend(theseKeys)
                if len(_sameKey_allKeys):
                    sameKey.keys = _sameKey_allKeys[-1].name  # just the last key pressed
                    sameKey.rt = _sameKey_allKeys[-1].rt
                    # a response ends the routine
                    continueRoutine = False
            
            # *differentKey* updates
            waitOnFlip = False
            if differentKey.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                differentKey.frameNStart = frameN  # exact frame index
                differentKey.tStart = t  # local t and not account for scr refresh
                differentKey.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(differentKey, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'differentKey.started')
                differentKey.status = STARTED
                # AllowedKeys looks like a variable named `different`
                if not type(different) in [list, tuple, np.ndarray]:
                    if not isinstance(different, str):
                        logging.error('AllowedKeys variable `different` is not string- or list-like.')
                        core.quit()
                    elif not ',' in different:
                        different = (different,)
                    else:
                        different = eval(different)
                # keyboard checking is just starting
                waitOnFlip = True
                win.callOnFlip(differentKey.clock.reset)  # t=0 on next screen flip
                win.callOnFlip(differentKey.clearEvents, eventType='keyboard')  # clear events on next screen flip
            if differentKey.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > differentKey.tStartRefresh + 2-frameTolerance:
                    # keep track of stop time/frame for later
                    differentKey.tStop = t  # not accounting for scr refresh
                    differentKey.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'differentKey.stopped')
                    differentKey.status = FINISHED
            if differentKey.status == STARTED and not waitOnFlip:
                theseKeys = differentKey.getKeys(keyList=list(different), waitRelease=False)
                _differentKey_allKeys.extend(theseKeys)
                if len(_differentKey_allKeys):
                    differentKey.keys = _differentKey_allKeys[-1].name  # just the last key pressed
                    differentKey.rt = _differentKey_allKeys[-1].rt
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
            for thisComponent in responseComponents:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # --- Ending Routine "response" ---
        for thisComponent in responseComponents:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        # Run 'End Routine' code from codeResponse
        if trialtype == comptype:
            condition = "same"
        else:
            condition = "different"
            
        
        # check responses
        if sameKey.keys in ['', [], None]:  # No response was made
            sameKey.keys = None
        thisExp.addData('sameKey.keys',sameKey.keys)
        if sameKey.keys != None:  # we had a response
            thisExp.addData('sameKey.rt', sameKey.rt)
        thisExp.nextEntry()
        # check responses
        if differentKey.keys in ['', [], None]:  # No response was made
            differentKey.keys = None
        thisExp.addData('differentKey.keys',differentKey.keys)
        if differentKey.keys != None:  # we had a response
            thisExp.addData('differentKey.rt', differentKey.rt)
        thisExp.nextEntry()
        # using non-slip timing so subtract the expected duration of this Routine (unless ended on request)
        if routineForceEnded:
            routineTimer.reset()
        else:
            routineTimer.addTime(-2.000000)
        
        # --- Prepare to start Routine "feedback" ---
        continueRoutine = True
        routineForceEnded = False
        # update component parameters for each repeat
        # keep track of which components have finished
        feedbackComponents = []
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
        while continueRoutine:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
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
        # the Routine "feedback" was not non-slip safe, so reset the non-slip timer
        routineTimer.reset()
        
        # --- Prepare to start Routine "crossEnd" ---
        continueRoutine = True
        routineForceEnded = False
        # update component parameters for each repeat
        # Run 'Begin Routine' code from code_end
        jittered_duration_cross = random()
        # keep track of which components have finished
        crossEndComponents = [fixcrossEnd]
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
            
            # *fixcrossEnd* updates
            if fixcrossEnd.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                fixcrossEnd.frameNStart = frameN  # exact frame index
                fixcrossEnd.tStart = t  # local t and not account for scr refresh
                fixcrossEnd.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(fixcrossEnd, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'fixcrossEnd.started')
                fixcrossEnd.setAutoDraw(True)
            if fixcrossEnd.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > fixcrossEnd.tStartRefresh + jittered_duration_cross-frameTolerance:
                    # keep track of stop time/frame for later
                    fixcrossEnd.tStop = t  # not accounting for scr refresh
                    fixcrossEnd.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'fixcrossEnd.stopped')
                    fixcrossEnd.setAutoDraw(False)
            
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
        
    # completed 1.0 repeats of 'learning_trials'
    
    
    # --- Prepare to start Routine "startTask" ---
    continueRoutine = True
    routineForceEnded = False
    # update component parameters for each repeat
    # Run 'Begin Routine' code from codeStartTask
    print("Start Experiment Screen. Press Space to continue.")
    
    
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
    thisExp.addData('spaceStartTask.keys',spaceStartTask.keys)
    if spaceStartTask.keys != None:  # we had a response
        thisExp.addData('spaceStartTask.rt', spaceStartTask.rt)
    thisExp.nextEntry()
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
            # add timestamp to datafile
            thisExp.timestampOnFlip(win, 'blankScreen.started')
            blankScreen.setAutoDraw(True)
        if blankScreen.status == STARTED:
            # is it time to stop? (based on global clock, using actual start)
            if tThisFlipGlobal > blankScreen.tStartRefresh + jittered_duration_blank-frameTolerance:
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
    # Run 'End Routine' code from codeBlank
    # print(f"Start: {round(blankScreen.tStart, 3)}, End: {round(blankScreen.tStop, 3)}, Duration: {round(blankScreen.tStop - blankScreen.tStart, 3)}")
    # the Routine "blank" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset()
    
    # set up handler to look after randomisation of conditions etc
    trials = data.TrialHandler(nReps=1.0, method='random', 
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
        
        # --- Prepare to start Routine "image1" ---
        continueRoutine = True
        routineForceEnded = False
        # update component parameters for each repeat
        firstImage.setPos(position1)
        firstImage.setImage(eval(trialtype))
        # keep track of which components have finished
        image1Components = [firstImage, fixcrossImage1]
        for thisComponent in image1Components:
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
        
        # --- Run Routine "image1" ---
        while continueRoutine and routineTimer.getTime() < 0.15:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # *firstImage* updates
            if firstImage.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                # keep track of start time/frame for later
                firstImage.frameNStart = frameN  # exact frame index
                firstImage.tStart = t  # local t and not account for scr refresh
                firstImage.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(firstImage, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'firstImage.started')
                firstImage.setAutoDraw(True)
            if firstImage.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > firstImage.tStartRefresh + 0.15-frameTolerance:
                    # keep track of stop time/frame for later
                    firstImage.tStop = t  # not accounting for scr refresh
                    firstImage.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'firstImage.stopped')
                    firstImage.setAutoDraw(False)
            
            # *fixcrossImage1* updates
            if fixcrossImage1.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                fixcrossImage1.frameNStart = frameN  # exact frame index
                fixcrossImage1.tStart = t  # local t and not account for scr refresh
                fixcrossImage1.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(fixcrossImage1, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'fixcrossImage1.started')
                fixcrossImage1.setAutoDraw(True)
            if fixcrossImage1.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > fixcrossImage1.tStartRefresh + 0.15-frameTolerance:
                    # keep track of stop time/frame for later
                    fixcrossImage1.tStop = t  # not accounting for scr refresh
                    fixcrossImage1.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'fixcrossImage1.stopped')
                    fixcrossImage1.setAutoDraw(False)
            
            # check for quit (typically the Esc key)
            if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
                core.quit()
            
            # check if all components have finished
            if not continueRoutine:  # a component has requested a forced-end of Routine
                routineForceEnded = True
                break
            continueRoutine = False  # will revert to True if at least one component still running
            for thisComponent in image1Components:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # --- Ending Routine "image1" ---
        for thisComponent in image1Components:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        # using non-slip timing so subtract the expected duration of this Routine (unless ended on request)
        if routineForceEnded:
            routineTimer.reset()
        else:
            routineTimer.addTime(-0.150000)
        
        # --- Prepare to start Routine "interstimulus_cross" ---
        continueRoutine = True
        routineForceEnded = False
        # update component parameters for each repeat
        # keep track of which components have finished
        interstimulus_crossComponents = [fixcrossITI]
        for thisComponent in interstimulus_crossComponents:
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
        
        # --- Run Routine "interstimulus_cross" ---
        while continueRoutine and routineTimer.getTime() < 0.5:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # *fixcrossITI* updates
            if fixcrossITI.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                fixcrossITI.frameNStart = frameN  # exact frame index
                fixcrossITI.tStart = t  # local t and not account for scr refresh
                fixcrossITI.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(fixcrossITI, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'fixcrossITI.started')
                fixcrossITI.setAutoDraw(True)
            if fixcrossITI.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > fixcrossITI.tStartRefresh + 0.5-frameTolerance:
                    # keep track of stop time/frame for later
                    fixcrossITI.tStop = t  # not accounting for scr refresh
                    fixcrossITI.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'fixcrossITI.stopped')
                    fixcrossITI.setAutoDraw(False)
            
            # check for quit (typically the Esc key)
            if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
                core.quit()
            
            # check if all components have finished
            if not continueRoutine:  # a component has requested a forced-end of Routine
                routineForceEnded = True
                break
            continueRoutine = False  # will revert to True if at least one component still running
            for thisComponent in interstimulus_crossComponents:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # --- Ending Routine "interstimulus_cross" ---
        for thisComponent in interstimulus_crossComponents:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        # using non-slip timing so subtract the expected duration of this Routine (unless ended on request)
        if routineForceEnded:
            routineTimer.reset()
        else:
            routineTimer.addTime(-0.500000)
        
        # --- Prepare to start Routine "image2" ---
        continueRoutine = True
        routineForceEnded = False
        # update component parameters for each repeat
        secondImage.setPos(position2)
        secondImage.setImage(eval(comptype))
        # keep track of which components have finished
        image2Components = [secondImage, fixcrossImage2]
        for thisComponent in image2Components:
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
        
        # --- Run Routine "image2" ---
        while continueRoutine and routineTimer.getTime() < 0.15:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # *secondImage* updates
            if secondImage.status == NOT_STARTED and tThisFlip >= 0-frameTolerance:
                # keep track of start time/frame for later
                secondImage.frameNStart = frameN  # exact frame index
                secondImage.tStart = t  # local t and not account for scr refresh
                secondImage.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(secondImage, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'secondImage.started')
                secondImage.setAutoDraw(True)
            if secondImage.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > secondImage.tStartRefresh + 0.15-frameTolerance:
                    # keep track of stop time/frame for later
                    secondImage.tStop = t  # not accounting for scr refresh
                    secondImage.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'secondImage.stopped')
                    secondImage.setAutoDraw(False)
            
            # *fixcrossImage2* updates
            if fixcrossImage2.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                fixcrossImage2.frameNStart = frameN  # exact frame index
                fixcrossImage2.tStart = t  # local t and not account for scr refresh
                fixcrossImage2.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(fixcrossImage2, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'fixcrossImage2.started')
                fixcrossImage2.setAutoDraw(True)
            if fixcrossImage2.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > fixcrossImage2.tStartRefresh + 0.15-frameTolerance:
                    # keep track of stop time/frame for later
                    fixcrossImage2.tStop = t  # not accounting for scr refresh
                    fixcrossImage2.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'fixcrossImage2.stopped')
                    fixcrossImage2.setAutoDraw(False)
            
            # check for quit (typically the Esc key)
            if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
                core.quit()
            
            # check if all components have finished
            if not continueRoutine:  # a component has requested a forced-end of Routine
                routineForceEnded = True
                break
            continueRoutine = False  # will revert to True if at least one component still running
            for thisComponent in image2Components:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # --- Ending Routine "image2" ---
        for thisComponent in image2Components:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        # using non-slip timing so subtract the expected duration of this Routine (unless ended on request)
        if routineForceEnded:
            routineTimer.reset()
        else:
            routineTimer.addTime(-0.150000)
        
        # --- Prepare to start Routine "response" ---
        continueRoutine = True
        routineForceEnded = False
        # update component parameters for each repeat
        sameKey.keys = []
        sameKey.rt = []
        _sameKey_allKeys = []
        differentKey.keys = []
        differentKey.rt = []
        _differentKey_allKeys = []
        # keep track of which components have finished
        responseComponents = [sameKey, differentKey]
        for thisComponent in responseComponents:
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
        
        # --- Run Routine "response" ---
        while continueRoutine and routineTimer.getTime() < 2.0:
            # get current time
            t = routineTimer.getTime()
            tThisFlip = win.getFutureFlipTime(clock=routineTimer)
            tThisFlipGlobal = win.getFutureFlipTime(clock=None)
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            # Run 'Each Frame' code from codeResponse
            if sameKey.
            
            # *sameKey* updates
            waitOnFlip = False
            if sameKey.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                sameKey.frameNStart = frameN  # exact frame index
                sameKey.tStart = t  # local t and not account for scr refresh
                sameKey.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(sameKey, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'sameKey.started')
                sameKey.status = STARTED
                # AllowedKeys looks like a variable named `same`
                if not type(same) in [list, tuple, np.ndarray]:
                    if not isinstance(same, str):
                        logging.error('AllowedKeys variable `same` is not string- or list-like.')
                        core.quit()
                    elif not ',' in same:
                        same = (same,)
                    else:
                        same = eval(same)
                # keyboard checking is just starting
                waitOnFlip = True
                win.callOnFlip(sameKey.clock.reset)  # t=0 on next screen flip
                win.callOnFlip(sameKey.clearEvents, eventType='keyboard')  # clear events on next screen flip
            if sameKey.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > sameKey.tStartRefresh + 2-frameTolerance:
                    # keep track of stop time/frame for later
                    sameKey.tStop = t  # not accounting for scr refresh
                    sameKey.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'sameKey.stopped')
                    sameKey.status = FINISHED
            if sameKey.status == STARTED and not waitOnFlip:
                theseKeys = sameKey.getKeys(keyList=list(same), waitRelease=False)
                _sameKey_allKeys.extend(theseKeys)
                if len(_sameKey_allKeys):
                    sameKey.keys = _sameKey_allKeys[-1].name  # just the last key pressed
                    sameKey.rt = _sameKey_allKeys[-1].rt
                    # a response ends the routine
                    continueRoutine = False
            
            # *differentKey* updates
            waitOnFlip = False
            if differentKey.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                differentKey.frameNStart = frameN  # exact frame index
                differentKey.tStart = t  # local t and not account for scr refresh
                differentKey.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(differentKey, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'differentKey.started')
                differentKey.status = STARTED
                # AllowedKeys looks like a variable named `different`
                if not type(different) in [list, tuple, np.ndarray]:
                    if not isinstance(different, str):
                        logging.error('AllowedKeys variable `different` is not string- or list-like.')
                        core.quit()
                    elif not ',' in different:
                        different = (different,)
                    else:
                        different = eval(different)
                # keyboard checking is just starting
                waitOnFlip = True
                win.callOnFlip(differentKey.clock.reset)  # t=0 on next screen flip
                win.callOnFlip(differentKey.clearEvents, eventType='keyboard')  # clear events on next screen flip
            if differentKey.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > differentKey.tStartRefresh + 2-frameTolerance:
                    # keep track of stop time/frame for later
                    differentKey.tStop = t  # not accounting for scr refresh
                    differentKey.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'differentKey.stopped')
                    differentKey.status = FINISHED
            if differentKey.status == STARTED and not waitOnFlip:
                theseKeys = differentKey.getKeys(keyList=list(different), waitRelease=False)
                _differentKey_allKeys.extend(theseKeys)
                if len(_differentKey_allKeys):
                    differentKey.keys = _differentKey_allKeys[-1].name  # just the last key pressed
                    differentKey.rt = _differentKey_allKeys[-1].rt
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
            for thisComponent in responseComponents:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # --- Ending Routine "response" ---
        for thisComponent in responseComponents:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)
        # Run 'End Routine' code from codeResponse
        if trialtype == comptype:
            condition = "same"
        else:
            condition = "different"
            
        
        # check responses
        if sameKey.keys in ['', [], None]:  # No response was made
            sameKey.keys = None
        thisExp.addData('sameKey.keys',sameKey.keys)
        if sameKey.keys != None:  # we had a response
            thisExp.addData('sameKey.rt', sameKey.rt)
        thisExp.nextEntry()
        # check responses
        if differentKey.keys in ['', [], None]:  # No response was made
            differentKey.keys = None
        thisExp.addData('differentKey.keys',differentKey.keys)
        if differentKey.keys != None:  # we had a response
            thisExp.addData('differentKey.rt', differentKey.rt)
        thisExp.nextEntry()
        # using non-slip timing so subtract the expected duration of this Routine (unless ended on request)
        if routineForceEnded:
            routineTimer.reset()
        else:
            routineTimer.addTime(-2.000000)
        
        # --- Prepare to start Routine "crossEnd" ---
        continueRoutine = True
        routineForceEnded = False
        # update component parameters for each repeat
        # Run 'Begin Routine' code from code_end
        jittered_duration_cross = random()
        # keep track of which components have finished
        crossEndComponents = [fixcrossEnd]
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
            
            # *fixcrossEnd* updates
            if fixcrossEnd.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
                # keep track of start time/frame for later
                fixcrossEnd.frameNStart = frameN  # exact frame index
                fixcrossEnd.tStart = t  # local t and not account for scr refresh
                fixcrossEnd.tStartRefresh = tThisFlipGlobal  # on global time
                win.timeOnFlip(fixcrossEnd, 'tStartRefresh')  # time at next scr refresh
                # add timestamp to datafile
                thisExp.timestampOnFlip(win, 'fixcrossEnd.started')
                fixcrossEnd.setAutoDraw(True)
            if fixcrossEnd.status == STARTED:
                # is it time to stop? (based on global clock, using actual start)
                if tThisFlipGlobal > fixcrossEnd.tStartRefresh + jittered_duration_cross-frameTolerance:
                    # keep track of stop time/frame for later
                    fixcrossEnd.tStop = t  # not accounting for scr refresh
                    fixcrossEnd.frameNStop = frameN  # exact frame index
                    # add timestamp to datafile
                    thisExp.timestampOnFlip(win, 'fixcrossEnd.stopped')
                    fixcrossEnd.setAutoDraw(False)
            
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
        
    # completed 1.0 repeats of 'trials'
    
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
